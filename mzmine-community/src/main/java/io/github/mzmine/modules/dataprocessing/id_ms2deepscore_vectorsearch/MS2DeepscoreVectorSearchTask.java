/*
 * Copyright (c) 2004-2025 The mzmine Development Team
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package io.github.mzmine.modules.dataprocessing.id_ms2deepscore_vectorsearch;

import ai.djl.MalformedModelException;
import ai.djl.repository.zoo.ModelNotFoundException;
import ai.djl.translate.TranslateException;
import io.github.mzmine.datamodel.MZmineProject;
import io.github.mzmine.datamodel.PolarityType;
import io.github.mzmine.datamodel.Scan;
import io.github.mzmine.datamodel.features.FeatureList;
import io.github.mzmine.datamodel.features.FeatureListRow;
import io.github.mzmine.modules.MZmineModule;
import io.github.mzmine.parameters.ParameterSet;
import io.github.mzmine.parameters.parametertypes.tolerances.MZTolerance;
import io.github.mzmine.taskcontrol.AbstractFeatureListTask;
import io.github.mzmine.taskcontrol.TaskStatus;
import io.github.mzmine.util.MemoryMapStorage;
import io.github.mzmine.util.scans.FragmentScanSelection;
import io.github.mzmine.util.scans.ScanUtils;
import io.github.mzmine.util.scans.similarity.SpectralSimilarity;
import io.github.mzmine.util.scans.similarity.impl.ms2deepscore.MS2DeepscoreModel;
import io.github.mzmine.datamodel.MassSpectrum;
import io.github.mzmine.datamodel.features.types.DataTypes;
import io.github.mzmine.datamodel.features.types.annotations.MS2DeepscoreMatchesType;
import io.github.mzmine.util.spectraldb.entry.DBEntryField;
import io.github.mzmine.util.spectraldb.entry.SpectralDBAnnotation;
import io.github.mzmine.util.spectraldb.entry.SpectralLibraryEntry;
import io.github.mzmine.util.spectraldb.entry.SpectralLibraryEntryFactory;
import io.github.mzmine.util.spectraldb.entry.SpectralLibrary;
import io.github.mzmine.project.ProjectService;
import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;

/**
 * Embeds selected MS2 scans using MS2Deepscore and queries a vector database for nearest neighbors.
 * Matches are attached as SpectralLibraryMatches to rows.
 */
public class MS2DeepscoreVectorSearchTask extends AbstractFeatureListTask {

  private static final Logger logger = Logger.getLogger(MS2DeepscoreVectorSearchTask.class.getName());

  private final @NotNull FeatureList[] featureLists;
  private final FragmentScanSelection scanSelector;
  private final File modelFile;
  private final File settingsFile;
  private final String dbUri;
  private final String collectionPos;
  private final String collectionNeg;
  private final int topK;
  private final double minScore;
  private final MZTolerance precursorTol;
  // keep parameter for future optional cosine confirmation refinement
  private final boolean confirmWithCosine;
  private int processedItems;

  public MS2DeepscoreVectorSearchTask(MZmineProject project, @NotNull FeatureList[] featureLists,
      ParameterSet parameters, @Nullable MemoryMapStorage storage, @NotNull Instant moduleCallDate,
      @NotNull Class<? extends MZmineModule> moduleClass) {
    super(storage, moduleCallDate, parameters, moduleClass);
    this.featureLists = featureLists;
    this.scanSelector = parameters.getParameter(MS2DeepscoreVectorSearchParameters.SPECTRA_SELECT)
        .createFragmentScanSelection(getMemoryMapStorage());
    this.modelFile = parameters.getValue(MS2DeepscoreVectorSearchParameters.MODEL_FILE);
    this.settingsFile = parameters.getValue(MS2DeepscoreVectorSearchParameters.SETTINGS_FILE);
    this.dbUri = parameters.getValue(MS2DeepscoreVectorSearchParameters.DB_URI);
    this.collectionPos = parameters.getValue(MS2DeepscoreVectorSearchParameters.COLLECTION_POS);
    this.collectionNeg = parameters.getValue(MS2DeepscoreVectorSearchParameters.COLLECTION_NEG);
    this.topK = parameters.getValue(MS2DeepscoreVectorSearchParameters.TOP_K);
    this.minScore = parameters.getValue(MS2DeepscoreVectorSearchParameters.MIN_SCORE);
    this.precursorTol = parameters.getValue(MS2DeepscoreVectorSearchParameters.PRECURSOR_TOL);
    this.confirmWithCosine = parameters.getValue(MS2DeepscoreVectorSearchParameters.CONFIRM_WITH_COSINE);
    totalItems = 0;
    processedItems = 0;
    for (FeatureList fl : featureLists) {
      // ensure separate MS2Deepscore column group is visible in the table
      fl.addRowType(DataTypes.get(MS2DeepscoreMatchesType.class));
      totalItems += fl.getNumberOfRows();
    }
  }

  @Override
  protected void process() {
    setStatus(TaskStatus.PROCESSING);
    try {
      final File[] modelFiles = resolveModelFiles();
      try (var model = new MS2DeepscoreModel(modelFiles[0], modelFiles[1]);
           var client = createClient()) {
        client.connect(dbUri);
        for (FeatureList fl : featureLists) {
          processFeatureList(fl, model, client);
        }
        setStatus(TaskStatus.FINISHED);
      }
    } catch (ModelNotFoundException | MalformedModelException | IOException e) {
      error("Model/DB error: " + e.getMessage(), e);
    }
  }

  private void processFeatureList(FeatureList featureList, MS2DeepscoreModel model,
      VectorDbClient client) {
    // Build an index of imported spectral libraries by DB# for enrichment
    final Map<String, SpectralLibraryEntry> libIndexByDbId = buildLibraryIndex();
    for (FeatureListRow row : featureList.getRows()) {
      if (isCanceled()) {
        return;
      }
      try {
        final Scan scan = pickOneScan(row);
        // must have a mass list with data points
        if (scan == null || scan.getMassList() == null
            || scan.getMassList().getNumberOfDataPoints() == 0) {
          processedItems++;
          continue;
        }
        // Predict embedding using original scan to preserve precursor m/z and polarity
        var embeddingNd = model.predictEmbedding(List.of(scan));
        float[] embedding = embeddingNd.toFloatArray();

        // Choose collection by ion mode; if unknown, query both collections and merge
        final PolarityType pol = ScanUtils.getPolarity(scan);
        final boolean queryBoth = pol == null || pol == PolarityType.UNKNOWN;

        // Optional precursor m/z scalar prefilter; if no precursor available, query by embedding only
        Map<String, Object> filters = null;
        final Double precursorMz = scan.getPrecursorMz();
        if (precursorMz != null) {
          double lo = precursorTol.getToleranceRange(precursorMz).lowerEndpoint();
          double hi = precursorTol.getToleranceRange(precursorMz).upperEndpoint();
          filters = new HashMap<>();
          filters.put("precursor_mz_between", new double[]{lo, hi});
        }

        // Query ANN
        List<VectorDbHit> hits;
        if (queryBoth) {
          List<VectorDbHit> posHits = client.query(collectionPos, embedding, topK, minScore, filters);
          List<VectorDbHit> negHits = client.query(collectionNeg, embedding, topK, minScore, filters);
          hits = new ArrayList<>(posHits.size() + negHits.size());
          hits.addAll(posHits);
          hits.addAll(negHits);
          hits.sort((a,b) -> Double.compare(b.score(), a.score()));
          if (hits.size() > topK) {
            hits = hits.subList(0, topK);
          }
        } else {
          final String collection = (pol == PolarityType.NEGATIVE) ? collectionNeg : collectionPos;
          hits = client.query(collection, embedding, topK, minScore, filters);
        }
        if (hits.isEmpty()) {
          processedItems++;
          continue;
        }

        // Convert hits to SpectralDBAnnotation
        List<SpectralDBAnnotation> annotations = new ArrayList<>(hits.size());
        for (VectorDbHit hit : hits) {
          SpectralLibraryEntry entry = hitToEntry(hit, libIndexByDbId);
          SpectralSimilarity sim = new SpectralSimilarity("MS2DeepscoreANN", hit.score(), 0,
              0d);
          SpectralDBAnnotation ann = new SpectralDBAnnotation(entry, sim, scan, null,
              scan.getPrecursorMz(), row.getAverageRT());
          annotations.add(ann);
        }

        // Attach to row under MS2Deepscore-specific column group for side-by-side comparison
        // with classical spectral search results.
        var existing = row.get(MS2DeepscoreMatchesType.class);
        List<SpectralDBAnnotation> ms2ds = existing == null ? new ArrayList<>() : new ArrayList<>(existing);
        ms2ds.addAll(annotations);
        row.set(MS2DeepscoreMatchesType.class, ms2ds);
      } catch (TranslateException e) {
        // Likely due to missing precursor m/z after internal filtering; ignore this row
        logger.log(Level.FINE, "Embedding skipped for row " + row.getID() + ": " + e.getMessage());
      } catch (IOException e) {
        logger.log(Level.WARNING, "Vector DB query failed for row " + row.getID(), e);
      } finally {
        processedItems++;
      }
    }
  }

  private @Nullable Scan pickOneScan(@NotNull FeatureListRow row) {
    List<Scan> scans = scanSelector.getAllFragmentSpectra(row);
    return scans.isEmpty() ? null : scans.getFirst();
  }

  private SpectralLibraryEntry hitToEntry(VectorDbHit hit, Map<String, SpectralLibraryEntry> libIndexByDbId) {
    // Prefer enriching from loaded MSP libraries using DB# (db_id) or ext_id
    final String dbId = extractDbId(hit);
    if (dbId != null && libIndexByDbId != null) {
      final SpectralLibraryEntry libEntry = libIndexByDbId.get(dbId);
      if (libEntry != null) {
        return libEntry;
      }
    }

    // Fallback: Build an entry with available metadata and optional embedded peaks
    Map<DBEntryField, Object> fields = new HashMap<>();
    fields.put(DBEntryField.ENTRY_ID, dbId != null ? dbId : hit.id());
    if (hit.precursorMz() != null) {
      fields.put(DBEntryField.PRECURSOR_MZ, hit.precursorMz());
    }
    io.github.mzmine.datamodel.DataPoint[] dataPoints = new io.github.mzmine.datamodel.DataPoint[]{};
    if (hit.metadata() != null) {
      putIfPresent(fields, DBEntryField.NAME, hit.metadata().get("name"));
      putIfPresent(fields, DBEntryField.FORMULA, hit.metadata().get("formula"));
      putIfPresent(fields, DBEntryField.SMILES, hit.metadata().get("smiles"));
      putIfPresent(fields, DBEntryField.INCHI, hit.metadata().get("inchi"));
      putIfPresent(fields, DBEntryField.INCHIKEY, hit.metadata().get("inchikey"));
      putIfPresent(fields, DBEntryField.POLARITY, hit.metadata().get("ionmode"));
      putIfPresent(fields, DBEntryField.COLLISION_ENERGY, hit.metadata().get("collision_energy"));
      putIfPresent(fields, DBEntryField.COMMENT, hit.metadata().get("source"));

      // Optional embedded peaks: array[[mz, intensity]]
      Object peaks = hit.metadata().get("peaks");
      if (peaks instanceof java.util.List<?> list && !list.isEmpty()) {
        java.util.ArrayList<io.github.mzmine.datamodel.DataPoint> dps = new java.util.ArrayList<>(list.size());
        for (Object o : list) {
          if (o instanceof java.util.List<?> pair && pair.size() >= 2) {
            try {
              double mz = Double.parseDouble(String.valueOf(pair.get(0)));
              double inten = Double.parseDouble(String.valueOf(pair.get(1)));
              dps.add(new io.github.mzmine.datamodel.impl.SimpleDataPoint(mz, inten));
            } catch (Exception ignored) {}
          }
        }
        if (!dps.isEmpty()) {
          dataPoints = dps.toArray(new io.github.mzmine.datamodel.DataPoint[0]);
        }
      }
    }
    return SpectralLibraryEntryFactory.create(getMemoryMapStorage(), fields, dataPoints);
  }

  private String extractDbId(VectorDbHit hit) {
    if (hit.metadata() != null) {
      Object db = hit.metadata().get("db_id");
      if (db != null) return String.valueOf(db);
      Object acc = hit.metadata().get("accession");
      if (acc != null) return String.valueOf(acc);
      Object ext = hit.metadata().get("ext_id");
      if (ext != null) return String.valueOf(ext);
    }
    return hit.id();
  }

  private Map<String, SpectralLibraryEntry> buildLibraryIndex() {
    try {
      final var project = ProjectService.getProjectManager().getCurrentProject();
      if (project == null) return Map.of();
      final List<SpectralLibrary> libs = project.getCurrentSpectralLibraries();
      if (libs == null || libs.isEmpty()) return Map.of();
      final Map<String, SpectralLibraryEntry> index = new HashMap<>();
      for (SpectralLibrary lib : libs) {
        for (SpectralLibraryEntry e : lib.getEntries()) {
          final Object id = e.getField(DBEntryField.ENTRY_ID).orElse(null);
          if (id != null) index.put(String.valueOf(id), e);
        }
      }
      return index;
    } catch (Exception ignored) {
      return Map.of();
    }
  }

  private static void putIfPresent(Map<DBEntryField, Object> fields, DBEntryField key,
      Object value) {
    if (value != null) {
      fields.put(key, value);
    }
  }

  private VectorDbClient createClient() {
    // Choose client strictly by URI scheme; default to Milvus if looks like milvus or http(s)
    if (dbUri != null && (dbUri.startsWith("milvus://") || dbUri.startsWith("http://") || dbUri.startsWith("https://"))) {
      return new MilvusVectorDbClient();
    }
    // Default to Milvus as requested; no SQL fallback
    return new MilvusVectorDbClient();
  }

  /**
   * Enforce unit_max intensity normalization on the mass list for DB-compatible embeddings.
   */
  private MassSpectrum normalizeSpectrumUnitMax(Scan scan) {
    MassSpectrum ms = scan.getMassList() != null ? scan.getMassList() : scan;
    int n = ms.getNumberOfDataPoints();
    double[] mz = new double[n];
    double[] intens = new double[n];
    double max = 0;
    for (int i = 0; i < n; i++) {
      mz[i] = ms.getMzValue(i);
      double v = ms.getIntensityValue(i);
      intens[i] = v;
      if (v > max) max = v;
    }
    if (max > 0) {
      for (int i = 0; i < n; i++) intens[i] = intens[i] / max;
    }
    return new io.github.mzmine.datamodel.impl.SimpleMassSpectrum(mz, intens);
  }

  @Override
  public String getTaskDescription() {
    return "MS2Deepscore vector search";
  }

  @Override
  protected @NotNull List<FeatureList> getProcessedFeatureLists() {
    return List.of(featureLists);
  }

  @Override
  public double getFinishedPercentage() {
    if (totalItems <= 0) {
      return 0;
    }
    return processedItems / (double) totalItems;
  }

  /**
   * Ensure model/settings exist. If missing, auto-download the dual-mode MS2Deepscore from Zenodo
   * (https://zenodo.org/records/16962640) into ~/.mzmine/models/ms2ds_dual/ and use those.
   */
  private File[] resolveModelFiles() throws IOException {
    final boolean haveLocal = modelFile != null && modelFile.exists() && settingsFile != null
        && settingsFile.exists();
    if (haveLocal) {
      return new File[]{modelFile, settingsFile};
    }

    final File targetDir = new File(System.getProperty("user.home"), 
        ".mzmine/models/ms2ds_dual");
    if (!targetDir.exists() && !targetDir.mkdirs()) {
      throw new IOException("Cannot create model directory " + targetDir);
    }
    final File targetModel = new File(targetDir, "ms2deepscore_model.pt");
    final File targetSettings = new File(targetDir, "settings.json");

    if (!targetModel.exists()) {
      downloadFile("https://zenodo.org/records/16962640/files/ms2deepscore_model.pt?download=1",
          targetModel);
    }
    if (!targetSettings.exists()) {
      downloadFile("https://zenodo.org/records/16962640/files/settings.json?download=1",
          targetSettings);
    }
    return new File[]{targetModel, targetSettings};
  }

  private void downloadFile(String url, File out) throws IOException {
    HttpClient client = HttpClient.newHttpClient();
    HttpRequest request = HttpRequest.newBuilder().uri(URI.create(url)).GET().build();
    try {
      HttpResponse<java.io.InputStream> resp = client.send(request,
          HttpResponse.BodyHandlers.ofInputStream());
      if (resp.statusCode() < 200 || resp.statusCode() >= 300) {
        throw new IOException("Download failed: " + resp.statusCode() + " for " + url);
      }
      try (java.io.InputStream in = resp.body();
           java.io.OutputStream os = new java.io.FileOutputStream(out)) {
        in.transferTo(os);
      }
    } catch (InterruptedException e) {
      Thread.currentThread().interrupt();
      throw new IOException("Interrupted while downloading " + url, e);
    }
  }
}


