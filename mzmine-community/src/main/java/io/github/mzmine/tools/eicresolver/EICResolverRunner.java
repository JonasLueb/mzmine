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

package io.github.mzmine.tools.eicresolver;

import com.google.common.collect.Range;
import io.github.mzmine.datamodel.features.ModularFeatureList;
import io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.GeneralResolverParameters;
import io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.Resolver;
import io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.ResolvingDimension;
import io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.minimumsearch.MinimumSearchFeatureResolverParameters;
import io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.savitzkygolay.SavitzkyGolayFeatureResolverParameters;
import io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.CwtRidgeResolverParameters;
import io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.WaveletResolverParameters;
import io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.AdvancedWaveletParameters;
import io.github.mzmine.project.impl.RawDataFileImpl;
import io.github.mzmine.util.CSVParsingUtils;
import io.github.mzmine.util.MemoryMapStorage;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.geom.Rectangle2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Stream;
import java.util.stream.Collectors;
import org.jetbrains.annotations.NotNull;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYShapeAnnotation;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.IntervalMarker;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.ui.Layer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import io.github.mzmine.datamodel.impl.SimpleScan;
import io.github.mzmine.datamodel.PolarityType;
import io.github.mzmine.datamodel.MassSpectrumType;
import io.github.mzmine.datamodel.featuredata.impl.SimpleIonTimeSeries;

/**
 * Standalone runner to evaluate chromatogram deconvolution resolvers on TSV EICs with labeled
 * ranges. It scans a directory for pairs of files "*_full.tsv" and "*_ranges.tsv".
 *
 * Usage:
 *   -DinputDir=/absolute/path -DoutDir=/absolute/output/path -DminPoints=3 -Diou=0.3
 * Or call main with args: inputDir [outDir]
 */
public class EICResolverRunner {

  public static void main(String[] args) throws Exception {
    final String inputDirArg = (args != null && args.length > 0 && args[0] != null && !args[0].isBlank())
        ? args[0]
        : System.getProperty("inputDir", "");
    if (inputDirArg == null || inputDirArg.isBlank()) {
      System.err.println("Please provide input directory via -DinputDir or first CLI arg.");
      return;
    }
    final Path inputDir = Paths.get(inputDirArg).toAbsolutePath().normalize();
    final String outDirArg = (args != null && args.length > 1 && args[1] != null && !args[1].isBlank())
        ? args[1]
        : System.getProperty("outDir", inputDir.resolve("resolver_eval").toString());
    final Path outDir = Paths.get(outDirArg).toAbsolutePath().normalize();
    final double iouThreshold = Double.parseDouble(System.getProperty("iou", "0.30"));
    final int minPoints = Integer.parseInt(System.getProperty("minPoints", "3"));
    final String format = System.getProperty("format", "png");
    final boolean outPng = format.equalsIgnoreCase("png") || format.equalsIgnoreCase("both");
    final boolean outHtml = format.equalsIgnoreCase("html") || format.equalsIgnoreCase("both");
    final double rtTol = Double.parseDouble(System.getProperty("rtTol", "0.01"));
    final boolean DEBUG = Boolean.parseBoolean(System.getProperty("debug", "false"));
    final boolean cwtAuto = Boolean.parseBoolean(System.getProperty("cwtAuto", "false"));
    final String resolverFilter = System.getProperty("resolver", "").trim();
    final String resolversFilter = System.getProperty("resolvers", "").trim();
    final String configId = System.getProperty("configId", "").trim();

    Files.createDirectories(outDir);
    System.out.printf(Locale.US, "Input: %s%nOutput: %s%nIoU>=%.2f, minPoints=%d, format=%s%n", inputDir, outDir,
        iouThreshold, minPoints, format);

    final List<Path> eicFiles;
    try (Stream<Path> s = Files.list(inputDir)) {
      eicFiles = s.filter(p -> p.getFileName().toString().endsWith("_full.tsv")).sorted().toList();
    }
    if (eicFiles.isEmpty()) {
      System.err.println("No *_full.tsv files found in " + inputDir);
      return;
    }

    // Create minimal environment for resolvers
    final RawDataFileImpl raw = new RawDataFileImpl("EIC_TS", inputDir.toString(), (MemoryMapStorage) null);
    final ModularFeatureList flist = new ModularFeatureList("EIC_FL", null, raw);

    // Prepare resolver parameter sets with sensible defaults
    List<ResolverSpec> resolvers = buildResolvers(flist, minPoints);
    if (!resolverFilter.isEmpty()) {
      resolvers = resolvers.stream().filter(r -> r.name.equalsIgnoreCase(resolverFilter)).toList();
    } else if (!resolversFilter.isEmpty()) {
      final List<String> allow = Arrays.stream(resolversFilter.split(",")).map(String::trim).filter(s -> !s.isEmpty()).toList();
      resolvers = resolvers.stream().filter(r -> allow.stream().anyMatch(a -> a.equalsIgnoreCase(r.name))).toList();
    }
    if (resolvers.isEmpty()) {
      System.err.println("No resolvers selected after filter. Use -Dresolver=MinimumSearch|SavitzkyGolay|CWT_Ridge or -Dresolvers=...");
      return;
    }

    final Path csvReport = outDir.resolve("evaluation_summary.csv");
    final StringBuilder report = new StringBuilder();
    report.append("eic\tresolver\ttp\tfp\tfn\tprecision\trecall\tf1\n");

    // aggregate per resolver: [tp, fp, fn]
    final Map<String, long[]> totals = new LinkedHashMap<>();
    resolvers.forEach(r -> totals.put(r.name, new long[3]));

    int savedPng = 0, savedHtml = 0;
    final List<String> eicNames = new ArrayList<>();
    for (Path fullTsv : eicFiles) {
      final String stem = fullTsv.getFileName().toString().replace("_full.tsv", "");
      final Path rangesTsv = inputDir.resolve(stem + "_ranges.tsv");
      // If ranges file is missing, treat as no features (empty labels)

      TsvEic eic = readEic(fullTsv.toFile());
      if (eic.x.length < 3) {
        System.err.println("Skipping small EIC (" + eic.x.length + " pts): " + stem);
        continue;
      }
      final List<Range<Double>> labels = Files.exists(rangesTsv) ? readRanges(rangesTsv.toFile()) : List.of();

      // Ensure scans monotonic in RT and construct dummy scans (not used for resolve(x,y), but kept for completeness)
      eic = ensureMonotonicRt(eic);

      // Debug: compute detrended/prezeroed and log stats
      final double[] detr = preprocessIntensities(eic.y);
      if (DEBUG) {
        // simple inline stats
        int zeros=0; double min=Double.POSITIVE_INFINITY, max=Double.NEGATIVE_INFINITY;
        for (double v: detr){ if(v==0) zeros++; min=Math.min(min,v); max=Math.max(max,v);}        
        System.out.printf(Locale.US, "DEBUG EIC=%s n=%d detr[min,max]=[%.3g,%.3g] zeros=%d\n", stem, eic.x.length, min, max, zeros);
        // optional debug TSV
        final Path dbg = outDir.resolve(stem + "_debug.tsv");
        try {
          final StringBuilder sb = new StringBuilder("rt\traw\tdetr\n");
          for (int i = 0; i < eic.x.length; i++) {
            sb.append(eic.x[i]).append('\t').append(eic.y[i]).append('\t').append(detr[i]).append('\n');
          }
          Files.writeString(dbg, sb.toString());
        } catch (Exception ignore) {}
      }

      final Path pngOut = outDir.resolve(stem + ".png");
      final Path htmlOut = outDir.resolve(stem + ".html");
      final Map<String, List<Range<Double>>> resultsByResolver = new LinkedHashMap<>();

      for (ResolverSpec spec : resolvers) {
        final List<Range<Double>> predicted = predictRanges(spec, eic, raw);
        resultsByResolver.put(spec.name, predicted);
        final Metrics m = evaluate(labels, predicted, iouThreshold, rtTol);
        report.append(stem).append('\t').append(spec.name).append('\t')
            .append(m.tp).append('\t').append(m.fp).append('\t').append(m.fn).append('\t')
            .append(format(m.precision)).append('\t').append(format(m.recall)).append('\t')
            .append(format(m.f1)).append('\n');
        final long[] agg = totals.get(spec.name);
        agg[0] += m.tp; agg[1] += m.fp; agg[2] += m.fn;
        if (DEBUG) {
          System.out.printf(Locale.US, "DEBUG EIC=%s resolver=%s predicted=%d\n", stem, spec.name, predicted.size());
        }
      }

      if (outPng) {
        saveChart(pngOut.toFile(), stem, eic, labels, resultsByResolver);
        savedPng++;
      }
      if (outHtml) {
        saveHtml(htmlOut.toFile(), stem, eic, labels, resultsByResolver);
        savedHtml++;
      }
      eicNames.add(stem);
    }

    Files.writeString(csvReport, report.toString());
    // build aggregate table html
    final StringBuilder aggTable = new StringBuilder();
    aggTable.append("<table border=\"1\" cellspacing=\"0\" cellpadding=\"4\"><thead><tr><th>Resolver</th><th>TP</th><th>FP</th><th>FN</th><th>Precision</th><th>Recall</th><th>F1</th></tr></thead><tbody>");
    for (Map.Entry<String,long[]> e : totals.entrySet()) {
      final long tp = e.getValue()[0], fp = e.getValue()[1], fn = e.getValue()[2];
      final double p = (tp+fp)==0?0:tp/(double)(tp+fp);
      final double r = (tp+fn)==0?0:tp/(double)(tp+fn);
      final double f1 = (p+r)==0?0:2*p*r/(p+r);
      aggTable.append("<tr><td>").append(e.getKey()).append("</td><td>")
          .append(tp).append("</td><td>").append(fp).append("</td><td>")
          .append(fn).append("</td><td>").append(format(p)).append("</td><td>")
          .append(format(r)).append("</td><td>").append(format(f1)).append("</td></tr>");
    }
    aggTable.append("</tbody></table>");

    // index html
    if (outHtml) {
      final Path index = outDir.resolve("index.html");
      final StringBuilder html = new StringBuilder();
      html.append("<!doctype html><html><head><meta charset=\"utf-8\"><title>EIC Resolver Results</title>")
          .append("<style>body{font-family:Sans-Serif} #viewer{width:100%;height:720px;border:1px solid #ccc} select{min-width:400px} pre{background:#f7f7f7;padding:8px;border:1px solid #ddd}</style>")
          .append("</head><body>");
      html.append("<h2>EIC Resolver Results</h2>");
      if (!configId.isEmpty()) {
        html.append("<p><b>Config:</b> ").append(configId).append("</p>");
      }
      html.append("<p>IoU=").append(iouThreshold).append(", rtTol=").append(rtTol).append(", minPoints=").append(minPoints).append("</p>");
      html.append("<p>CSV: <code>").append(csvReport.getFileName()).append("</code></p>");
      html.append("<h3>Aggregate metrics</h3>").append(aggTable.toString());
      // show resolver parameters
      html.append("<h3>Resolver parameters</h3><pre>").append(parameterSummary(minPoints, iouThreshold, rtTol)).append("</pre>");
      html.append("<label for=\"eic\">Select EIC: </label><select id=\"eic\">");
      for (String n : eicNames) html.append("<option value=\"").append(n).append(".html\">").append(n).append("</option>");
      html.append("</select>");
      html.append("<iframe id=\"viewer\"></iframe>\n<script>\nconst sel=document.getElementById('eic');\nfunction load(){const v=sel.value; document.getElementById('viewer').src=encodeURI(v);}\nsel.addEventListener('change',load);\nload();\n</script>");
      html.append("</body></html>");
      Files.writeString(index, html.toString());
    }

    // parameter summary
    final Path params = outDir.resolve("resolver_parameters.txt");
    Files.writeString(params, parameterSummary(minPoints, iouThreshold, rtTol));

    System.out.printf(Locale.US, "Saved %d PNG%s, %d HTML%s, report: %s, params: %s%n",
        savedPng, savedPng == 1 ? "" : "s", savedHtml, savedHtml == 1 ? "" : "s", csvReport.getFileName(), params.getFileName());
  }

  private static @NotNull List<ResolverSpec> buildResolvers(ModularFeatureList flist, int minPoints) {
    final List<ResolverSpec> list = new ArrayList<>();

    // Minimum search
    final MinimumSearchFeatureResolverParameters minParams = new MinimumSearchFeatureResolverParameters();
    minParams.setParameter(GeneralResolverParameters.dimension, ResolvingDimension.RETENTION_TIME);
    minParams.setParameter(GeneralResolverParameters.MIN_NUMBER_OF_DATAPOINTS, minPoints);
    // Permissive defaults with -D overrides
    final double msMinAbs = Double.parseDouble(System.getProperty("minAbs", "20"));
    final double msMinRel = Double.parseDouble(System.getProperty("minRel", "0.01"));
    final double msRatio = Double.parseDouble(System.getProperty("ratio", "1.3"));
    final double msThresh = Double.parseDouble(System.getProperty("msThresh", "0.3"));
    final double msSearch = Double.parseDouble(System.getProperty("searchRange", "0.03"));
    final double msDurMin = Double.parseDouble(System.getProperty("msDurMin", "0.0"));
    final double msDurMax = Double.parseDouble(System.getProperty("msDurMax", "10.0"));
    minParams.setParameter(io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.minimumsearch.MinimumSearchFeatureResolverParameters.MIN_ABSOLUTE_HEIGHT, msMinAbs);
    minParams.setParameter(io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.minimumsearch.MinimumSearchFeatureResolverParameters.MIN_RELATIVE_HEIGHT, msMinRel);
    minParams.setParameter(io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.minimumsearch.MinimumSearchFeatureResolverParameters.MIN_RATIO, msRatio);
    minParams.setParameter(io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.minimumsearch.MinimumSearchFeatureResolverParameters.CHROMATOGRAPHIC_THRESHOLD_LEVEL, msThresh);
    minParams.setParameter(io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.minimumsearch.MinimumSearchFeatureResolverParameters.SEARCH_RT_RANGE, msSearch);
    minParams.setParameter(io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.minimumsearch.MinimumSearchFeatureResolverParameters.PEAK_DURATION, com.google.common.collect.Range.closed(msDurMin, msDurMax));
    list.add(new ResolverSpec("MinimumSearch", Objects.requireNonNull(minParams.getResolver(minParams, flist))));

    // Savitzky-Golay
    final SavitzkyGolayFeatureResolverParameters sgParams = new SavitzkyGolayFeatureResolverParameters();
    sgParams.setParameter(GeneralResolverParameters.dimension, ResolvingDimension.RETENTION_TIME);
    sgParams.setParameter(GeneralResolverParameters.MIN_NUMBER_OF_DATAPOINTS, minPoints);
    // Permissive defaults with -D overrides
    final double sgMin = Double.parseDouble(System.getProperty("sgMin", "20"));
    final double sgDeriv = Double.parseDouble(System.getProperty("sgDeriv", "0.3"));
    final double sgDurMin = Double.parseDouble(System.getProperty("sgDurMin", "0.0"));
    final double sgDurMax = Double.parseDouble(System.getProperty("sgDurMax", "10.0"));
    sgParams.setParameter(SavitzkyGolayFeatureResolverParameters.MIN_PEAK_HEIGHT, sgMin);
    sgParams.setParameter(SavitzkyGolayFeatureResolverParameters.DERIVATIVE_THRESHOLD_LEVEL, sgDeriv);
    sgParams.setParameter(SavitzkyGolayFeatureResolverParameters.PEAK_DURATION, com.google.common.collect.Range.closed(sgDurMin, sgDurMax));
    list.add(new ResolverSpec("SavitzkyGolay", Objects.requireNonNull(sgParams.getResolver(sgParams, flist))));

    // Wavelet (peak detector) – processed like other resolvers
    final WaveletResolverParameters waveParams = new WaveletResolverParameters();
    waveParams.setParameter(GeneralResolverParameters.dimension, ResolvingDimension.RETENTION_TIME);
    waveParams.setParameter(GeneralResolverParameters.MIN_NUMBER_OF_DATAPOINTS, minPoints);
    // Permissive defaults with -D overrides
    final double waveSnr = Double.parseDouble(System.getProperty("waveSnr", "1.5"));
    final String waveTopToEdgeStr = System.getProperty("waveTopToEdge", "").trim();
    final double waveMinHeight = Double.parseDouble(System.getProperty("waveMinHeight", "0"));
    final String waveNoise = System.getProperty("waveNoise", "MAD").trim();
    final String waveScales = System.getProperty("waveScales", "").trim();
    final String waveKernel = System.getProperty("waveKernel", "").trim();
    final String waveNoiseWindow = System.getProperty("waveNoiseWindow", "").trim();
    waveParams.setParameter(WaveletResolverParameters.snr, waveSnr);
    if (!waveTopToEdgeStr.isEmpty()) {
      try { waveParams.setParameter(WaveletResolverParameters.topToEdge, true, Double.parseDouble(waveTopToEdgeStr)); } catch (Exception ignore) {}
    }
    waveParams.setParameter(WaveletResolverParameters.minHeight, waveMinHeight);
    {
      WaveletResolverParameters.NoiseCalculation nc = WaveletResolverParameters.NoiseCalculation.STANDARD_DEVIATION;
      final String w = waveNoise.toLowerCase(Locale.ROOT);
      if (w.startsWith("mad") || w.contains("median")) nc = WaveletResolverParameters.NoiseCalculation.MEDIAN_ABSOLUTE_DEVIATION;
      waveParams.setParameter(WaveletResolverParameters.noiseCalculation, nc);
    }
    // Optional advanced parameters (scales, kernel, noise window)
    if (!waveScales.isEmpty() || !waveKernel.isEmpty() || !waveNoiseWindow.isEmpty()) {
      final AdvancedWaveletParameters adv = new AdvancedWaveletParameters();
      if (!waveScales.isEmpty()) {
        adv.setParameter(AdvancedWaveletParameters.scales, true, waveScales);
      }
      if (!waveKernel.isEmpty()) {
        try { adv.setParameter(AdvancedWaveletParameters.WAVELET_KERNEL_RADIUS_FACTOR, true, Double.parseDouble(waveKernel)); } catch (Exception ignore) {}
      }
      if (!waveNoiseWindow.isEmpty()) {
        try { adv.setParameter(AdvancedWaveletParameters.LOCAL_NOISE_WINDOW_FACTOR, true, Double.parseDouble(waveNoiseWindow)); } catch (Exception ignore) {}
      }
      waveParams.setParameter(WaveletResolverParameters.advancedParameters, true);
      waveParams.getParameter(WaveletResolverParameters.advancedParameters).setEmbeddedParameters(adv);
    }
    list.add(new ResolverSpec("Wavelet", Objects.requireNonNull(waveParams.getResolver(waveParams, flist))));

    // CWT ridge
    final CwtRidgeResolverParameters cwtParams = new CwtRidgeResolverParameters();
    cwtParams.setParameter(GeneralResolverParameters.dimension, ResolvingDimension.RETENTION_TIME);
    cwtParams.setParameter(GeneralResolverParameters.MIN_NUMBER_OF_DATAPOINTS, minPoints);
    final double cwtSnr = Double.parseDouble(System.getProperty("cwtSnr", "0.8"));
    final double cwtBoundary = Double.parseDouble(System.getProperty("cwtBoundary", "0.05"));
    final int cwtPersist = Integer.parseInt(System.getProperty("cwtPersist", "1"));
    final double cwtMin = Double.parseDouble(System.getProperty("cwtMin", "0.0"));
    final boolean cwtAutoFlag = Boolean.parseBoolean(System.getProperty("cwtAuto", "false"));
    final String cwtQuant = System.getProperty("cwtQuant", "").trim();
    final String cwtQuantQ = System.getProperty("cwtQuantQ", "").trim();
    cwtParams.setParameter(io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.CwtRidgeResolverParameters.snrCwt, cwtSnr);
    cwtParams.setParameter(io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.CwtRidgeResolverParameters.boundaryThresholdFactor, cwtBoundary);
    cwtParams.setParameter(io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.CwtRidgeResolverParameters.ridgeMinPersistence, (double) cwtPersist);
    cwtParams.setParameter(io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.CwtRidgeResolverParameters.minApexHeight, cwtMin);
    cwtParams.setParameter(io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.CwtRidgeResolverParameters.autoTune, cwtAutoFlag);
    if (!cwtQuant.isEmpty()) {
      try {
        final var mode = io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.CwtRidgeResolverParameters.QuantBaselineMethod.valueOf(cwtQuant);
        cwtParams.setParameter(io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.CwtRidgeResolverParameters.quantBaselineMethod, mode);
        if (mode == io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.CwtRidgeResolverParameters.QuantBaselineMethod.LOW_QUANTILE) {
          final double q = cwtQuantQ.isEmpty() ? 0.1 : Double.parseDouble(cwtQuantQ);
          cwtParams.setParameter(io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.CwtRidgeResolverParameters.quantBaselineQuantile, true, q);
        } else {
          // disable optional quantile when NONE
          cwtParams.setParameter(io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.CwtRidgeResolverParameters.quantBaselineQuantile, false, 0.1);
        }
      } catch (Exception ignored) {}
    } else if (!cwtQuantQ.isEmpty()) {
      // Quantile provided without explicit mode → assume LOW_QUANTILE
      try {
        cwtParams.setParameter(io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.CwtRidgeResolverParameters.quantBaselineMethod,
            io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.CwtRidgeResolverParameters.QuantBaselineMethod.LOW_QUANTILE);
        final double q = Double.parseDouble(cwtQuantQ);
        cwtParams.setParameter(io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.CwtRidgeResolverParameters.quantBaselineQuantile, true, q);
      } catch (Exception ignored) {}
    }
    list.add(new ResolverSpec("CWT_Ridge", Objects.requireNonNull(cwtParams.getResolver(cwtParams, flist))));

    return list;
  }

  private static TsvEic ensureMonotonicRt(TsvEic eic) {
    // If RTs look like seconds (large), convert to minutes.
    final double median = percentile(eic.x, 0.5);
    final double[] xCopy = Arrays.copyOf(eic.x, eic.x.length);
    final double[] yCopy = Arrays.copyOf(eic.y, eic.y.length);
    if (median > 100.0) {
      for (int i = 0; i < xCopy.length; i++) xCopy[i] = xCopy[i] / 60.0;
    }
    // sort by x just in case input has slight disorder
    final int n = xCopy.length;
    final Integer[] idx = new Integer[n];
    for (int i = 0; i < n; i++) idx[i] = i;
    Arrays.sort(idx, Comparator.comparingDouble(i -> xCopy[i]));
    final double[] sx = new double[n];
    final double[] sy = new double[n];
    for (int i = 0; i < n; i++) { sx[i] = xCopy[idx[i]]; sy[i] = yCopy[idx[i]]; }
    return new TsvEic(sx, sy);
  }

  private static double percentile(double[] a, double q) {
    final double[] c = Arrays.copyOf(a, a.length);
    Arrays.sort(c);
    final double pos = q * (c.length - 1);
    final int lo = (int) Math.floor(pos), hi = (int) Math.ceil(pos);
    if (lo == hi) return c[lo];
    final double t = pos - lo;
    return c[lo] * (1 - t) + c[hi] * t;
  }

  private static TsvEic readEic(File file) throws IOException {
    final char sep = CSVParsingUtils.autoDetermineSeparatorDefaultFallback(file);
    try (BufferedReader br = Files.newBufferedReader(file.toPath())) {
      final List<String[]> rows = CSVParsingUtils.readData(br, String.valueOf(sep));
      if (rows.isEmpty()) return new TsvEic(new double[0], new double[0]);
      final String[] header = rows.get(0);
      int xIdx = guessIndex(header, "rt", "time", "x");
      int yIdx = guessIndex(header, "intensity", "y", "height", "signal");
      int startRow = 1;
      // heuristic if header lookup fails or if table has >2 columns
      if (xIdx < 0 || yIdx < 0) {
        startRow = 0; // treat first row as data
        final int cols = rows.stream().mapToInt(a -> a.length).max().orElse(0);
        // build numeric columns
        final List<double[]> colsData = new ArrayList<>();
        for (int c = 0; c < cols; c++) {
          final double[] col = new double[rows.size() - startRow];
          int k = 0;
          for (int r = startRow; r < rows.size(); r++) {
            final String[] row = rows.get(r);
            final String s = (c < row.length) ? safe(row[c]) : null;
            col[k++] = s == null ? Double.NaN : parseOrNaN(s);
          }
          colsData.add(col);
        }
        // find monotonic increasing column as X
        xIdx = -1;
        for (int c = 0; c < colsData.size(); c++) {
          if (isMostlyIncreasing(colsData.get(c)) && span(colsData.get(c)) > 1e-6) { xIdx = c; break; }
        }
        if (xIdx < 0 && cols >= 1) xIdx = 0;
        // find Y as column with max stddev (excluding x)
        double bestStd = -1; yIdx = -1;
        for (int c = 0; c < colsData.size(); c++) {
          if (c == xIdx) continue;
          double st = stddev(colsData.get(c));
          if (st > bestStd) { bestStd = st; yIdx = c; }
        }
        if (yIdx < 0 && cols >= 2) yIdx = 1;
      }

      final List<double[]> vals = new ArrayList<>(rows.size());
      for (int r = startRow; r < rows.size(); r++) {
        final String[] row = rows.get(r);
        if (row.length <= Math.max(xIdx, yIdx)) continue;
        final String xStr = safe(row[xIdx]);
        final String yStr = safe(row[yIdx]);
        if (xStr == null || yStr == null) continue;
        try {
          final double x = Double.parseDouble(xStr);
          final double y = Double.parseDouble(yStr);
          vals.add(new double[]{x, y});
        } catch (NumberFormatException ignored) { }
      }
      final double[] x = new double[vals.size()];
      final double[] y = new double[vals.size()];
      for (int i = 0; i < vals.size(); i++) {
        x[i] = vals.get(i)[0];
        y[i] = vals.get(i)[1];
      }
      return new TsvEic(x, y);
    } catch (Exception e) {
      throw new IOException("Failed to parse EIC TSV: " + file, e);
    }
  }

  private static boolean isMostlyIncreasing(double[] a) {
    int valid = 0, inc = 0;
    double prev = Double.NaN;
    for (double v : a) {
      if (Double.isNaN(v)) continue;
      if (!Double.isNaN(prev)) { valid++; if (v > prev) inc++; }
      prev = v;
    }
    return valid > 0 && inc >= 0.9 * valid; // 90% increasing
  }

  private static double span(double[] a) {
    double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;
    for (double v : a) if (!Double.isNaN(v)) { min = Math.min(min, v); max = Math.max(max, v); }
    return (min == Double.POSITIVE_INFINITY) ? 0 : (max - min);
  }

  private static double stddev(double[] a) {
    double sum = 0, sum2 = 0; int n = 0;
    for (double v : a) if (!Double.isNaN(v)) { sum += v; sum2 += v * v; n++; }
    if (n <= 1) return 0;
    double mean = sum / n; return Math.sqrt(Math.max(0, sum2 / n - mean * mean));
  }

  private static double parseOrNaN(String s) {
    try { return Double.parseDouble(s); } catch (Exception ex) { return Double.NaN; }
  }

  private static List<Range<Double>> readRanges(File file) throws IOException {
    final char sep = CSVParsingUtils.autoDetermineSeparatorDefaultFallback(file);
    try (BufferedReader br = Files.newBufferedReader(file.toPath())) {
      final List<String[]> rows = CSVParsingUtils.readData(br, String.valueOf(sep));
      if (rows.isEmpty()) return List.of();
      final String[] header = rows.get(0);
      int sIdx = guessIndex(header, "start", "rt_start", "begin", "left");
      int eIdx = guessIndex(header, "end", "rt_end", "right");
      int startRow = 1;
      if (sIdx < 0 || eIdx < 0) {
        // assume two columns without header
        sIdx = 0; eIdx = 1; startRow = 0;
      }
      final List<Range<Double>> ranges = new ArrayList<>();
      for (int r = startRow; r < rows.size(); r++) {
        final String[] row = rows.get(r);
        if (row.length <= Math.max(sIdx, eIdx)) continue;
        final String sStr = safe(row[sIdx]);
        final String eStr = safe(row[eIdx]);
        if (sStr == null || eStr == null) continue;
        try {
          double a = Double.parseDouble(sStr);
          double b = Double.parseDouble(eStr);
          if (a > b) { double t = a; a = b; b = t; }
          // assume units are same as EIC; convert seconds->minutes when clearly seconds
          if (Math.max(a, b) > 100.0) { a /= 60.0; b /= 60.0; }
          ranges.add(Range.closed(a, b));
        } catch (NumberFormatException ignored) { }
      }
      return ranges;
    } catch (Exception e) {
      throw new IOException("Failed to parse ranges TSV: " + file, e);
    }
  }

  private static int guessIndex(String[] header, String... keys) {
    final Map<String, Integer> map = new HashMap<>();
    for (int i = 0; i < header.length; i++) {
      if (header[i] == null) continue;
      map.put(header[i].trim().toLowerCase(Locale.ROOT), i);
    }
    for (String k : keys) {
      final Integer idx = map.get(k.toLowerCase(Locale.ROOT));
      if (idx != null) return idx;
    }
    // try contains
    for (String k : keys) {
      for (Map.Entry<String, Integer> e : map.entrySet()) {
        if (e.getKey().contains(k.toLowerCase(Locale.ROOT))) return e.getValue();
      }
    }
    return -1;
  }

  private static String safe(String s) {
    if (s == null) return null;
    s = s.trim();
    return s.isEmpty() ? null : s;
  }

  private static void saveChart(File out, String title, TsvEic eic, List<Range<Double>> labels,
      Map<String, List<Range<Double>>> resultsByResolver) throws IOException {
    final XYSeries series = new XYSeries("EIC");
    for (int i = 0; i < eic.x.length; i++) series.add(eic.x[i], eic.y[i]);
    final XYSeriesCollection ds = new XYSeriesCollection(series);

    final JFreeChart chart = ChartFactory.createXYLineChart(title, "RT (min)", "Intensity", ds,
        PlotOrientation.VERTICAL, true, false, false);
    final XYPlot plot = chart.getXYPlot();
    plot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);

    // Labeled ranges as background markers
    for (Range<Double> r : labels) {
      final IntervalMarker m = new IntervalMarker(r.lowerEndpoint(), r.upperEndpoint());
      m.setPaint(new Color(0, 128, 0, 40));
      m.setOutlinePaint(new Color(0, 128, 0, 120));
      m.setOutlineStroke(new BasicStroke(1f));
      plot.addDomainMarker(m, Layer.BACKGROUND);
    }

    // Predicted ranges for each resolver as translucent overlays (stacked)
    final LinkedHashMap<String, Color> colorMap = new LinkedHashMap<>();
    colorMap.put("MinimumSearch", new Color(0, 0, 255, 50));
    colorMap.put("SavitzkyGolay", new Color(255, 0, 0, 50));
    colorMap.put("CWT_Ridge", new Color(255, 165, 0, 50));
    final Color[] fallback = new Color[]{new Color(128, 0, 128, 50), new Color(0, 128, 128, 50)};
    int fb = 0;
    final double yMax = Arrays.stream(eic.y).max().orElse(1d);
    for (Map.Entry<String, List<Range<Double>>> e : resultsByResolver.entrySet()) {
      final Color c = colorMap.containsKey(e.getKey()) ? colorMap.get(e.getKey()) : fallback[(fb++) % fallback.length];
      for (Range<Double> r : e.getValue()) {
        final Rectangle2D rect = new Rectangle2D.Double(r.lowerEndpoint(), 0, r.upperEndpoint() - r.lowerEndpoint(), yMax);
        final XYShapeAnnotation ann = new XYShapeAnnotation(rect, new BasicStroke(0f), new Color(0, 0, 0, 0), c);
        plot.addAnnotation(ann);
      }
    }

    // Subtitle with resolvers shown
    final String legend = resultsByResolver.keySet().stream().collect(Collectors.joining(", "));
    chart.addSubtitle(new org.jfree.chart.title.TextTitle("Resolvers shown: " + legend));

    ChartUtils.saveChartAsPNG(out, chart, 1200, 600);
  }

  private static void saveHtml(File out, String title, TsvEic eic, List<Range<Double>> labels,
      Map<String, List<Range<Double>>> resultsByResolver) throws IOException {
    final String plotly = "https://cdn.plot.ly/plotly-2.27.0.min.js";
    final StringBuilder html = new StringBuilder();
    html.append("<!doctype html><html><head><meta charset=\"utf-8\"><title>")
        .append(title).append("</title><script src=\"").append(plotly).append("\"></script></head><body>");
    html.append("<div id=\"plot\" style=\"width:1200px;height:600px;\"></div>");
    html.append("<div style=\"font-family:Sans-Serif;margin-top:8px\">");
    html.append("<label><input type=\"checkbox\" id=\"labels\" checked> Show labels</label> ");
    for (String name : resultsByResolver.keySet()) {
      html.append("<label style=\"margin-left:12px\"><input type=\"checkbox\" class=\"resolver\" data-name=\"")
          .append(name).append("\" checked> ")
          .append(name).append("</label>");
    }
    html.append("</div>");
    // counts
    html.append("<div style=\"font-family:Sans-Serif;margin-top:6px;color:#555\" id=\"counts\"></div>");
    html.append("<script>\n");
    html.append("const x=").append(toJsonArray(eic.x)).append(";\n");
    html.append("const y=").append(toJsonArray(eic.y)).append(";\n");
    html.append("const labels=").append(toJsonRanges(labels)).append(";\n");
    html.append("const preds={");
    int k=0; int sz=resultsByResolver.size();
    for (Map.Entry<String,List<Range<Double>>> e : resultsByResolver.entrySet()) {
      html.append(quote(e.getKey())).append(": ").append(toJsonRanges(e.getValue()));
      if (++k<sz) html.append(',');
    }
    html.append("};\n");
    html.append("const colors={MinimumSearch:'rgba(0,0,255,0.25)',SavitzkyGolay:'rgba(255,0,0,0.25)',CWT_Ridge:'rgba(255,165,0,0.25)'};\n");
    html.append("function shapes(showLabels,visible){const s=[]; if(showLabels){for(const r of labels){s.push({type:'rect',xref:'x',yref:'paper',x0:r[0],x1:r[1],y0:0,y1:1,fillcolor:'rgba(0,128,0,0.20)',line:{width:0}});} } for(const [name,rs] of Object.entries(preds)){ if(!visible.has(name)) continue; const c=colors[name]||'rgba(128,0,128,0.25)'; for(const r of rs){s.push({type:'rect',xref:'x',yref:'paper',x0:r[0],x1:r[1],y0:0,y1:1,fillcolor:c,line:{width:0}});} } return s; }\n");
    html.append("function rgbaToSolid(rgba){if(!rgba) return '#000'; const m=rgba.match(/rgba\\((\\d+),(\\d+),(\\d+),/); return m?`rgb(${m[1]},${m[2]},${m[3]})`:'#000';}\n");
    html.append("function maskTrace(name){const rs=preds[name]||[]; const y2=new Array(y.length).fill(null); let j=0; for(const r of rs){const x0=r[0],x1=r[1]; while(j<x.length && x[j]<x0) j++; let k=j; while(k<x.length && x[k]<=x1){ y2[k]=y[k]; k++; } } return {x:x,y:y2,type:'scatter',mode:'lines',name:name,line:{color:rgbaToSolid(colors[name])}};}\n");
    html.append("function buildData(visible){const arr=[{x:x,y:y,type:'scatter',mode:'lines',name:'EIC',line:{color:'#000'}}]; for(const name of Object.keys(preds)){ if(visible.has(name)) arr.push(maskTrace(name)); } return arr;}\n");
    html.append("const vis=new Set(Object.keys(preds));\n");
    html.append("let layout={title:").append(quote(title)).append(",xaxis:{title:'RT (min)'},yaxis:{title:'Intensity'},shapes:shapes(true,vis)};\n");
    html.append("function updateCounts(){const el=document.getElementById('counts'); const nLab=labels.length; let s='Labels: '+nLab+' | '; for(const [n,rs] of Object.entries(preds)){ s+= n+': '+rs.length+' '; } el.textContent=s.trim(); }\n");
    html.append("Plotly.newPlot('plot',buildData(vis),layout); updateCounts();\n");
    html.append("document.getElementById('labels').addEventListener('change',ev=>{layout.shapes=shapes(ev.target.checked,vis); Plotly.react('plot',buildData(vis),layout);});\n");
    html.append("for(const cb of document.querySelectorAll('.resolver')){cb.addEventListener('change',ev=>{const n=ev.target.dataset.name; if(ev.target.checked) vis.add(n); else vis.delete(n); layout.shapes=shapes(document.getElementById('labels').checked,vis); Plotly.react('plot',buildData(vis),layout);});}\n");
    html.append("</script></body></html>");
    Files.writeString(out.toPath(), html.toString());
  }

  private static String quote(String s){ return "\""+s.replace("\"","\\\"")+"\""; }
  private static String toJsonArray(double[] a){ StringBuilder sb=new StringBuilder("["); for(int i=0;i<a.length;i++){ if(i>0) sb.append(','); sb.append(Double.toString(a[i])); } return sb.append(']').toString(); }
  private static String toJsonRanges(List<Range<Double>> rs){ StringBuilder sb=new StringBuilder("["); for(int i=0;i<rs.size();i++){ if(i>0) sb.append(','); Range<Double> r=rs.get(i); sb.append('[').append(r.lowerEndpoint()).append(',').append(r.upperEndpoint()).append(']'); } return sb.append(']').toString(); }

  /**
   * Build an IonTimeSeries and call the resolver via AbstractResolver.resolve(series, storage)
   * to mimic MZmine's normal pathway where resolvers split a series into sub-series.
   */
  private static List<Range<Double>> predictRanges(ResolverSpec spec, TsvEic eic, RawDataFileImpl raw) {
    try {
      // build scans
      final int n = eic.x.length;
      final double[] iy = "CWT_Ridge".equals(spec.name) ? Arrays.copyOf(eic.y, n) : preprocessIntensities(eic.y);
      final List<SimpleScan> scans = new ArrayList<>(n);
      for (int i = 0; i < n; i++) {
        // simple centroid scan with single (mz,intensity) point; mz can be constant, not used here
        final double[] mz = new double[]{0d};
        final double[] inten = new double[]{iy[i]};
        final SimpleScan sc = new SimpleScan(raw, i+1, 1, (float)eic.x[i], null, mz, inten,
            MassSpectrumType.CENTROIDED, PolarityType.UNKNOWN, "EIC", null);
        scans.add(sc);
      }
      final SimpleIonTimeSeries series = new SimpleIonTimeSeries(raw.getMemoryMapStorage(), new double[n], iy, scans);
      final var sub = spec.resolver.resolve(series, raw.getMemoryMapStorage());
      final List<Range<Double>> out = new ArrayList<>();
      for (var s : sub) {
        if (s.getNumberOfValues() == 0) continue;
        final double start = s.getRetentionTime(0);
        final double end = s.getRetentionTime(s.getNumberOfValues()-1);
        out.add(Range.closed(Math.min(start,end), Math.max(start,end)));
      }
      // Prefer series-based resolution. If nothing found, fallback to direct xy-based.
      if (!out.isEmpty()) return out;
      return spec.resolver.resolve(eic.x, eic.y);
    } catch (Throwable t) {
      // fallback to direct resolve(x,y)
      return spec.resolver.resolve(eic.x, eic.y);
    }
  }

  private static double quantile(double[] a, double q) {
    final double[] c = Arrays.copyOf(a, a.length);
    Arrays.sort(c);
    final double pos = q * (c.length - 1);
    final int lo = (int) Math.floor(pos), hi = (int) Math.ceil(pos);
    if (lo == hi) return c[lo];
    final double t = pos - lo;
    return c[lo] * (1 - t) + c[hi] * t;
  }

  private static double[] preprocessIntensities(double[] y) {
    final int n = y.length;
    final double[] base = new double[n];
    final double frac = Double.parseDouble(System.getProperty("baselineFrac", "0.05"));
    final int w = Math.max(3, (int) Math.round(n * frac));
    final int half = w / 2;
    // simple moving average baseline
    for (int i = 0; i < n; i++) {
      int a = Math.max(0, i - half), b = Math.min(n - 1, i + half);
      double sum = 0; int cnt = 0;
      for (int j = a; j <= b; j++) { sum += y[j]; cnt++; }
      base[i] = cnt > 0 ? sum / cnt : 0d;
    }
    final double[] detr = new double[n];
    for (int i = 0; i < n; i++) detr[i] = Math.max(0d, y[i] - base[i]);
    // pre-zero small intensities to emulate background removal
    final double preQ = Double.parseDouble(System.getProperty("preZeroQ", "0.85"));
    final double thr = quantile(detr, preQ);
    for (int i = 0; i < n; i++) if (detr[i] < thr) detr[i] = 0d;
    return detr;
  }

  private static Metrics evaluate(List<Range<Double>> labels, List<Range<Double>> preds,
      double iouThr, double rtTol) {
    final boolean[] labelUsed = new boolean[labels.size()];
    int tp = 0, fp = 0;
    for (Range<Double> p0 : preds) {
      final Range<Double> p = expandRange(p0, rtTol);
      int bestIdx = -1; double bestIoU = 0d;
      for (int i = 0; i < labels.size(); i++) {
        if (labelUsed[i]) continue;
        final double iou = iou(expandRange(labels.get(i), rtTol), p);
        if (iou > bestIoU) { bestIoU = iou; bestIdx = i; }
      }
      if (bestIdx >= 0 && bestIoU >= iouThr) {
        tp++; labelUsed[bestIdx] = true;
      } else {
        fp++;
      }
    }
    int fn = 0; for (boolean u : labelUsed) if (!u) fn++;
    final double precision = tp + fp == 0 ? 0 : ((double) tp) / (tp + fp);
    final double recall = tp + fn == 0 ? 0 : ((double) tp) / (tp + fn);
    final double f1 = precision + recall == 0 ? 0 : 2 * precision * recall / (precision + recall);
    return new Metrics(tp, fp, fn, precision, recall, f1);
  }

  private static Range<Double> expandRange(Range<Double> r, double tol) {
    return Range.closed(r.lowerEndpoint() - tol, r.upperEndpoint() + tol);
  }

  private static double iou(Range<Double> a, Range<Double> b) {
    final double left = Math.max(a.lowerEndpoint(), b.lowerEndpoint());
    final double right = Math.min(a.upperEndpoint(), b.upperEndpoint());
    if (right <= left) return 0d;
    final double inter = right - left;
    final double uni = (a.upperEndpoint() - a.lowerEndpoint()) + (b.upperEndpoint() - b.lowerEndpoint()) - inter;
    if (uni <= 0) return 0d;
    return inter / uni;
  }

  private static String format(double v) {
    return new DecimalFormat("0.###").format(v);
  }

  private record TsvEic(double[] x, double[] y) {}

  private record ResolverSpec(String name, Resolver resolver) {
    @Override public String toString() { return name; }
  }

  private record Metrics(int tp, int fp, int fn, double precision, double recall, double f1) {}

  private static String parameterSummary(int minPoints, double iou, double rtTol) {
    final StringBuilder sb = new StringBuilder();
    sb.append("Evaluation settings\n");
    sb.append("- IoU: ").append(iou).append('\n');
    sb.append("- rtTol (min): ").append(rtTol).append('\n');
    sb.append("- minPoints: ").append(minPoints).append('\n');
    sb.append('\n');
    sb.append("Resolver defaults / overrides (-D...)\n");
    sb.append("MinimumSearch: minAbs=").append(System.getProperty("minAbs", "50"))
        .append(", minRel=").append(System.getProperty("minRel", "0.01"))
        .append(", ratio=").append(System.getProperty("ratio", "1.3"))
        .append(", thresh=").append(System.getProperty("msThresh", "0.5"))
        .append(", searchRange=").append(System.getProperty("searchRange", "0.03")).append('\n');
    sb.append("SavitzkyGolay: minHeight=").append(System.getProperty("sgMin", "50"))
        .append(", derivThresh=").append(System.getProperty("sgDeriv", "0.5")).append('\n');
    sb.append("CWT_Ridge: snr=").append(System.getProperty("cwtSnr", "2.0"))
        .append(", boundaryFraction=").append(System.getProperty("cwtBoundary", "0.05"))
        .append(", minApexHeight=").append(System.getProperty("cwtMin", "0.0"))
        .append(", ridgeMinPersistence=").append(System.getProperty("cwtPersist", "1")).append('\n');
    return sb.toString();
  }
}


