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

import io.github.mzmine.parameters.Parameter;
import io.github.mzmine.parameters.impl.SimpleParameterSet;
import io.github.mzmine.parameters.parametertypes.IntegerParameter;
import io.github.mzmine.parameters.parametertypes.DoubleParameter;
import io.github.mzmine.parameters.parametertypes.StringParameter;
import io.github.mzmine.parameters.parametertypes.filenames.FileNameParameter;
import io.github.mzmine.parameters.parametertypes.filenames.FileSelectionType;
import io.github.mzmine.parameters.parametertypes.BooleanParameter;
import io.github.mzmine.parameters.parametertypes.selectors.FeatureListsParameter;
import io.github.mzmine.parameters.parametertypes.tolerances.MZToleranceParameter;
import io.github.mzmine.modules.dataprocessing.filter_scan_merge_select.SpectraMergeSelectParameter;
import io.github.mzmine.parameters.parametertypes.combowithinput.MsLevelFilter;
import io.github.mzmine.parameters.parametertypes.combowithinput.MsLevelFilterParameter;

public class MS2DeepscoreVectorSearchParameters extends SimpleParameterSet {

  public static final FeatureListsParameter PEAK_LISTS = new FeatureListsParameter();

  public static final MsLevelFilterParameter MS_LEVEL = new MsLevelFilterParameter(
      new MsLevelFilter.Options[]{MsLevelFilter.Options.MS1, MsLevelFilter.Options.MS2},
      new MsLevelFilter(MsLevelFilter.Options.MS2));

  public static final SpectraMergeSelectParameter SPECTRA_SELECT = SpectraMergeSelectParameter.createSpectraLibrarySearchDefaultNoMSn();

  public static final FileNameParameter MODEL_FILE = new FileNameParameter(
      "MS2Deepscore model (.pt)",
      "Path to the PyTorch model file from the dual-mode release.",
      FileSelectionType.OPEN, true);

  public static final FileNameParameter SETTINGS_FILE = new FileNameParameter(
      "Model settings (settings.json)",
      "Path to the model settings JSON from the same release.",
      FileSelectionType.OPEN, true);

  public static final StringParameter DB_URI = new StringParameter("Vector DB URI",
      "Connection string or file path to the vector database (e.g., Milvus/SQLite).", "");

  public static final StringParameter COLLECTION_POS = new StringParameter("Collection (positive)",
      "Collection for positive ion mode (e.g., ms2ds_lc_esi_positive).", "ms2ds_lc_esi_positive");

  public static final StringParameter COLLECTION_NEG = new StringParameter("Collection (negative)",
      "Collection for negative ion mode (e.g., ms2ds_lc_esi_negative).", "ms2ds_lc_esi_negative");

  public static final IntegerParameter TOP_K = new IntegerParameter("Top K",
      "Number of nearest neighbors to retrieve per query.", 20);

  public static final DoubleParameter MIN_SCORE = new DoubleParameter("Min similarity",
      "Minimum cosine similarity to retain from ANN results.",
      java.text.NumberFormat.getNumberInstance(), 0.7);

  public static final MZToleranceParameter PRECURSOR_TOL = new MZToleranceParameter(
      "Precursor m/z filter",
      "Optional candidate prefilter in the vector DB using precursor m/z.", 0.01, 25);

  public static final BooleanParameter CONFIRM_WITH_COSINE = new BooleanParameter(
      "Confirm with cosine", "Optionally fetch hit spectra and compute conventional cosine.",
      false);

  public MS2DeepscoreVectorSearchParameters() {
    super(new Parameter[]{PEAK_LISTS, MS_LEVEL, SPECTRA_SELECT, MODEL_FILE, SETTINGS_FILE, DB_URI,
        COLLECTION_POS, COLLECTION_NEG, TOP_K, MIN_SCORE, PRECURSOR_TOL, CONFIRM_WITH_COSINE});
    // set sensible defaults for the local dev environment
    DB_URI.setValue("milvus://localhost:19530?db=mz_connect_mvp");
    MODEL_FILE.setValue(new java.io.File("build/generated/ms2ds_dual/ms2deepscore_model_java.pt"));
    SETTINGS_FILE.setValue(new java.io.File("build/generated/ms2ds_dual/settings.json"));
  }
}


