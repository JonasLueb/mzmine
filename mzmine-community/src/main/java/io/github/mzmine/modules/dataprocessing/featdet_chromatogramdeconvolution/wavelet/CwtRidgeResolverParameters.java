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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet;

import io.github.mzmine.datamodel.features.ModularFeatureList;
import io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.FeatureResolverSetupDialog;
import io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.GeneralResolverParameters;
import io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.Resolver;
// type is referenced in getResolver
import io.github.mzmine.parameters.ParameterSet;
import io.github.mzmine.parameters.parametertypes.AdvancedParametersParameter;
import io.github.mzmine.parameters.parametertypes.BooleanParameter;
import io.github.mzmine.parameters.parametertypes.ComboParameter;
import io.github.mzmine.parameters.parametertypes.DoubleParameter;
import io.github.mzmine.parameters.parametertypes.OptionalParameter;
import io.github.mzmine.util.ExitCode;
import java.text.DecimalFormat;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

/**
 * Parameters for the CWT ridge-persistence based resolver. Designed to coexist with the existing
 * WaveletPeakDetector.
 */
public class CwtRidgeResolverParameters extends GeneralResolverParameters {

  public static final DoubleParameter snrCwt = new DoubleParameter("CWT SNR threshold",
      "Minimum CWT-domain SNR for maxima (relative to MAD at smallest scale).",
      new DecimalFormat("#.#"), 3d, 0d, Double.MAX_VALUE);

  public static final DoubleParameter minApexHeight = new DoubleParameter("Minimum apex height",
      "Minimum apex intensity in the raw EIC.", new DecimalFormat("0.###E0"), 1e3, 0d,
      Double.MAX_VALUE);

  public static final DoubleParameter boundaryThresholdFactor = new DoubleParameter(
      "Boundary threshold (fallback)",
      "Fallback boundary as a fraction of apex height if CWT zero-crossings are not found.",
      new DecimalFormat("0.00"), 0.05, 0d, 1d);

  public static final DoubleParameter ridgeMinPersistence = new DoubleParameter(
      "Min ridge length (scales)",
      "Minimum number of adjacent scales a maxima must persist to be considered a ridge.",
      new DecimalFormat("#"), 3d, 1d, 100d);

  public static final OptionalParameter<DoubleParameter> ridgeNeighborhood = new OptionalParameter<>(
      new DoubleParameter("Ridge neighborhood (points)",
          "Max lateral movement of the ridge per scale (in points).",
          new DecimalFormat("#"), 2d, 0d, 1000d));

  public static final BooleanParameter enableDriftAwareDetrend = new BooleanParameter(
      "Enable drift-aware detrend (detection only)",
      "If enabled, apply a gentle LOESS detrend only if pronounced low-frequency drift is detected.",
      false);

  public enum QuantBaselineMethod {
    LOW_QUANTILE,
    NONE
  }

  public static final ComboParameter<QuantBaselineMethod> quantBaselineMethod = new ComboParameter<>(
      "Quantification baseline", "Method to estimate local baseline for area integration.",
      QuantBaselineMethod.values(), QuantBaselineMethod.LOW_QUANTILE);

  public static final OptionalParameter<DoubleParameter> quantBaselineQuantile = new OptionalParameter<>(
      new DoubleParameter("Baseline quantile",
          "Quantile for LOW_QUANTILE baseline within the peak window (0..1)",
          new DecimalFormat("0.00"), 0.1, 0d, 1d));

  public static final BooleanParameter autoTune = new BooleanParameter("Auto-tune parameters",
      "Enable per-chromatogram heuristic tuning of CWT parameters (scales, SNR, persistence, baseline, apex).",
      false);

  public static final AdvancedParametersParameter<AdvancedWaveletParameters> advanced =
      new AdvancedParametersParameter<>(new AdvancedWaveletParameters());

  public CwtRidgeResolverParameters() {
    super(PEAK_LISTS, dimension, groupMS2Parameters, snrCwt, minApexHeight,
        boundaryThresholdFactor, ridgeMinPersistence, ridgeNeighborhood, enableDriftAwareDetrend,
        quantBaselineMethod, quantBaselineQuantile, autoTune, MIN_NUMBER_OF_DATAPOINTS, SUFFIX,
        handleOriginal, advanced);
  }

  @Nullable
  @Override
  public Resolver getResolver(ParameterSet parameterSet, ModularFeatureList flist) {
    return new CwtRidgeResolver(parameterSet, flist);
  }

  @NotNull
  @Override
  public ExitCode showSetupDialog(boolean valueCheckRequired) {
    final FeatureResolverSetupDialog dialog = new FeatureResolverSetupDialog(valueCheckRequired,
        this, null);
    dialog.showAndWait();
    return dialog.getExitCode();
  }
}


