package io.github.mzmine.modules.dataprocessing.id_ms2deepscore_vectorsearch;

import io.github.mzmine.util.spectraldb.entry.SpectralDBAnnotation;
import io.github.mzmine.util.scans.similarity.SpectralSimilarity;
import java.util.Objects;

/**
 * Helper methods to calculate rescored MS2Deepscore heuristics while retaining the original ANN
 * similarity score.
 */
public final class MS2DeepscoreScoringUtils {

  private MS2DeepscoreScoringUtils() {
  }

  /**
   * Compute a heuristic score that combines the ANN similarity with classical spectrum support
   * metrics. The individual terms are intentionally bounded between 0 and 1 so chemists can reason
   * about their contributions.
   */
  public static double computeRescore(final SpectralDBAnnotation annotation) {
    if (annotation == null) {
      return 0d;
    }

    final SpectralSimilarity sim = annotation.getSimilarity();
    if (sim == null) {
      return 0d;
    }

    final double sAnn = clamp01(sim.getScore());
    final double matchedSignals = saturatingSupport(sim.getOverlap());
    final double explained = clamp01(sim.getExplainedLibraryIntensity());

    // base ANN weight keeps ranking familiar, support nudges trustworthy matches upwards
    double rescored = 0.60 * sAnn + 0.20 * matchedSignals + 0.15 * explained;

    // penalise spectra with no real overlap/explained intensity
    if (sim.getOverlap() <= 0 || explained < 0.05) {
      rescored -= 0.30;
    }

    return clamp01(rescored);
  }

  private static double saturatingSupport(final int overlap) {
    if (overlap <= 0) {
      return 0d;
    }
    // Approaches 1.0 around 6–8 matched signals.
    return 1d - Math.exp(-(double) overlap / 6d);
  }

  private static double clamp01(final double value) {
    return Math.max(0d, Math.min(1d, value));
  }
}

