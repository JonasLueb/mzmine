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

import com.google.common.collect.Range;
import io.github.mzmine.datamodel.Scan;
import io.github.mzmine.datamodel.featuredata.IonTimeSeries;
import io.github.mzmine.datamodel.features.ModularFeatureList;
import io.github.mzmine.modules.MZmineModule;
import io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.AbstractResolver;
import io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.GeneralResolverParameters;
import io.github.mzmine.parameters.ParameterSet;
import io.github.mzmine.util.MathUtils;
import io.github.mzmine.util.RangeUtils;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
// no logger used at the moment
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.jetbrains.annotations.NotNull;

/**
 * Ridge-persistence CWT resolver. Uses Ricker (Mexican hat) wavelet, adaptive MAD thresholds at the
 * smallest scale, traces ridges across adjacent scales, merges maxima across scales to single peak
 * candidates, and derives bounds by CWT zero-crossings/curvature with fallback.
 */
public class CwtRidgeResolver extends AbstractResolver {

  private final double[] scales;
  private final double snrCwt;
  private final double minApexHeight;
  private final double boundaryThresholdFraction;
  private final int minPoints;
  private final int ridgeMinPersistence;
  private final int ridgeNeighborhoodMaxMove;
  private final CwtRidgeResolverParameters.QuantBaselineMethod quantBaseline;
  private final double quantBaselineQuantile;

  public CwtRidgeResolver(@NotNull final ParameterSet parameters,
      @NotNull final ModularFeatureList flist) {
    super(parameters, flist);

    final var adv = parameters.getParameter(CwtRidgeResolverParameters.advanced);
    final String scaleStr = adv.getValueOrDefault(AdvancedWaveletParameters.scales,
        AdvancedWaveletParameters.DEFAULT_SCALES);
    boolean auto = parameters.getValue(CwtRidgeResolverParameters.autoTune);
    double[] tmpScales = Arrays.stream(scaleStr.split(",")).map(String::trim)
        .mapToDouble(Double::parseDouble).sorted().toArray();
    double tmpSnr = parameters.getValue(CwtRidgeResolverParameters.snrCwt);
    double tmpMinApex = parameters.getValue(CwtRidgeResolverParameters.minApexHeight);
    double tmpBoundary = parameters.getValue(CwtRidgeResolverParameters.boundaryThresholdFactor);
    this.minPoints = parameters.getValue(GeneralResolverParameters.MIN_NUMBER_OF_DATAPOINTS);
    int tmpPersist = (int) Math.round(parameters.getValue(CwtRidgeResolverParameters.ridgeMinPersistence));
    this.ridgeNeighborhoodMaxMove = parameters.getEmbeddedParameterValueIfSelectedOrElse(
        CwtRidgeResolverParameters.ridgeNeighborhood, 2d).intValue();
    CwtRidgeResolverParameters.QuantBaselineMethod tmpQuantBaseline = parameters.getValue(CwtRidgeResolverParameters.quantBaselineMethod);
    double tmpQuantQ = parameters.getEmbeddedParameterValueIfSelectedOrElse(
        CwtRidgeResolverParameters.quantBaselineQuantile, 0.1);

    if (auto) {
      // Lazy auto-tuning placeholders; compute from a minimal synthetic window at runtime in resolve()
      this.scales = tmpScales; // will be refined in resolve(double[],double[])
      this.snrCwt = Math.max(0.8, tmpSnr);
      this.minApexHeight = Math.max(0.0, tmpMinApex);
      this.boundaryThresholdFraction = Math.max(0.03, Math.min(0.1, tmpBoundary));
      this.ridgeMinPersistence = Math.max(1, tmpPersist);
      this.quantBaseline = tmpQuantBaseline;
      this.quantBaselineQuantile = tmpQuantQ;
    } else {
      this.scales = tmpScales;
      this.snrCwt = tmpSnr;
      this.minApexHeight = tmpMinApex;
      this.boundaryThresholdFraction = tmpBoundary;
      this.ridgeMinPersistence = tmpPersist;
      this.quantBaseline = tmpQuantBaseline;
      this.quantBaselineQuantile = tmpQuantQ;
    }
  }

  @Override
  public @NotNull Class<? extends MZmineModule> getModuleClass() {
    return CwtRidgeResolverModule.class;
  }

  @Override
  public @NotNull <T extends IonTimeSeries<? extends Scan>> List<T> resolve(@NotNull T series,
      io.github.mzmine.util.MemoryMapStorage storage) {
    final List<T> sub = super.resolve(series, storage);
    final List<T> adjusted = new ArrayList<>(sub.size());
    for (T s : sub) {
      final int n = s.getNumberOfValues();
      final double[] intens = new double[n];
      for (int i = 0; i < n; i++) intens[i] = s.getIntensity(i);
      double baseline = 0d;
      if (quantBaseline == CwtRidgeResolverParameters.QuantBaselineMethod.LOW_QUANTILE) {
        baseline = MathUtils.calcQuantile(intens, quantBaselineQuantile);
      }
      final double[] adjustedIntens = new double[n];
      for (int i = 0; i < n; i++) adjustedIntens[i] = Math.max(0d, intens[i] - baseline);
      @SuppressWarnings("unchecked")
      final T replaced = (T) s.copyAndReplace(storage, adjustedIntens);
      adjusted.add(replaced);
    }
    return adjusted;
  }

  @Override
  public @NotNull List<Range<Double>> resolve(double[] x, double[] y) {
    if (x == null || y == null || x.length != y.length || x.length < 5) {
      return List.of();
    }
    // If auto tuning enabled, derive CWT parameters from current series
    if (generalParameters.getValue(CwtRidgeResolverParameters.autoTune)) {
      final Auto auto = tuneFromSeries(x, y);
      // overwrite local runtime params
      System.arraycopy(auto.scales, 0, this.scales, 0, Math.min(this.scales.length, auto.scales.length));
      // snrCwt, boundary, persist, apex, baseline quant are applied via local variables below
      final double snr = auto.snrCwt;
      final double boundary = auto.boundaryFraction;
      final int persist = auto.ridgeMinPersistence;
      final double minApex = auto.minApexHeight;
      final CwtRidgeResolverParameters.QuantBaselineMethod qbm = auto.useLowQuantileBaseline ? CwtRidgeResolverParameters.QuantBaselineMethod.LOW_QUANTILE : CwtRidgeResolverParameters.QuantBaselineMethod.NONE;
      final double q = auto.quantileQ;
      return resolveWithParams(x, y, auto.scales, snr, minApex, boundary, persist, qbm, q);
    }

    // 1) Optional gentle detrend for detection only (not for quant), if strong drift detected
    final double[] yDetect = maybeDetrendForDetection(x, y);

    // 2) Compute CWT at all scales
    final Map<Double, double[]> cwt = computeCwt(yDetect, scales);
    if (cwt.isEmpty()) return List.of();

    // 3) MAD at smallest scale -> adaptive threshold
    final double[] smallest = cwt.get(scales[0]);
    final double noiseSigma = Math.max(estimateMadSigma(smallest), 1e-9);
    final double absoluteCoeffThreshold = snrCwt * noiseSigma;

    // 4) Local maxima per scale above threshold
    final Map<Integer, List<Integer>> maximaByIndexScale = new HashMap<>();
    for (double scale : scales) {
      final double[] coeff = cwt.get(scale);
      final List<Integer> maxima = findLocalMaxima(coeff, absoluteCoeffThreshold);
      for (int idx : maxima) {
        List<Integer> list = maximaByIndexScale.computeIfAbsent(idx, k -> new ArrayList<>());
        list.add((int) Math.round(scale));
      }
    }

    if (maximaByIndexScale.isEmpty()) return List.of();

    // 5) Ridge tracing across scales (greedy linking with neighborhood constraint)
    final List<Ridge> ridges = traceRidges(cwt, maximaByIndexScale, ridgeNeighborhoodMaxMove);

    // 6) Filter by ridge persistence and apex height in raw EIC
    final List<Integer> apexIndicesRaw = new ArrayList<>();
    for (Ridge r : ridges) {
      if (r.length() >= ridgeMinPersistence) {
        final int apexIdx = r.argmaxIndex(cwt);
        if (apexIdx >= 0 && y[apexIdx] >= minApexHeight) {
          apexIndicesRaw.add(apexIdx);
        }
      }
    }
    if (apexIndicesRaw.isEmpty()) return List.of();

    // 6b) De-duplicate very close apex indices (suppress neighbors around stronger apex)
    final List<int[]> apexWithScore = new ArrayList<>(); // [index, score*1e6 to sort as int]
    for (int idx : apexIndicesRaw) {
      double best = 0d;
      for (double[] coeff : cwt.values()) {
        final double v = Math.abs(coeff[idx]);
        if (v > best) best = v;
      }
      apexWithScore.add(new int[]{idx, (int) Math.round(best * 1_000_000d)});
    }
    apexWithScore.sort((a, b) -> Integer.compare(b[1], a[1])); // by score desc
    final boolean[] taken = new boolean[y.length];
    final int suppress = Math.max(1, ridgeNeighborhoodMaxMove);
    final List<Integer> apexIndices = new ArrayList<>();
    for (int[] cand : apexWithScore) {
      final int idx = cand[0];
      boolean conflict = false;
      final int a = Math.max(0, idx - suppress), b = Math.min(y.length - 1, idx + suppress);
      for (int i = a; i <= b; i++) { if (taken[i]) { conflict = true; break; } }
      if (!conflict) {
        apexIndices.add(idx);
        for (int i = a; i <= b; i++) taken[i] = true;
      }
    }
    if (apexIndices.isEmpty()) return List.of();
    apexIndices.sort(Integer::compareTo);

    // 7) Derive bounds via nearest CWT zero-crossings at characteristic scale; fallback to raw threshold
    final List<Range<Double>> ranges = new ArrayList<>();
    final boolean[] used = new boolean[y.length];
    for (int apex : apexIndices) {
      if (used[apex]) continue;
      int charScale = selectCharacteristicScale(cwt, apex);
      int left = findZeroCrossingLeft(cwt.get((double) charScale), apex);
      int right = findZeroCrossingRight(cwt.get((double) charScale), apex, y.length);

      if (left < 0 || right < 0 || right <= left) {
        // fallback using boundary threshold on raw signal
        final double thr = y[apex] * boundaryThresholdFraction;
        left = walkLeftByThreshold(y, apex, thr);
        right = walkRightByThreshold(y, apex, thr);
      }
      if (left < 0 || right < 0 || right <= left) continue;

      // refine to contiguous-above-threshold, unimodal support around the apex
      final int[] lr = trimBoundsUnimodal(y, apex, left, right, boundaryThresholdFraction, minPoints);
      if (lr == null) continue;
      left = lr[0]; right = lr[1];

      // Apex refinement: re-center on raw maximum within current bounds and re-derive bounds
      final int apexRef = argmax(y, left, right);
      if (apexRef != apex) {
        apex = apexRef;
        charScale = selectCharacteristicScale(cwt, apex);
        int l2 = findZeroCrossingLeft(cwt.get((double) charScale), apex);
        int r2 = findZeroCrossingRight(cwt.get((double) charScale), apex, y.length);
        if (l2 < 0 || r2 < 0 || r2 <= l2) {
          final double thr2 = y[apex] * boundaryThresholdFraction;
          l2 = walkLeftByThreshold(y, apex, thr2);
          r2 = walkRightByThreshold(y, apex, thr2);
        }
        if (l2 >= 0 && r2 >= 0 && r2 > l2) {
          final int[] lr2 = trimBoundsUnimodal(y, apex, l2, r2, boundaryThresholdFraction, minPoints);
          if (lr2 != null) { left = lr2[0]; right = lr2[1]; }
        }
      }
      // Secondary-peak splitting: if a clear valley separates two local maxima, split here
      final List<int[]> subSegments = splitBySecondaryPeaks(y, left, right, noiseSigma, minPoints);
      if (subSegments.size() <= 1) {
        // scale-aware min width and prominence check to avoid spurious narrow/noisy peaks
        final int minWidth = Math.max(minPoints, Math.max(1, (int)Math.round(1.2 * Math.max(1,charScale))));
        final double base = localBaseline(y, left, right);
        final double prominence = Math.max(0d, y[apex] - base);
        // slope gate: require sufficient edge slope unless CWT SNR is very strong
        final double slope = localMaxSlope(x, y, left, right);
        final double slopeMin = 2.0 * noiseSigma / Math.max(1e-9, medianStep(x));
        if ((right - left + 1) < minWidth || prominence < 2.0 * noiseSigma || slope < slopeMin) {
          // Fallback: accept if CWT coefficient is strong enough
          final double coeffSNR = Math.abs(cwt.get((double) charScale)[apex]) / Math.max(1e-9, noiseSigma);
          if (!(coeffSNR >= 3.5 && (right - left + 1) >= minPoints)) continue;
        }
        boolean ok = true; for (int i = left; i <= right; i++) if (used[i]) { ok = false; break; }
        if (!ok) continue; for (int i = left; i <= right; i++) used[i] = true; ranges.add(Range.closed(x[left], x[right]));
      } else {
        for (int[] seg : subSegments) {
          int la = seg[0], rb = seg[1]; if (rb <= la) continue;
          // recompute around sub-apex
          final int ap = argmax(y, la, rb);
          int cs2 = selectCharacteristicScale(cwt, ap);
          int l3 = findZeroCrossingLeft(cwt.get((double) cs2), ap);
          int r3 = findZeroCrossingRight(cwt.get((double) cs2), ap, y.length);
          if (l3 < 0 || r3 < 0 || r3 <= l3) { final double thr3 = y[ap] * boundaryThresholdFraction; l3 = walkLeftByThreshold(y, ap, thr3); r3 = walkRightByThreshold(y, ap, thr3); }
          if (l3 < 0 || r3 < 0 || r3 <= l3) continue;
          final int[] lr3 = trimBoundsUnimodal(y, ap, Math.max(la,l3), Math.min(rb,r3), boundaryThresholdFraction, minPoints);
          if (lr3 == null) continue; la = lr3[0]; rb = lr3[1];
          final int minWidth = Math.max(minPoints, Math.max(1, (int)Math.round(1.2 * Math.max(1,cs2))));
          final double base = localBaseline(y, la, rb);
          final double prominence = Math.max(0d, y[ap] - base);
          final double slope = localMaxSlope(x, y, la, rb);
          final double slopeMin = 2.0 * noiseSigma / Math.max(1e-9, medianStep(x));
          if ((rb - la + 1) < minWidth || prominence < 2.0 * noiseSigma || slope < slopeMin) {
            final double coeffSNR = Math.abs(cwt.get((double) cs2)[ap]) / Math.max(1e-9, noiseSigma);
            if (!(coeffSNR >= 3.5 && (rb - la + 1) >= minPoints)) continue;
          }
          boolean ok = true; for (int i = la; i <= rb; i++) if (used[i]) { ok = false; break; }
          if (!ok) continue; for (int i = la; i <= rb; i++) used[i] = true; ranges.add(Range.closed(x[la], x[rb]));
        }
      }
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
    }

    // Merge connected/proximal ranges conservatively (also merge small gaps up to ~1 median step)
    ranges.sort(Comparator.comparing(Range::lowerEndpoint));
    final LinkedList<Range<Double>> merged = new LinkedList<>();
    // median RT step
    double step = 0d;
    if (x.length > 1) {
      final double[] dx = new double[x.length - 1];
      for (int i = 1; i < x.length; i++) dx[i - 1] = x[i] - x[i - 1];
      Arrays.sort(dx);
      step = dx[dx.length / 2];
    }
    final double gapTol = Math.max(step * 0.75, 1e-12);
    for (Range<Double> r : ranges) {
      if (merged.isEmpty()) { merged.add(r); continue; }
      Range<Double> last = merged.getLast();
      final boolean overlap = last.isConnected(r) && RangeUtils.rangeLength(last.intersection(r)) > 0d;
      final boolean smallGap = r.lowerEndpoint() - last.upperEndpoint() <= gapTol;
      if (overlap || smallGap) {
        merged.removeLast();
        merged.add(Range.closed(last.lowerEndpoint(),
            Double.max(last.upperEndpoint(), r.upperEndpoint())));
      } else {
        merged.add(r);
      }
    }
    return merged;
  }

  private List<Range<Double>> resolveWithParams(double[] x, double[] y, double[] scales, double snr,
      double minApex, double boundary, int persist,
      CwtRidgeResolverParameters.QuantBaselineMethod qbm, double q) {
    final double[] yDetect = maybeDetrendForDetection(x, y);
    final Map<Double, double[]> cwt = computeCwt(yDetect, scales);
    if (cwt.isEmpty()) return List.of();
    final double[] smallest = cwt.get(scales[0]);
    final double noiseSigma = Math.max(estimateMadSigma(smallest), 1e-9);
    final double thr = snr * noiseSigma;
    final Map<Integer, List<Integer>> maximaByIndexScale = new HashMap<>();
    for (double s : scales) {
      final double[] coeff = cwt.get(s);
      final List<Integer> maxima = findLocalMaxima(coeff, thr);
      for (int idx : maxima) maximaByIndexScale.computeIfAbsent(idx, k -> new ArrayList<>()).add((int) Math.round(s));
    }
    if (maximaByIndexScale.isEmpty()) return List.of();
    final List<Ridge> ridges = traceRidges(cwt, maximaByIndexScale, ridgeNeighborhoodMaxMove);
    final Map<Integer,Integer> apexToLen = new HashMap<>();
    final List<Integer> apexIndicesRaw = new ArrayList<>();
    for (Ridge r : ridges) {
      if (r.length() >= persist) {
        final int i = r.argmaxIndex(cwt);
        if (i >= 0 && y[i] >= minApex) { apexIndicesRaw.add(i); apexToLen.put(i, r.length()); }
      }
    }
    if (apexIndicesRaw.isEmpty()) return List.of();
    // temporarily override baseline behavior
    final var oldQbm = this.quantBaseline; final var oldQ = this.quantBaselineQuantile;
    // build ranges (reuse existing implementation by temporarily setting fields)
    final List<Integer> apexIndices = new ArrayList<>();
    { // simple dedup
      final boolean[] used = new boolean[y.length];
      final List<int[]> scored = new ArrayList<>();
      for (int idx : apexIndicesRaw) {
        double best = 0d; for (double[] coeff : cwt.values()) best = Math.max(best, Math.abs(coeff[idx]));
        scored.add(new int[]{idx, (int) Math.round(best * 1_000_000d)});
      }
      scored.sort((a,b)->Integer.compare(b[1],a[1]));
      final int suppress = Math.max(1, ridgeNeighborhoodMaxMove);
      for (int[] s : scored) {
        int idx = s[0]; boolean conflict=false; int a=Math.max(0,idx-suppress), b=Math.min(y.length-1,idx+suppress);
        for (int i=a;i<=b;i++){ if(used[i]){conflict=true; break;} }
        if(!conflict){ apexIndices.add(idx); for(int i=a;i<=b;i++) used[i]=true; }
      }
    }
    // derive bounds similar to default flow
    final class PeakCand { int apex, left, right, bin; double score; }
    final List<PeakCand> cands = new ArrayList<>();
    for (int apex : apexIndices) {
      int cs = selectCharacteristicScale(cwt, apex);
      int left = findZeroCrossingLeft(cwt.get((double) cs), apex);
      int right = findZeroCrossingRight(cwt.get((double) cs), apex, y.length);
      if (left < 0 || right < 0 || right <= left) {
        final double thrL = y[apex] * boundary; left = walkLeftByThreshold(y, apex, thrL); right = walkRightByThreshold(y, apex, thrL);
      }
      if (left < 0 || right < 0 || right <= left) continue;
      // refine to contiguous-above-threshold, unimodal support around the apex
      final int[] lr = trimBoundsUnimodal(y, apex, left, right, boundary, minPoints);
      if (lr == null) continue;
      left = lr[0]; right = lr[1];

      // Apex refinement in auto path as well
      final int apexRef = argmax(y, left, right);
      if (apexRef != apex) {
        apex = apexRef;
        cs = selectCharacteristicScale(cwt, apex);
        int l2 = findZeroCrossingLeft(cwt.get((double) cs), apex);
        int r2 = findZeroCrossingRight(cwt.get((double) cs), apex, y.length);
        if (l2 < 0 || r2 < 0 || r2 <= l2) {
          final double thr2 = y[apex] * boundary;
          l2 = walkLeftByThreshold(y, apex, thr2);
          r2 = walkRightByThreshold(y, apex, thr2);
        }
        if (l2 >= 0 && r2 >= 0 && r2 > l2) {
          final int[] lr2 = trimBoundsUnimodal(y, apex, l2, r2, boundary, minPoints);
          if (lr2 != null) { left = lr2[0]; right = lr2[1]; }
        }
      }
      // Secondary-peak split inside candidate window for auto path
      final List<int[]> subSegments = splitBySecondaryPeaks(y, left, right, noiseSigma, minPoints);
      if (subSegments.size() <= 1) {
        // scale-aware min width and prominence gate
        final int minWidth = Math.max(minPoints, Math.max(1, (int)Math.round(1.2 * Math.max(1,cs))));
        final double base = localBaseline(y, left, right);
        final double prominence = Math.max(0d, y[apex] - base);
        if ((right - left + 1) < minWidth || prominence < 2.0 * noiseSigma) {
          final double coeffSNR = Math.abs(cwt.get((double) cs)[apex]) / Math.max(1e-9, noiseSigma);
          if (!(coeffSNR >= 3.5 && (right - left + 1) >= minPoints)) continue;
        }
        final double coeffApex = Math.abs(cwt.get((double) cs)[apex]) / noiseSigma; final int len = apexToLen.getOrDefault(apex, 1);
        final PeakCand pc = new PeakCand(); pc.apex = apex; pc.left = left; pc.right = right; pc.score = coeffApex + 0.5 * len; pc.bin = (int)Math.floor((x[apex]-x[0]) / 1.0); cands.add(pc);
      } else {
        for (int[] seg : subSegments) {
          int la = seg[0], rb = seg[1]; if (rb <= la) continue; final int ap = argmax(y, la, rb);
          int cs2 = selectCharacteristicScale(cwt, ap);
          int l3 = findZeroCrossingLeft(cwt.get((double) cs2), ap);
          int r3 = findZeroCrossingRight(cwt.get((double) cs2), ap, y.length);
          if (l3 < 0 || r3 < 0 || r3 <= l3) { final double thr3 = y[ap] * boundary; l3 = walkLeftByThreshold(y, ap, thr3); r3 = walkRightByThreshold(y, ap, thr3); }
          if (l3 < 0 || r3 < 0 || r3 <= l3) continue;
          final int[] lr3 = trimBoundsUnimodal(y, ap, Math.max(la,l3), Math.min(rb,r3), boundary, minPoints);
          if (lr3 == null) continue; la = lr3[0]; rb = lr3[1];
          final int minWidth = Math.max(minPoints, Math.max(1, (int)Math.round(1.2 * Math.max(1,cs2))));
          final double base = localBaseline(y, la, rb);
          final double prominence = Math.max(0d, y[ap] - base);
          if ((rb - la + 1) < minWidth || prominence < 2.0 * noiseSigma) {
            final double coeffSNR = Math.abs(cwt.get((double) cs2)[ap]) / Math.max(1e-9, noiseSigma);
            if (!(coeffSNR >= 3.5 && (rb - la + 1) >= minPoints)) continue;
          }
          final double coeffA = Math.abs(cwt.get((double) cs2)[ap]) / noiseSigma; final int lenA = apexToLen.getOrDefault(ap, 1);
          final PeakCand pc = new PeakCand(); pc.apex = ap; pc.left = la; pc.right = rb; pc.score = coeffA + 0.5 * lenA; pc.bin = (int)Math.floor((x[ap]-x[0]) / 1.0); cands.add(pc);
        }
      }
    }
    if (cands.isEmpty()) return List.of();

    // Per-minute top-K and min-separation to curb FPs (auto only)
    final int topK = 8;
    // adapt min separation to sampling: at least ~3 median steps, but not below 0.03 min
    double dt = 0d; for (int i=1;i<x.length;i++) dt += (x[i]-x[i-1]); dt = (x.length>1)? dt/(x.length-1) : 0.01;
    final double minSep = Math.max(0.03, 3.0 * dt);
    // group by bin
    final Map<Integer,List<PeakCand>> byBin = new HashMap<>();
    for (PeakCand pc : cands) byBin.computeIfAbsent(pc.bin, k->new ArrayList<>()).add(pc);
    final List<PeakCand> filtered = new ArrayList<>();
    for (Map.Entry<Integer,List<PeakCand>> e : byBin.entrySet()) {
      final List<PeakCand> binList = e.getValue();
      binList.sort((a,b)->Double.compare(b.score,a.score));
      // adapt K to density: allow up to min(topK, 1 + binList.size()/3)
      final int maxK = Math.min(topK, 1 + binList.size()/3);
      int kept = 0; for (PeakCand pc : binList) { if (kept >= maxK) break; filtered.add(pc); kept++; }
    }
    // global separation
    filtered.sort((a,b)->Double.compare(b.score,a.score));
    final List<PeakCand> kept = new ArrayList<>();
    for (PeakCand pc : filtered) {
      boolean ok = true; for (PeakCand k : kept) { if (Math.abs(x[pc.apex]-x[k.apex]) < minSep) { ok = false; break; } }
      if (ok) kept.add(pc);
    }

    kept.sort((a,b)->Double.compare(a.left,b.left));
    final List<Range<Double>> ranges = new ArrayList<>();
    for (PeakCand pc : kept) ranges.add(Range.closed(x[pc.left], x[pc.right]));
    return ranges;
  }

  private static Auto tuneFromSeries(double[] x, double[] y) {
    final int n = x.length;
    final double[] dx = new double[Math.max(1,n-1)]; for(int i=1;i<n;i++) dx[i-1]=Math.max(1e-12,x[i]-x[i-1]); Arrays.sort(dx); final double dt=dx[dx.length/2];
    final int w=Math.max(3,Math.min(9,n/200)); final double[] smooth=new double[n]; final int half=w/2; for(int i=0;i<n;i++){int a=Math.max(0,i-half),b=Math.min(n-1,i+half); double s=0; int c=0; for(int j=a;j<=b;j++){s+=y[j];c++;} smooth[i]=c>0?s/c:0;}
    double[] widths=new double[Math.max(1,n/200+1)]; int wc=0; for(int i=1;i<n-1;i++){ if(smooth[i]>smooth[i-1]&&smooth[i]>=smooth[i+1]&&smooth[i]>0){ double hm=smooth[i]*0.5; int l=i; while(l>0&&smooth[l]>=hm&&smooth[l]<=smooth[Math.min(i,l+1)]) l--; int r=i; while(r<n-1&&smooth[r]>=hm&&smooth[r]<=smooth[Math.max(i,r-1)]) r++; widths[wc++]=Math.max(dt,x[Math.min(r,n-1)]-x[Math.max(l,0)]); if(wc>=widths.length) break; } }
    if(wc==0){ widths[wc++]=dt*4;} widths=Arrays.copyOf(widths,wc); Arrays.sort(widths); final double wMed=widths[wc/2];
    final double sMed=Math.max(1,Math.round(wMed/Math.max(1e-12,2.5*dt))); final double[] scales=new double[]{Math.max(1,sMed-2),Math.max(1,sMed-1),sMed,sMed+1,sMed+2,sMed+4,sMed+6};
    final double[] kernel= CwtRidgeResolver.WaveletKernels.ricker((int)Math.max(3,Math.round(10*scales[0])),scales[0]);
    final double[] coeff=convolve(y,kernel); final double noiseSigma=Math.max(estimateMadSigma(coeff),1e-9);
    int maxima=0; double sumSnr=0; for(int i=1;i<coeff.length-1;i++){ if(coeff[i]>coeff[i-1]&&coeff[i]>=coeff[i+1]&&coeff[i]>0){ maxima++; sumSnr+=coeff[i]/noiseSigma; }}
    final double peaksPerMin=maxima/Math.max(1e-9,(x[n-1]-x[0])); double snr=(sumSnr>0&&maxima>0)?Math.min(3.0,Math.max(0.8,(sumSnr/maxima)*0.6)):1.0; if(peaksPerMin>30) snr+=0.5; if(peaksPerMin<5) snr=Math.max(0.8,snr-0.3);
    int persist=sMed>=6?2:1; if(peaksPerMin>30) persist=Math.max(2,persist);
    final double boundary=Math.max(0.03,Math.min(0.08,0.05+(2-Math.min(2,sMed/6.0))*0.01));
    final double p95=io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.CwtRidgeResolver.percentile(y,0.95); final double minApex=Math.max(0.0,Math.min(p95,noiseSigma*3.0));
    final double p10=io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.CwtRidgeResolver.percentile(y,0.10), p50=io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.wavelet.CwtRidgeResolver.percentile(y,0.50);
    final boolean useLowQ=(p50>0&&(p10/p50)>0.2); final double q=0.1;
    return new Auto(scales,snr,persist,boundary,minApex,useLowQ,q);
  }

  private record Auto(double[] scales, double snrCwt, int ridgeMinPersistence, double boundaryFraction,
                      double minApexHeight, boolean useLowQuantileBaseline, double quantileQ){}

  private static Map<Double, double[]> computeCwt(double[] y, double[] scales) {
    final Map<Double, double[]> cwt = new TreeMap<>();
    for (double s : scales) {
      final double[] kernel = WaveletKernels.ricker((int) Math.max(3, Math.round(10 * s)), s);
      cwt.put(s, convolve(y, kernel));
    }
    return cwt;
  }

  private static double[] convolve(double[] signal, double[] kernel) {
    final int n = signal.length, m = kernel.length, half = m / 2;
    final double[] out = new double[n];
    for (int i = 0; i < n; i++) {
      double sum = 0d;
      for (int j = 0; j < m; j++) {
        int idx = i - half + j;
        if (idx >= 0 && idx < n) sum += signal[idx] * kernel[j];
      }
      out[i] = sum;
    }
    return out;
  }

  private static double estimateMadSigma(double[] data) {
    final Median med = new Median();
    final double m = med.evaluate(data);
    final double[] dev = new double[data.length];
    for (int i = 0; i < data.length; i++) dev[i] = Math.abs(data[i] - m);
    return 1.4826 * med.evaluate(dev);
  }

  private static List<Integer> findLocalMaxima(double[] a, double thr) {
    final List<Integer> idx = new ArrayList<>();
    for (int i = 1; i < a.length - 1; i++) {
      if (a[i] > a[i - 1] && a[i] >= a[i + 1] && a[i] > thr) idx.add(i);
    }
    return idx;
  }

  private static class Ridge {
    // map scaleIndex->position index
    final List<int[]> path = new ArrayList<>();
    int length() { return path.size(); }
    int argmaxIndex(Map<Double, double[]> cwt) {
      double best = -Double.MAX_VALUE; int bestPos = -1;
      for (int[] sp : path) {
        final double scale = sp[0];
        final int pos = sp[1];
        final double v = Math.abs(Objects.requireNonNull(cwt.get((double) scale))[pos]);
        if (v > best) { best = v; bestPos = pos; }
      }
      return bestPos;
    }
  }

  private static List<Ridge> traceRidges(Map<Double, double[]> cwt,
      Map<Integer, List<Integer>> maximaByIndexScale, int maxMove) {
    final double[] scaleValues = cwt.keySet().stream().mapToDouble(Double::doubleValue).toArray();
    final List<Ridge> ridges = new ArrayList<>();
    // seed ridges at maxima of smallest scale
    final double smallest = scaleValues[0];
    final List<Integer> seeds = new ArrayList<>();
    for (Map.Entry<Integer, List<Integer>> e : maximaByIndexScale.entrySet()) {
      if (e.getValue().contains((int) Math.round(smallest))) seeds.add(e.getKey());
    }
    seeds.sort(Integer::compareTo);

    for (int seed : seeds) {
      final Ridge r = new Ridge();
      int prevPos = seed;
      r.path.add(new int[]{(int) Math.round(smallest), seed});
      for (int si = 1; si < scaleValues.length; si++) {
        final int scl = (int) Math.round(scaleValues[si]);
        int best = -1; double bestVal = -Double.MAX_VALUE;
        for (int move = -maxMove; move <= maxMove; move++) {
          int pos = prevPos + move;
          final List<Integer> mx = maximaByIndexScale.get(pos);
          if (mx != null && mx.contains(scl)) {
            final double val = Math.abs(cwt.get(scaleValues[si])[pos]);
            if (val > bestVal) { bestVal = val; best = pos; }
          }
        }
        if (best == -1) break; // ridge ends
        r.path.add(new int[]{scl, best});
        prevPos = best;
      }
      ridges.add(r);
    }
    return ridges;
  }

  private static int selectCharacteristicScale(Map<Double, double[]> cwt, int index) {
    double best = -Double.MAX_VALUE; int bestScale = -1;
    for (Map.Entry<Double, double[]> e : cwt.entrySet()) {
      final double v = Math.abs(e.getValue()[index]);
      if (v > best) { best = v; bestScale = (int) Math.round(e.getKey()); }
    }
    return bestScale;
  }

  private static int findZeroCrossingLeft(double[] coeff, int center) {
    for (int i = center - 1; i >= 1; i--) {
      if (coeff[i] == 0d || Math.signum(coeff[i]) != Math.signum(coeff[i - 1])) return i;
    }
    return -1;
  }

  private static int findZeroCrossingRight(double[] coeff, int center, int n) {
    for (int i = center + 1; i < n - 1; i++) {
      if (coeff[i] == 0d || Math.signum(coeff[i]) != Math.signum(coeff[i + 1])) return i;
    }
    return -1;
  }

  private static int walkLeftByThreshold(double[] y, int center, double thr) {
    int i = center; while (i > 0 && y[i - 1] >= thr && y[i - 1] <= y[i]) i--; return i;
  }

  private static int walkRightByThreshold(double[] y, int center, double thr) {
    int i = center; while (i < y.length - 1 && y[i + 1] >= thr && y[i + 1] <= y[i]) i++; return i;
  }

  private double[] maybeDetrendForDetection(double[] x, double[] y) {
    if (!generalParameters.getValue(CwtRidgeResolverParameters.enableDriftAwareDetrend)) {
      return y;
    }
    // detect drift via low-frequency energy ratio using coarse LOESS proxy: compare 5th and 95th percentile windowed means
    final double[] ys = Arrays.copyOf(y, y.length);
    final int w = Math.max(3, y.length / 50);
    final double[] smooth = new double[y.length];
    int half = w / 2;
    for (int i = 0; i < y.length; i++) {
      int a = Math.max(0, i - half), b = Math.min(y.length - 1, i + half);
      double sum = 0d; int cnt = 0;
      for (int j = a; j <= b; j++) { sum += ys[j]; cnt++; }
      smooth[i] = cnt > 0 ? sum / cnt : 0d;
    }
    // heuristic drift metric
    double p5 = percentile(smooth, 0.05), p95 = percentile(smooth, 0.95);
    boolean drift = (p95 - p5) > 0.1 * Math.max(1e-9, percentile(ys, 0.95));
    if (!drift) return y;
    final double[] out = new double[y.length];
    for (int i = 0; i < y.length; i++) out[i] = Math.max(0d, y[i] - smooth[i]);
    return out;
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

  private static final class WaveletKernels {
    static double[] ricker(int size, double scale) {
      final int m = size | 1; // odd
      final int half = m / 2;
      final double[] k = new double[m];
      for (int i = -half, j = 0; j < m; i++, j++) {
        final double t = i / scale;
        k[j] = (1.0 - t * t) * Math.exp(-0.5 * t * t);
      }
      // zero-mean
      double mean = 0; for (double v : k) mean += v; mean /= m; for (int j = 0; j < m; j++) k[j] -= mean;
      return k;
    }
  }

  private static int argmax(double[] y, int left, int right) {
    int a = Math.max(0, left), b = Math.min(y.length - 1, right);
    int idx = a; double best = y[a];
    for (int i = a + 1; i <= b; i++) if (y[i] > best) { best = y[i]; idx = i; }
    return idx;
  }

  private static double medianStep(double[] x){
    if (x.length < 2) return 0.01;
    final double[] dx = new double[x.length - 1];
    for (int i = 1; i < x.length; i++) dx[i - 1] = Math.max(1e-12, x[i] - x[i - 1]);
    Arrays.sort(dx); return dx[dx.length / 2];
  }

  /**
   * Estimate maximum absolute slope within [left,right] using central differences.
   */
  private static double localMaxSlope(double[] x, double[] y, int left, int right){
    int a = Math.max(1, left), b = Math.min(y.length - 2, right);
    double best = 0d;
    for (int i = a; i <= b; i++) {
      final double dx = Math.max(1e-12, x[i+1] - x[i-1]);
      final double dy = y[i+1] - y[i-1];
      best = Math.max(best, Math.abs(dy / dx));
    }
    return best;
  }

  /**
   * Split [left,right] into sub-segments if multiple local maxima with a clear valley exist.
   * Criterion: any valley whose depth exceeds max(1.5*sigma, 0.2*localProminence) triggers a split.
   * Ensures each sub-segment has at least minPoints.
   */
  private static List<int[]> splitBySecondaryPeaks(double[] y, int left, int right, double sigma, int minPoints) {
    final List<int[]> out = new ArrayList<>();
    if (right - left + 1 < Math.max(minPoints * 2, 7)) { out.add(new int[]{left, right}); return out; }
    // find local maxima indices in [left,right]
    final List<Integer> peaks = new ArrayList<>();
    for (int i = Math.max(left + 1, 1); i <= Math.min(right - 1, y.length - 2); i++) {
      if (y[i] > y[i - 1] && y[i] >= y[i + 1]) peaks.add(i);
    }
    if (peaks.size() <= 1) { out.add(new int[]{left, right}); return out; }
    // evaluate valleys between adjacent peaks
    final List<Integer> cuts = new ArrayList<>();
    for (int pi = 0; pi < peaks.size() - 1; pi++) {
      final int p1 = peaks.get(pi), p2 = peaks.get(pi + 1);
      if (p2 - p1 + 1 < minPoints) continue;
      int valley = p1; double minV = y[p1];
      for (int i = p1 + 1; i <= p2; i++) { if (y[i] < minV) { minV = y[i]; valley = i; } }
      final double drop1 = Math.max(0d, y[p1] - minV);
      final double drop2 = Math.max(0d, y[p2] - minV);
      final double thr = Math.max(1.5 * sigma, 0.2 * Math.min(y[p1], y[p2]));
      if (drop1 >= thr && drop2 >= thr) cuts.add(valley);
    }
    if (cuts.isEmpty()) { out.add(new int[]{left, right}); return out; }
    // build segments using cuts, enforcing minPoints
    int segStart = left;
    for (int c : cuts) {
      if (c - segStart + 1 >= minPoints) out.add(new int[]{segStart, c});
      segStart = c + 1;
    }
    if (right - segStart + 1 >= minPoints) out.add(new int[]{segStart, right});
    if (out.isEmpty()) out.add(new int[]{left, right});
    return out;
  }

  /**
   * Compute a low-quantile local baseline within [left,right] to avoid bias from tails.
   */
  private static double localBaseline(double[] y, int left, int right) {
    final int a = Math.max(0, left);
    final int b = Math.min(y.length - 1, right);
    if (a > b) return 0d;
    final double[] w = Arrays.copyOfRange(y, a, b + 1);
    return MathUtils.calcQuantile(w, 0.10);
  }

  /**
   * Trim [left,right] to the maximal contiguous segment around apex that is
   * (1) above a dynamic threshold and (2) unimodal: intensity must monotonically approach the apex
   * from both sides. Ensures at least minPoints samples; returns null if this cannot be achieved.
   */
  private static int[] trimBoundsToContiguousAboveThreshold(double[] y, int apex, int left, int right,
      double thr, int minPoints) {
    int l = apex;
    while (l > left && y[l - 1] >= thr && y[l - 1] <= y[l]) l--;
    int r = apex;
    while (r < right && y[r + 1] >= thr && y[r + 1] <= y[r]) r++;
    // If too short, try minimal symmetric expansion while staying above threshold
    while ((r - l + 1) < minPoints) {
      boolean expanded = false;
      if (l > left && y[l - 1] >= thr) { l--; expanded = true; }
      if ((r - l + 1) < minPoints && r < right && y[r + 1] >= thr) { r++; expanded = true; }
      if (!expanded) break;
    }
    if ((r - l + 1) < minPoints) return null;
    return new int[]{l, r};
  }

  /**
   * High-level wrapper selecting a dynamic threshold using local baseline and apex prominence,
   * then delegating to contiguous + unimodality trimming.
   */
  private static int[] trimBoundsUnimodal(double[] y, int apex, int left, int right,
      double boundaryFraction, int minPoints) {
    final double base = localBaseline(y, left, right);
    final double prominence = Math.max(0d, y[apex] - base);
    final double thrProm = base + 0.10 * prominence; // 10% of prominence above baseline
    final double thrFrac = boundaryFraction * y[apex];
    final double thr = Math.max(thrProm, thrFrac);
    return trimBoundsToContiguousAboveThreshold(y, apex, left, right, thr, minPoints);
  }
}


