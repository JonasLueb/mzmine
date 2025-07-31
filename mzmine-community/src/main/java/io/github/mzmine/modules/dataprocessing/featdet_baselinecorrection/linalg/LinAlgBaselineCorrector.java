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

package io.github.mzmine.modules.dataprocessing.featdet_baselinecorrection.linalg;

import io.github.mzmine.datamodel.features.FeatureList;
import io.github.mzmine.modules.dataprocessing.featdet_baselinecorrection.BaselineCorrectionParameters;
import io.github.mzmine.modules.dataprocessing.featdet_baselinecorrection.BaselineCorrector;
import io.github.mzmine.modules.dataprocessing.featdet_baselinecorrection.UnivariateBaselineCorrector;
import io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.minimumsearch.MinimumSearchFeatureResolver;
import io.github.mzmine.parameters.ParameterSet;
import io.github.mzmine.util.MemoryMapStorage;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class LinAlgBaselineCorrector extends UnivariateBaselineCorrector {

  private final int polynomialDegree;
  private final double tolerance;
  private final int maxIterations;

  public LinAlgBaselineCorrector() {
    this(null, 50, 3, 1e-3, 100, "baseline", null);
  }

  public LinAlgBaselineCorrector(@Nullable MemoryMapStorage storage, double samplePercentage,
      int polynomialDegree, double tolerance, int maxIterations, @NotNull String suffix,
      @Nullable MinimumSearchFeatureResolver resolver) {
    super(storage, samplePercentage, suffix, resolver);
    this.polynomialDegree = polynomialDegree;
    this.tolerance = tolerance;
    this.maxIterations = maxIterations;
  }

  @Override
  public @NotNull String getName() {
    return "Linear Algebra baseline corrector";
  }

  @Override
  public @NotNull Class<? extends ParameterSet> getParameterSetClass() {
    return LinAlgBaselineCorrectorParameters.class;
  }

  @Override
  public BaselineCorrector newInstance(ParameterSet parameters, MemoryMapStorage storage,
      FeatureList flist) {
    final ParameterSet embedded = parameters.getParameter(
        BaselineCorrectionParameters.correctionAlgorithm).getEmbeddedParameters();
    
    // LinAlg algorithm doesn't use peak removal - it has built-in peak suppression
    final MinimumSearchFeatureResolver resolver = null;

    return new LinAlgBaselineCorrector(storage,
        embedded.getValue(LinAlgBaselineCorrectorParameters.samplePercentage),
        embedded.getValue(LinAlgBaselineCorrectorParameters.polynomialDegree),
        embedded.getValue(LinAlgBaselineCorrectorParameters.threshold),
        embedded.getValue(LinAlgBaselineCorrectorParameters.iterations),
        parameters.getValue(BaselineCorrectionParameters.suffix), resolver);
  }

  @Override
  protected UnivariateFunction initializeFunction(double[] rtValues, double[] intensityValues) {
    
    // Implementation of PeakUtils baseline algorithm
    // Based on: https://bitbucket.org/lucashnegri/peakutils/src/master/peakutils/baseline.py
    
    final int n = intensityValues.length;
    final int order = polynomialDegree + 1;
    
    // Create normalized x values to avoid numerical issues
    final double maxIntensity = Math.abs(getMaxAbsValue(intensityValues));
    final double cond = Math.pow(maxIntensity, 1.0 / order);
    
    // Create normalized x coordinates
    final double[] x = new double[n];
    for (int i = 0; i < n; i++) {
      x[i] = cond * i / (n - 1);
    }
    
    // Create Vandermonde matrix for polynomial fitting
    final RealMatrix vanderMatrix = createVandermondeMatrix(x, order);
    
    // Compute pseudo-inverse using SVD
    final SingularValueDecomposition svd = new SingularValueDecomposition(vanderMatrix);
    final DecompositionSolver solver = svd.getSolver();
    final RealMatrix vanderPinv = solver.getInverse();
    
    // Initialize coefficients and working arrays
    RealVector coeffs = new ArrayRealVector(order, 1.0);
    final double[] y = intensityValues.clone(); // Working copy
    final double[] baseline = new double[n];
    
    // Iterative fitting
    for (int iter = 0; iter < maxIterations; iter++) {
      // Fit polynomial to current y values
      final RealVector yVector = new ArrayRealVector(y);
      final RealVector coeffsNew = vanderPinv.operate(yVector);
      
      // Check convergence
      final RealVector coeffDiff = coeffsNew.subtract(coeffs);
      final double coeffNorm = coeffs.getNorm();
      final double convergence = coeffNorm > 0 ? coeffDiff.getNorm() / coeffNorm : coeffDiff.getNorm();
      
      if (convergence < tolerance) {
        break;
      }
      
      // Update coefficients
      coeffs = coeffsNew;
      
      // Compute new baseline
      final RealVector baselineVector = vanderMatrix.operate(coeffs);
      for (int i = 0; i < n; i++) {
        baseline[i] = baselineVector.getEntry(i);
      }
      
      // Apply minimum operation: y = min(original_y, baseline)
      for (int i = 0; i < n; i++) {
        y[i] = Math.min(intensityValues[i], baseline[i]);
      }
    }
    
    // Create final baseline function
    final double[] finalCoeffs = coeffs.toArray();
    final double rtMin = rtValues[0];
    final double rtMax = rtValues[rtValues.length - 1];
    
    return rtValue -> {
      // Normalize RT value to [0, cond] range
      final double normalizedRt = cond * (rtValue - rtMin) / (rtMax - rtMin);
      
      // Evaluate polynomial
      double result = 0.0;
      double rtPower = 1.0;
      for (int i = order - 1; i >= 0; i--) {
        result += finalCoeffs[i] * rtPower;
        rtPower *= normalizedRt;
      }
      
      return result;
    };
  }
  
  private double getMaxAbsValue(double[] values) {
    double max = 0.0;
    for (double value : values) {
      max = Math.max(max, Math.abs(value));
    }
    return max;
  }
  
  private RealMatrix createVandermondeMatrix(double[] x, int order) {
    final int n = x.length;
    final double[][] matrix = new double[n][order];
    
    for (int i = 0; i < n; i++) {
      double xPower = 1.0;
      // Fill from highest degree to lowest (reverse order for Vandermonde)
      for (int j = order - 1; j >= 0; j--) {
        matrix[i][j] = xPower;
        xPower *= x[i];
      }
    }
    
    return new Array2DRowRealMatrix(matrix);
  }

  public @NotNull String getDescription() {
    return "PeakUtils-style iterative polynomial baseline correction with built-in peak suppression. "
        + "Iteratively fits polynomials while automatically suppressing peaks using min(data, baseline) operation. "
        + "No pre-processing peak removal needed - algorithm handles peaks internally.";
  }
}