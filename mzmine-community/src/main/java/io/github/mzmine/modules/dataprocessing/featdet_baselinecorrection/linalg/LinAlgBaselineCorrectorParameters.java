/*
 * Copyright (c) 2004-2025 The MZmine Development Team
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

import io.github.mzmine.main.MZmineCore;
import io.github.mzmine.modules.dataprocessing.featdet_baselinecorrection.AbstractBaselineCorrectorParameters;
import io.github.mzmine.parameters.parametertypes.DoubleParameter;
import io.github.mzmine.parameters.parametertypes.IntegerParameter;

public class LinAlgBaselineCorrectorParameters extends AbstractBaselineCorrectorParameters {

  public static final IntegerParameter polynomialDegree = new IntegerParameter("Polynomial degree",
      """
      Degree of the polynomial that will estimate the data baseline.
      A low degree may fail to detect all the baseline present, while a high
      degree may make the data too oscillatory, especially at the edges.
      Typical values range from 2 to 6.
      """, 3, 1, 10);

  public static final DoubleParameter threshold = new DoubleParameter("Tolerance", """
      Tolerance to use when comparing the difference between the current
      fit coefficients and the ones from the last iteration.
      The iteration procedure will stop when the difference between them is lower than this value.
      Typical values range from 1e-6 to 1e-2.
      """, MZmineCore.getConfiguration().getScoreFormat(), 1e-3);

  public static final IntegerParameter iterations = new IntegerParameter("Max iterations",
      """
      Maximum number of iterations to perform.
      Higher values allow for better convergence but increase computation time.
      Typical values range from 10 to 200.
      """, 100, 1, 500);

  public LinAlgBaselineCorrectorParameters() {
    // Note: applyPeakRemoval is intentionally excluded for LinAlg algorithm
    // because it has built-in iterative peak suppression
    super(samplePercentage.cloneParameter(), polynomialDegree, threshold, iterations);
  }
}