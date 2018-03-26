using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

using Meta.Numerics.Analysis;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Describes the result of a linear logistic regression fit.
    /// </summary>
    public sealed class LinearLogisticRegressionResult : FitResult {

        // We need a goodness-of-fit measurement

        internal LinearLogisticRegressionResult (IReadOnlyList<double> x, IReadOnlyList<bool> y) {

            Debug.Assert(x != null);
            Debug.Assert(y != null);
            Debug.Assert(x.Count == y.Count);

            // check size of data set
            int n = x.Count;
            if (n < 3) throw new InsufficientDataException();

            // The linear logistic model is:
            //   p_i = \sigma(t_i) \quad t_i = a + b x_i
            // So the log likelihood of the data set under the model is:
            //   \ln L = \sum_{{\rm true} i} \ln p_i + \sum_{{\rm false} i} \ln (1 - p_i)
            //         = \sum_{{\rm true} i} \ln \sigma(t_i) + \sum_{{\rm false} i} \ln (1 - \sigma(t_i))
            // Taking derivatives:
            //   \frac{\partial L}{\partial a} = \sum_{{\rm true} i} \frac{\sigma'(t_i)}{\sigma(t_i)}
            //     + \sum_{{\rm false} i} \frac{-\sigma'(t_i)}{1 - \sigma(t_i)}
            //   \frac{\partial L}{\partial b} = \sum_{{\rm true} i} \frac{\sigma'(t_i)}{\sigma(t_i)} x_i
            //     + \sum_{{\rm false} i} \frac{-\sigma'(t_i)}{1 - \sigma(t_i)} x_i
            // Using \sigma(t) = \frac{1}{1 + e^{-t}}, we can derive:
            //   \frac{\sigma'(t)}{\sigma(t)} = \sigma(-t)
            //   \frac{\sigma'(t)}{1 - \sigma(t)} = \sigma(t)
            // So this becomes
            //   \frac{\partial L}{\partial a} = \sum_i \pm \sigma(\mp t_i)
            //   \frac{\partial L}{\partial b} = \sum_i \pm \sigma(\mp t_i) x_i
            // where the upper sign is for true values and the lower sign is for false values.
            // Find the simultaneous zeros of these equations to obtain the likelihood-maximizing a, b.

            // To get the curvature matrix, we need the second derivatives.
            //   \frac{\partial^2 L}{\partial a^2} = - \sum_i \sigma'(\mp t_i)
            //   \frac{\partial^2 L}{\partial a \partial b} = - \sum_i \sigma'(\mp t_i) x_i
            //   \frac{\partial^2 L}{\partial b^2} = - \sum_i \sigma'(\mp t_i) x_i^2

            // We need an initial guess at the parameters. Begin with the Ansatz of the logistic model:
            //    \frac{p}{1-p} = e^{\alpha + \beta x}
            // Differentiate and do some algebra to get:
            //    \frac{\partial p}{\partial x} = \beta p ( 1 - p)
            // Evaluating at means, and noting that p (1 - p) = var(y) and that, in a development around the means,
            //    cov(p, x) = \frac{\partial p}{\partial x} var(x)
            // we get
            //    \beta = \frac{cov(y, x)}{var(x) var(y)}
            // This approximation gets the sign right, but it looks like it usually gets the magnitude quite wrong.
            // The problem with the approach is that var(y) = p (1 - p) assumes y are chosen with fixed p, but they aren't.
            // We need to re-visit this analysis.

            double xMean, yMean, xxSum, yySum, xySum;
            Bivariate.ComputeBivariateMomentsUpToTwo(x, y.Select(z => z ? 1.0 : 0.0), out n, out xMean, out yMean, out xxSum, out yySum, out xySum);
            double p = yMean;
            double b0 = xySum / xxSum / yySum * n;
            double a0 = Math.Log(p / (1.0 - p)) - b0 * xMean;

            Func<IReadOnlyList<double>, IReadOnlyList<double>> J = (IReadOnlyList<double> a) => {
                double dLda = 0.0;
                double dLdb = 0.0;
                for (int i = 0; i < n; i++) {
                    double t = a[0] + a[1] * x[i];
                    if (y[i]) {
                        double s = Sigma(-t);
                        dLda += s;
                        dLdb += s * x[i];
                    } else {
                        double s = Sigma(t);
                        dLda -= s;
                        dLdb -= s * x[i];
                    }
                }
                return (new double[] { dLda, dLdb });
            };

            ColumnVector b = MultiFunctionMath.FindZero(J, new double[] { a0, b0 });

            SymmetricMatrix C = new SymmetricMatrix(2);
            for (int i = 0; i < n; i++) {
                double t = b[0] + b[1] * x[i];
                if (y[i]) t = -t;
                double e = Math.Exp(-t);
                double sp = e / MoreMath.Sqr(1.0 + e);
                C[0, 0] += sp;
                C[0, 1] += sp * x[i];
                C[1, 1] += sp * x[i] * x[i];
            }
            CholeskyDecomposition CD = C.CholeskyDecomposition();
            if (CD == null) throw new DivideByZeroException();
            C = CD.Inverse();

            best = b;
            covariance = C;

        }

        private static double Sigma (double t) {
            return (1.0 / (1.0 + Math.Exp(-t)));
        }

        private ColumnVector best;
        private SymmetricMatrix covariance;

        internal override ParameterCollection CreateParameters () {
            return (new ParameterCollection(new string[] { nameof(Intercept), nameof(Slope) }, best, covariance));
        }

        /// <summary>
        /// Gets the intercept of the regression model.
        /// </summary>
        public UncertainValue Intercept {
            get {
                return (new UncertainValue(best[0], Math.Sqrt(covariance[0, 0])));
            }
        }

        /// <summary>
        /// Gets the slope of the regression model.
        /// </summary>
        public UncertainValue Slope {
            get {
                return (new UncertainValue(best[1], Math.Sqrt(covariance[1, 1])));
            }
        }

        /// <summary>
        /// Computes the outcome probability for the given value.
        /// </summary>
        /// <param name="x">The value of the independent variable.</param>
        /// <returns>The probability of a <see langword="true"/> outcome, as predicted by the model.</returns>
        public UncertainValue Predict (double x) {

            double t = best[0] + best[1] * x;
            double vt = covariance[0, 0] + 2.0 * covariance[0, 1] * x + covariance[1, 1] * x * x;

            double e = Math.Exp(-t);
            double p = 1.0 / (1.0 + e);
            double dp = e * p * p * Math.Sqrt(vt);

            return (new UncertainValue(p, dp));
        }

    }

}
