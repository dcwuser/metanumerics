using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents the result of a polynomial regression fit.
    /// </summary>
    public sealed class PolynomialRegressionResult : GeneralLinearRegressionResult {

        internal PolynomialRegressionResult (IReadOnlyList<double> x, IReadOnlyList<double> y, int degree) : base() {

            Debug.Assert(x != null);
            Debug.Assert(y != null);
            Debug.Assert(x.Count == y.Count);
            Debug.Assert(degree >= 0);

            m = degree;
            n = x.Count;
            if (n < (m + 1)) throw new InsufficientDataException();

            // Construct the n X m design matrix X_{ij} = x_{i}^{j}
            RectangularMatrix X = new RectangularMatrix(n, m + 1);
            ColumnVector Y = new ColumnVector(n);
            for (int i = 0; i < n; i++) {
                double x_i = x[i];
                X[i, 0] = 1.0;
                for (int j = 1; j <= m; j++) {
                    X[i, j] = X[i, j - 1] * x_i;
                }
                double y_i = y[i];
                Y[i] = y_i;
            }

            // Use X = QR to solve X b = y and compute C
            QRDecomposition.SolveLinearSystem(X, Y, out b, out C);

            // Compute mean and total sum of squares.
            // This could be done inside loop above, but this way we get to re-use code from Univariate.
            double yMean;
            Univariate.ComputeMomentsUpToSecond(y, out n, out yMean, out SST);

            // Compute residuals
            SSR = 0.0;
            SSF = 0.0;
            ColumnVector yHat = X * b;
            residuals = new List<double>(n);
            for (int i = 0; i < n; i++) {
                double z = y[i] - yHat[i];
                residuals.Add(z);
                SSR += z * z;
                SSF += MoreMath.Sqr(yHat[i] - yMean);
            }
            sigma2 = SSR / (n - (m + 1));

            // Scale up C by \sigma^2
            // (It sure would be great to be able to overload *=.)
            for (int i = 0; i <= m; i++) {
                for (int j = i; j <= m; j++) {
                    C[i, j] = C[i, j] * sigma2;
                }
            }

        }

        private int n, m;

        private readonly ColumnVector b;

        private readonly SymmetricMatrix C;

        private readonly double SST, SSF, SSR;

        private readonly double sigma2;

        private readonly List<double> residuals;

        /// <inheritdoc/>
        public override UncertainValue Intercept {
            get {
                return (Coefficient(0));
            }
        }

        /// <summary>
        /// Gets the coefficient of x<sup>k</sup>.
        /// </summary>
        /// <param name="k">The degree of the term, with must lie between 0 and the degree of the model polynomial.</param>
        /// <returns>The best-fit value the coefficient, with uncertainty.</returns>
        /// <remarks>
        /// <para>Note that the <paramref name="k"/>=0 term is the intercept.</para>
        /// </remarks>
        public UncertainValue Coefficient (int k) {
            return (this.Parameters[k].Estimate);
        }

        /// <summary>
        /// Computes the predicted y at a new value of x.
        /// </summary>
        /// <param name="x">The new value of x.</param>
        /// <returns>The value predicted by the model, with uncertainty.</returns>
        public UncertainValue Predict (double x) {

            double[] vStore = new double[this.Parameters.Count];
            double xk = 1.0;
            vStore[0] = xk;
            for (int k = 1; k < vStore.Length; k++) {
                xk *= x;
                vStore[k] = xk;
            }
            ColumnVector v = new ColumnVector(vStore, vStore.Length);

            double y = v.Transpose() * b;

            double vCv = v.Transpose() * C * v;

            double dy = Math.Sqrt(sigma2 + vCv);

            return (new UncertainValue(y, dy));

        }

        /// <summary>
        /// Gets the residuals 
        /// </summary>
        /// <value>A list of the differences between each measured and predicted value.</value>
        public IReadOnlyList<double> Residuals {
            get {
                return (residuals);
            }
        }

        internal override OneWayAnovaResult CreateAnova () {
            // Use sums-of-squares to do ANOVA
            AnovaRow fit = new AnovaRow(SSF, m);
            AnovaRow residual = new AnovaRow(SSR, n - (m + 1));
            AnovaRow total = new AnovaRow(SST, n - 1);
            OneWayAnovaResult anova = new OneWayAnovaResult(fit, residual, total);
            return (anova);
        }

        internal override ParameterCollection CreateParameters () {
            string[] names = new string[m + 1];
            names[0] = "1";
            if (m > 0) names[1] = "x";
            for (int i = 2; i <= m; i++) names[i] = $"x^{i}";
            ParameterCollection parameters = new ParameterCollection(names, b, C);
            return (parameters);
        }

    }
}
