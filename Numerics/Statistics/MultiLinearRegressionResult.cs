using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Describes the result of a multiple linear regression fit.
    /// </summary>
    public sealed class MultiLinearRegressionResult : GeneralLinearRegressionResult {

        internal MultiLinearRegressionResult(IReadOnlyList<double> yColumn, IReadOnlyList<IReadOnlyList<double>> xColumns, IReadOnlyList<string> xNames) : base() {

            Debug.Assert(yColumn != null);
            Debug.Assert(xColumns != null);
            Debug.Assert(xColumns.Count > 0);
            Debug.Assert(xNames.Count == xColumns.Count);

            n = yColumn.Count;
            m = xColumns.Count;
            if (n <= m) throw new InsufficientDataException();

            // Compute the design matrix X.
            interceptIndex = -1;
            RectangularMatrix X = new RectangularMatrix(n, m);
            for (int c = 0; c < m; c++) {
                IReadOnlyList<double> xColumn = xColumns[c];
                if (xColumn == null) {
                    Debug.Assert(xNames[c] == "Intercept");
                    Debug.Assert(interceptIndex < 0);
                    for (int r = 0; r < n; r++) {
                        X[r, c] = 1.0;
                    }
                    interceptIndex = c;
                } else {
                    Debug.Assert(xNames[c] != null);
                    if (xColumn.Count != n) throw new DimensionMismatchException();
                    for (int r = 0; r < n; r++) {
                        X[r, c] = xColumn[r];
                    }
                }
            }
            Debug.Assert(interceptIndex >= 0);
            ColumnVector v = new ColumnVector(yColumn);

            // Use X = QR to solve X b = y and compute C.
            QRDecomposition.SolveLinearSystem(X, v, out b, out C);

            // For ANOVA, we will need mean and variance of y
            int yn;
            double ym;
            Univariate.ComputeMomentsUpToSecond(yColumn, out yn, out ym, out SST);

            // Compute residuals
            SSR = 0.0;
            SSF = 0.0;
            ColumnVector yHat = X * b;
            residuals = new List<double>(n);
            for (int i = 0; i < n; i++) {
                double z = yColumn[i] - yHat[i];
                residuals.Add(z);
                SSR += z * z;
                SSF += MoreMath.Sqr(yHat[i] - ym);
            }
            sigma2 = SSR / (n - m);

            // Scale up C by \sigma^2
            // (It sure would be great to be able to overload *=.)
            for (int i = 0; i < m; i++) {
                for (int j = i; j < m; j++) {
                    C[i, j] = C[i, j] * sigma2;
                }
            }

            names = xNames;

        }

        private int n, m;
        private readonly IReadOnlyList<string> names;
        private readonly ColumnVector b;
        private readonly SymmetricMatrix C;
        private readonly double SST, SSF, SSR;
        private readonly double sigma2;
        private readonly List<double> residuals;
        private readonly int interceptIndex;

        /// <inheritdoc/>
        public override UncertainValue Intercept {
            get {
                return (new UncertainValue(b[interceptIndex], Math.Sqrt(C[interceptIndex, interceptIndex])));
            }
        }

        /// <summary>
        /// Gets the coefficient of the input with the given index.
        /// </summary>
        /// <param name="k">The index of the input.</param>
        /// <returns>The best-fit value the coefficient, with uncertainty.</returns>
        public UncertainValue CoefficientOf (int k) {
            return (new UncertainValue(b[k], Math.Sqrt(C[k, k])));
        }

        /// <summary>
        /// Gets the coefficient of the input with the given name.
        /// </summary>
        /// <param name="name">The name of the input column.</param>
        /// <returns>The best-fit value the coefficient, with uncertainty.</returns>
        public UncertainValue CoefficientOf (string name) {
            return (this.Parameters[name].Estimate);
        }

        /// <summary>
        /// Computes the predicted value of y for the given inputs.
        /// </summary>
        /// <param name="x">The inputs, in the same order as the input data columns.</param>
        /// <returns>The predicted value, with uncertainty, of the output.</returns>
        public UncertainValue Predict (params double[] x) {
            return (Predict((IReadOnlyList<double>) x));
        }

        // Add override that accepts parameters

        /// <summary>
        /// Computes the predicted y at a new value of x.
        /// </summary>
        /// <param name="x">The new value of x.</param>
        /// <returns>The value predicted by the model, with uncertainty.</returns>
        public UncertainValue Predict (IReadOnlyList<double> x) {

            if (x == null) throw new ArgumentNullException(nameof(x));
            if (x.Count + 1 != b.Dimension) throw new DimensionMismatchException();

            ColumnVector v = new ColumnVector(x.Count + 1);
            int xIndex = 0;
            for (int vIndex = 0; vIndex < v.Dimension; vIndex++) {
                if (vIndex == interceptIndex) {
                    v[vIndex] = 1.0;
                } else {
                    v[vIndex] = x[xIndex];
                    xIndex++;
                }
            }

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
            AnovaRow fit = new AnovaRow(SSF, m - 1);
            AnovaRow residual = new AnovaRow(SSR, n - m);
            AnovaRow total = new AnovaRow(SST, n - 1);
            OneWayAnovaResult anova = new OneWayAnovaResult(fit, residual, total);
            return (anova);
        }

        internal override ParameterCollection CreateParameters () {
            return (new ParameterCollection(names, b, C));
        }

    }
}
