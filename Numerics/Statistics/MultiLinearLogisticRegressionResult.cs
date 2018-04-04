using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Analysis;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Describes the result of a linear logistic regression fit.
    /// </summary>
    public sealed class MultiLinearLogisticRegressionResult : FitResult {

        internal MultiLinearLogisticRegressionResult (IReadOnlyList<bool> yColumn, IReadOnlyList<IReadOnlyList<double>> xColumns, IReadOnlyList<string> xNames) {

            Debug.Assert(yColumn != null);
            Debug.Assert(xColumns != null);
            Debug.Assert(xNames != null);
            Debug.Assert(xColumns.Count == xNames.Count);

            int n = yColumn.Count;
            int m = xColumns.Count;
            if (n <= m) throw new InsufficientDataException();

            interceptIndex = -1;
            for (int c = 0; c < m; c++) {
                IReadOnlyList<double> xColumn = xColumns[c];
                if (xColumn == null) {
                    Debug.Assert(interceptIndex < 0);
                    Debug.Assert(xNames[c] == "Intercept");
                    interceptIndex = c;
                }
                else
                {
                    if (xColumn.Count != n) throw new DimensionMismatchException();
                }
            }
            Debug.Assert(interceptIndex >= 0);


            // Define the log likelihood as a function of the parameter set
            Func<IReadOnlyList<double>, double> logLikelihood = (IReadOnlyList<double> a) => {
                Debug.Assert(a != null);
                Debug.Assert(a.Count == m);

                double L = 0.0;
                for (int k = 0; k < n; k++) {
                    double t = 0.0;
                    for (int i = 0; i < m; i++) {
                        if (i == interceptIndex) {
                            t += a[i];
                        } else {
                            t += a[i] * xColumns[i][k];
                        }
                    }
                    double ez = Math.Exp(t);

                    if (yColumn[k]) {
                        L -= MoreMath.LogOnePlus(1.0 / ez);
                    } else {
                        L -= MoreMath.LogOnePlus(ez);
                    }
                }
                return (L);
            };

            // We need  a better starting value.
            double[] start = new double[m];
            //double[] start = new double[] { -1.5, +2.5, +0.5 };

            // Search out the likelihood-maximizing parameter set.
            MultiExtremum maximum = MultiFunctionMath.FindLocalMaximum(logLikelihood, start);
            b = maximum.Location;
            CholeskyDecomposition CD = maximum.HessianMatrix.CholeskyDecomposition();
            if (CD == null) throw new DivideByZeroException();
            C = CD.Inverse();

            names = xNames;

        }

        private readonly int interceptIndex;
        private readonly ColumnVector b;
        private readonly SymmetricMatrix C;
        private readonly IReadOnlyList<string> names;

        /// <summary>
        /// Gets the intercept of the regression model.
        /// </summary>
        public UncertainValue Intercept {
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
        /// Computes the predicted outcome for the given inputs.
        /// </summary>
        /// <param name="x">The inputs, in the same order as the input data columns.</param>
        /// <returns>The probability of a <see langword="true"/> outcome, as predicted by the model.</returns>
        public UncertainValue Predict (params double[] x) {
            return (Predict((IReadOnlyList<double>) x));
        }

        /// <summary>
        /// Computes the predicted outcome for the given inputs.
        /// </summary>
        /// <param name="x">The value of the independent variable.</param>
        /// <returns>The probability of a <see langword="true"/> outcome, as predicted by the model.</returns>
        public UncertainValue Predict (IReadOnlyList<double> x) {
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (x.Count + 1 != b.Dimension) throw new DimensionMismatchException();

            ColumnVector v = new ColumnVector(b.Dimension);
            int xIndex = 0;
            for (int i = 0; i < v.Dimension; i++) {
                if (i == interceptIndex) {
                    v[i] = 1.0;
                } else {
                    v[i] = x[xIndex];
                    xIndex++;
                }
            }

            double t = v.Transpose() * b;
            double vt = v.Transpose() * C * v;

            double e = Math.Exp(-t);
            double p = 1.0 / (1.0 + e);
            double dp = e * p * p * Math.Sqrt(vt);

            return (new UncertainValue(p, dp));
            //return (1.0 / (1.0 + Math.Exp(-t)));
        }

        internal override ParameterCollection CreateParameters () {
            return (new ParameterCollection(names, b, C));
        }

    }

}
