using System;
using System.Collections.Generic;
using System.Diagnostics;


using Meta.Numerics.Analysis;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Describes the result of a fit to a non-linear function.
    /// </summary>
    public sealed class NonlinearRegressionResult : FitResult {

        internal NonlinearRegressionResult(
            IReadOnlyList<double> x, IReadOnlyList<double> y,
            Func<IReadOnlyList<double>, double, double> function,
            IReadOnlyList<double> start, IReadOnlyList<string> names) {

            Debug.Assert(x != null);
            Debug.Assert(y != null);
            Debug.Assert(function != null);
            Debug.Assert(start != null);
            Debug.Assert(names != null);
            Debug.Assert(x.Count == y.Count);
            Debug.Assert(start.Count > 0);
            Debug.Assert(names.Count == start.Count);

            int n = x.Count;
            int d = start.Count;
            if (n <= d) throw new InsufficientDataException();

            MultiExtremum min = MultiFunctionMath.FindLocalMinimum((IReadOnlyList<double> a) => {
                double ss = 0.0;
                for (int i = 0; i < n; i++) {
                    double r = y[i] - function(a, x[i]);
                    ss += r * r;
                }
                return (ss);
            }, start);

            CholeskyDecomposition cholesky = min.HessianMatrix.CholeskyDecomposition();
            if (cholesky == null) throw new DivideByZeroException();

            b = min.Location;
            C = cholesky.Inverse();
            C = (2.0 * min.Value / (n - d)) * C;

            residuals = new List<double>(n);
            for (int i = 0; i < n; i++) {
                double r = y[i] - function(b, x[i]);
                residuals.Add(r);
            }

            this.names = names;
            this.function = function;
        }

        private readonly IReadOnlyList<string> names;
        private readonly ColumnVector b;
        private readonly SymmetricMatrix C;
        private readonly Func<IReadOnlyList<double>, double, double> function;
        private readonly List<double> residuals;

        /// <summary>
        /// Gets the residuals from the fit.
        /// </summary>
        public IReadOnlyList<double> Residuals {
            get {
                return (residuals);
            }
        }

        /// <summary>
        /// Predicts the value for a new x based on the fitted model.
        /// </summary>
        /// <param name="x">The value of the independent variable.</param>
        /// <returns>The predicted value of the dependent variable.</returns>
        public double Predict (double x) {
            return (function(b, x));
        }

        internal override ParameterCollection CreateParameters () {
            return (new ParameterCollection(names, b, C));
        }

    }
}
