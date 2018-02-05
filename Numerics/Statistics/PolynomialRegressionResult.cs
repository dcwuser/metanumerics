using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents the result of a multiple linear regression fit.
    /// </summary>
    public sealed class MultiLinearRegressionResult : GeneralLinearRegressionResult {

        internal MultiLinearRegressionResult (ParameterCollection parameters, OneWayAnovaResult anova, Sample residuals) : base(parameters, anova, residuals) {

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
        public UncertainValue Predict (IReadOnlyList<double> x) {

            if (x == null) throw new ArgumentNullException(nameof(x));
            if (x.Count != this.Parameters.Best.Dimension) throw new DimensionMismatchException();

            ColumnVector v = new ColumnVector(x);

            double y = v.Transpose() * this.Parameters.Best;

            double sigma2 = this.Anova.Residual.SumOfSquares / this.Anova.Residual.DegreesOfFreedom;

            double vCv = v.Transpose() * this.Parameters.Covariance * v;

            double dy = Math.Sqrt(sigma2 + vCv);

            return (new UncertainValue(y, dy));

        }


    }

    /// <summary>
    /// Represents the result of a polynomial regression fit.
    /// </summary>
    public sealed class PolynomialRegressionResult : GeneralLinearRegressionResult {

        internal PolynomialRegressionResult (ParameterCollection parameters, OneWayAnovaResult anova, Sample residuals) : base(parameters, anova, residuals) {

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

            double y = v.Transpose() * this.Parameters.Best;

            double sigma2 = this.Anova.Residual.SumOfSquares / this.Anova.Residual.DegreesOfFreedom;

            double vCv = v.Transpose() * this.Parameters.Covariance * v;

            double dy = Math.Sqrt(sigma2 + vCv);

            return (new UncertainValue(y, dy));

        }

    }
}
