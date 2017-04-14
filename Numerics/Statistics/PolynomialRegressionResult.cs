using System;
using System.Collections.Generic;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {


    /// <summary>
    /// Represents the result of a regression fit.
    /// </summary>
    public class RegressionResult : FitResult {

        internal RegressionResult (ColumnVector b, SymmetricMatrix C, OneWayAnovaResult anova, Sample residuals) : base(b, C, anova.Result) {
            this.anova = anova;
            this.residuals = residuals;
            //residuals.IsReadOnly = true;
        }

        private OneWayAnovaResult anova;

        private Sample residuals;

        /// <summary>
        /// Gets r<sup>2</sup> for the regression. 
        /// </summary>
        public virtual double RSquared {
            get {
                return (this.Anova.Factor.SumOfSquares / this.Anova.Total.SumOfSquares);
            }
        }


        /// <summary>
        /// Gets an F-test for the regression.
        /// </summary>
        public virtual TestResult F {
            get {
                return (this.Anova.Factor.Result);
            }
        }

        /// <summary>
        /// Gets an ANOVA analysis of the regression.
        /// </summary>
        public virtual OneWayAnovaResult Anova {
            get {
                return (anova);
            }
        }

        /// <summary>
        /// Gets the residuals.
        /// </summary>
        public virtual Sample Residuals {
            get {
                return (residuals);
            }
        }

    }

    /// <summary>
    /// Represents the result of a multiple linear regression fit.
    /// </summary>
    public class MultiLinearRegressionResult : RegressionResult {

        internal MultiLinearRegressionResult (ColumnVector b, SymmetricMatrix C, OneWayAnovaResult anova, Sample residuals) : base(b, C, anova, residuals) {

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
            return (this.Parameter(k));
        }


        /// <summary>
        /// Computes the predicted y at a new value of x.
        /// </summary>
        /// <param name="x">The new value of x.</param>
        /// <returns>The value predicted by the model, with uncertainty.</returns>
        public UncertainValue Predict (IList<double> x) {

            ColumnVector v = new ColumnVector(x);

            double y = v.Transpose() * this.Parameters;

            double sigma2 = this.Anova.Residual.SumOfSquares / this.Anova.Residual.DegreesOfFreedom;

            double vCv = v.Transpose() * this.CovarianceMatrix * v;

            double dy = Math.Sqrt(sigma2 + vCv);

            return (new UncertainValue(y, dy));

        }


    }

    /// <summary>
    /// Represents the result of a polynomial regression fit.
    /// </summary>
    public class PolynomialRegressionResult : RegressionResult {

        internal PolynomialRegressionResult (ColumnVector b, SymmetricMatrix C, OneWayAnovaResult anova, Sample residuals) : base(b, C, anova, residuals) {

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
            return (this.Parameter(k));
        }


        /// <summary>
        /// Computes the predicted y at a new value of x.
        /// </summary>
        /// <param name="x">The new value of x.</param>
        /// <returns>The value predicted by the model, with uncertainty.</returns>
        public UncertainValue Predict (double x) {

            double[] vStore = new double[this.Dimension];
            double xk = 1.0;
            vStore[0] = xk;
            for (int k = 1; k < vStore.Length; k++) {
                xk *= x;
                vStore[k] = xk;
            }
            ColumnVector v = new ColumnVector(vStore, vStore.Length);

            double y = v.Transpose() * this.Parameters;

            double sigma2 = this.Anova.Residual.SumOfSquares / this.Anova.Residual.DegreesOfFreedom;

            double vCv = v.Transpose() * this.CovarianceMatrix * v;

            double dy = Math.Sqrt(sigma2 + vCv);

            return (new UncertainValue(y, dy));

        }

    }
}
