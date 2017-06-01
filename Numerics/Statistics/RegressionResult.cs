using System;
using System.Diagnostics;
using System.Collections.Generic;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents the result of a regression fit.
    /// </summary>
    public class RegressionResult {

        internal RegressionResult (ParameterCollection parameters, OneWayAnovaResult anova, Sample residuals) {
            Debug.Assert(parameters != null);
            Debug.Assert(anova != null);
            Debug.Assert(residuals != null);
            this.parameters = parameters;
            this.anova = anova;
            this.residuals = residuals;
        }

        private readonly ParameterCollection parameters;

        private readonly OneWayAnovaResult anova;

        private readonly Sample residuals;

        /// <summary>
        /// Gets the parameters of the regression.
        /// </summary>
        public virtual ParameterCollection Parameters {
            get {
                return (parameters);
            }
        }

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

}
