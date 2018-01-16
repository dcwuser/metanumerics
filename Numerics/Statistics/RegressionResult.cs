using System;
using System.Diagnostics;
using System.Collections.Generic;

namespace Meta.Numerics.Statistics {



    public class RegressionResult {

        internal RegressionResult (ParameterCollection parameters, Sample residuals) {
            Debug.Assert(parameters != null);
            Debug.Assert(residuals != null);
            this.parameters = parameters;
            this.residuals = residuals;
        }

        private readonly ParameterCollection parameters;

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
        /// Gets the residuals.
        /// </summary>
        public virtual Sample Residuals {
            get {
                return (residuals);
            }
        }

    }

    /// <summary>
    /// Represents the result of a general linear regression.
    /// </summary>
    public class GeneralLinearRegressionResult : RegressionResult {

        internal GeneralLinearRegressionResult (ParameterCollection parameters, OneWayAnovaResult anova, Sample residuals) : base(parameters, residuals) {
            Debug.Assert(anova != null);
            this.anova = anova;
        }

        private readonly OneWayAnovaResult anova;

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

    }

}
