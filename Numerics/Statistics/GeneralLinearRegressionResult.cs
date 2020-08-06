using System;
using System.Diagnostics;
using System.Collections.Generic;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Describes the result of any generalized linear regression.
    /// </summary>
    public abstract class GeneralLinearRegressionResult : ResidualsResult {

        internal GeneralLinearRegressionResult () : base() {
            this.anova = new Lazy<OneWayAnovaResult>(CreateAnova);
        }

        private readonly Lazy<OneWayAnovaResult> anova;

        /// <summary>
        /// Gets an estimate of the intercept.
        /// </summary>
        public abstract UncertainValue Intercept { get; }

        /// <summary>
        /// Gets r<sup>2</sup> for the regression. 
        /// </summary>
        public virtual double RSquared {
            get {
                return anova.Value.Factor.SumOfSquares / anova.Value.Total.SumOfSquares;
            }
        }


        /// <summary>
        /// Gets an F-test for the regression.
        /// </summary>
        public virtual TestResult F {
            get {
                return anova.Value.Factor.Result;
            }
        }

        /// <summary>
        /// Gets an ANOVA analysis of the regression.
        /// </summary>
        public virtual OneWayAnovaResult Anova {
            get {
                return anova.Value;
            }
        }

        /// <inheritdoc />
        public override double SumOfSquaredResiduals {
            get {
                return anova.Value.Residual.SumOfSquares;
            }
        }

        internal abstract OneWayAnovaResult CreateAnova ();

    }

}
