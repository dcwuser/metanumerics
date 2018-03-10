using System;
using System.Diagnostics;
using System.Collections.Generic;

namespace Meta.Numerics.Statistics {

    // Rename to RegressionResult or FitResult

    /// <summary>
    /// The base class of all fit results.
    /// </summary>
    public class BaseFitResult {

        internal BaseFitResult (ParameterCollection parameters, double logLikelihood) {
            Debug.Assert(parameters != null);
            this.parameters = parameters;
            this.logLikelihood = logLikelihood;
        }

        private readonly ParameterCollection parameters;

        private readonly double logLikelihood;

        /// <summary>
        /// Gets the parameters of the regression.
        /// </summary>
        public virtual ParameterCollection Parameters {
            get {
                return (parameters);
            }
        }

        /// <summary>
        /// Gets the log-likelihood value of the fit.
        /// </summary>
        /// <see href="https://en.wikipedia.org/wiki/Likelihood_function"/>
        public virtual double LogLikelihood {
            get {
                return (logLikelihood);
            }
        }

    }

    // Rename to NonlinearRegressionResult or GeneralRegressionResult
    // Make residuals into a list

    /// <summary>
    /// Represents the result of a regression fit.
    /// </summary>
    public class RegressionResult : BaseFitResult {

        internal RegressionResult (ParameterCollection parameters, double logLikelyhood, List<double> residuals) : base(parameters, logLikelyhood) {
            Debug.Assert(residuals != null);
            this.residuals = residuals;
        }

        private readonly List<double> residuals;

        /// <summary>
        /// Gets the residuals.
        /// </summary>
        public virtual IReadOnlyList<double> Residuals {
            get {
                return (residuals);
            }
        }

    }

    /// <summary>
    /// Represents the result of a general linear regression.
    /// </summary>
    public class GeneralLinearRegressionResult : RegressionResult {

        internal GeneralLinearRegressionResult (ParameterCollection parameters, OneWayAnovaResult anova, List<double> residuals) : base(parameters, Double.NaN, residuals) {
            Debug.Assert(anova != null);
            this.anova = anova;
        }

        private readonly OneWayAnovaResult anova;

        /// <summary>
        /// Gets an estimate of the intercept.
        /// </summary>
        public virtual UncertainValue Intercept {
            get {
                return (Parameters["Intercept"].Estimate);
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

    }

}
