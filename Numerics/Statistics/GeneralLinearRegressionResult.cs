using System;
using System.Diagnostics;
using System.Collections.Generic;

namespace Meta.Numerics.Statistics {

    // Rename to RegressionResult or FitResult

    /// <summary>
    /// The base class of all fit results.
    /// </summary>
    public class BaseFitResult {

        internal BaseFitResult ()
        {
            this.parameters = new Lazy<ParameterCollection>(CreateParameters);
            this.logLikelihood = Double.NaN;
        }

        internal BaseFitResult (ParameterCollection parameters, double logLikelihood) {
            Debug.Assert(parameters != null);
            this.parameters = new Lazy<ParameterCollection>(() => parameters);
            this.logLikelihood = logLikelihood;
        }

        //private readonly ParameterCollection parameters;

        private readonly Lazy<ParameterCollection> parameters;

        private readonly double logLikelihood;

        /// <summary>
        /// Gets the parameters of the regression.
        /// </summary>
        public virtual ParameterCollection Parameters {
            get {
                return (parameters.Value);
            }
        }

        internal virtual ParameterCollection CreateParameters ()
        {
            throw new InvalidOperationException();
        }

        /*
        /// <summary>
        /// Gets the log-likelihood value of the fit.
        /// </summary>
        /// <see href="https://en.wikipedia.org/wiki/Likelihood_function"/>
        public virtual double LogLikelihood {
            get {
                return (logLikelihood);
            }
        }
        */
    }

    /// <summary>
    /// Represents the result of a general linear regression.
    /// </summary>
    public abstract class GeneralLinearRegressionResult : BaseFitResult {

        internal GeneralLinearRegressionResult () : base() {
            this.anova = new Lazy<OneWayAnovaResult>(CreateAnova);
        }

        private readonly Lazy<OneWayAnovaResult> anova;

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
                return (anova.Value);
            }
        }

        internal abstract OneWayAnovaResult CreateAnova ();

    }

}
