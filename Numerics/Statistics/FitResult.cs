using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

using Meta.Numerics;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// The base class of results for all fits.
    /// </summary>
    public abstract class FitResult {

        internal FitResult () {
            this.parameters = new Lazy<ParameterCollection>(CreateParameters);
            //this.logLikelihood = Double.NaN;
        }

        private readonly Lazy<ParameterCollection> parameters;

        /// <summary>
        /// Gets the parameters of the regression.
        /// </summary>
        public virtual ParameterCollection Parameters {
            get {
                return (parameters.Value);
            }
        }

        internal abstract ParameterCollection CreateParameters ();

        /*
        /// <summary>
        /// Gets the log-likelihood of the fit.
        /// </summary>
        /// <see href="https://en.wikipedia.org/wiki/Likelihood_function"/>
        public virtual double LogLikelihood {
            get {
                return (logLikelihood);
            }
        }
        */
    }

}

