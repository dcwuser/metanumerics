using System;
using System.Diagnostics;
using System.Collections.Generic;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Describes the result of a regression.
    /// </summary>
    public abstract class RegressionResult : FitResult {

        internal RegressionResult () : base() {

        }

        private readonly IReadOnlyList<double> residuals;
        private readonly double sigma;

        /// <summary>
        /// Gets the residuals from the fit.
        /// </summary>
        public IReadOnlyList<double> Residuals {
            get {
                return (residuals);
            }
        }

        /// <summary>
        /// Gets standard error of the regression.
        /// </summary>
        public virtual double StandardError {
            get {
                return (sigma);
            }
        }

    }

}
