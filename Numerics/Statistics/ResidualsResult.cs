using System;
using System.Collections.Generic;
using System.Text;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Describes the result of a fit with residuals.
    /// </summary>
    public abstract class ResidualsResult : FitResult {

        internal ResidualsResult () : base() { }

        /// <summary>
        /// Gets a list of the residuals of the fit.
        /// </summary>
        /// <value>A read-only list, in the same order as the original data, of the difference between each measured and predicted value.</value>
        public abstract IReadOnlyList<double> Residuals { get; }

        /// <summary>
        /// Gets the sum of the squares of all residuals.
        /// </summary>
        /// <value>The sum of the squares of all <see cref="Residuals"/>.</value>
        public abstract double SumOfSquaredResiduals { get; }

        // We could implement this, but everyone who inherits would override, so we don't.

    }
}
