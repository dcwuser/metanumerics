using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Statistics.Distributions
{
    /// <summary>
    /// Represents the result of a fit to a distribution.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public abstract class DistributionFitResult<T> : BaseFitResult
    {
        internal DistributionFitResult(T distribution, TestResult goodnessOfFit) : base()
        {
            this.Distribution = distribution;
            this.GoodnessOfFit = goodnessOfFit;
        }

        /// <summary>
        /// Gets the best-fit distribution.
        /// </summary>
        public virtual T Distribution { get; private set; }

        /// <summary>
        /// Gets a test of the quality of the fit of the distribution to the data.
        /// </summary>
        public virtual TestResult GoodnessOfFit { get; private set; }

    }
}
