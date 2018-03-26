using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Statistics.Distributions
{
    /// <summary>
    /// Represents the result of a fit to a distribution.
    /// </summary>
    /// <typeparam name="T">The type of distribution fit.</typeparam>
    public abstract class DistributionFitResult<T> : FitResult
    {
        internal DistributionFitResult(T distribution, TestResult goodnessOfFit) : base()
        {
            Debug.Assert(distribution != null);
            Debug.Assert(goodnessOfFit != null);
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
        /// <remarks>
        /// <para>The null hypothesis of the test is that the sample is drawn from
        /// the given distribution.</para>
        /// </remarks>
        public virtual TestResult GoodnessOfFit { get; private set; }

    }
}
