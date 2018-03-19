using System;
using System.Collections.Generic;

namespace Meta.Numerics.Statistics.Distributions
{
    /// <summary>
    /// Represents the result of a fit to the exponential distribution.
    /// </summary>
    public sealed class ExponentialFitResult : DistributionFitResult<ExponentialDistribution>
    {
        internal ExponentialFitResult(UncertainValue mu, TestResult goodnessOfFit) : base(new ExponentialDistribution(mu.Value), goodnessOfFit)
        {
            this.mu = mu;
        }

        private readonly UncertainValue mu;

        /// <summary>
        /// Gets an estimate, with uncertainty, of the mean of the exponential distribution.
        /// </summary>
        public UncertainValue Mean
        {
            get
            {
                return (mu);
            }
        }

        internal override ParameterCollection CreateParameters()
        {
            ParameterCollection parameters = new ParameterCollection(
                nameof(Mean), mu.Value, MoreMath.Sqr(mu.Uncertainty)
            );
            return (parameters);
        }

    }
}
