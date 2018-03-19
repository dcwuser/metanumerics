using System;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {
    /// <summary>
    /// Represents the result of a sample to a normal distribution.
    /// </summary>
    public sealed class NormalFitResult : DistributionFitResult<NormalDistribution> {

        internal NormalFitResult (UncertainValue mu, UncertainValue sigma, NormalDistribution distribution, TestResult goodnessOfFit) : base(distribution, goodnessOfFit) {
            this.mu = mu;
            this.sigma = sigma;
        }

        private readonly UncertainValue mu;
        private readonly UncertainValue sigma;

        internal override ParameterCollection CreateParameters () {
            return (new ParameterCollection(
                nameof(Mean), mu.Value, MoreMath.Sqr(mu.Uncertainty),
                nameof(StandardDeviation), sigma.Value, MoreMath.Sqr(sigma.Uncertainty),
                0.0
            ));
        }

        /// <summary>
        /// Gets an estimate, with uncertainty, of the mean.
        /// </summary>
        public UncertainValue Mean {
            get {
                return (mu);
            }
        }

        /// <summary>
        /// Gets and estimate, with uncertainty, of the standard deviation.
        /// </summary>
        public UncertainValue StandardDeviation {
            get {
                return (sigma);
            }
        }

    }
}
