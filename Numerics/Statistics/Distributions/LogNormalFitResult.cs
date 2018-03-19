using System;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Contains the result of a fit to a log-normal distribution.
    /// </summary>
    public sealed class LognormalFitResult : DistributionFitResult<LognormalDistribution> {

        internal LognormalFitResult (UncertainValue mu, UncertainValue sigma, LognormalDistribution distribution, TestResult goodnessOfFit) : base(distribution, goodnessOfFit) {
            this.mu = mu;
            this.sigma = sigma;
        }

        private readonly UncertainValue mu;
        private readonly UncertainValue sigma;

        internal override ParameterCollection CreateParameters () {
            return (new ParameterCollection(
                nameof(Mu), mu.Value, MoreMath.Sqr(mu.Uncertainty),
                nameof(Sigma), sigma.Value, MoreMath.Sqr(sigma.Uncertainty),
                0.0
            ));
        }

        /// <summary>
        /// Gets an estimate, with uncertainty, of the mean.
        /// </summary>
        public UncertainValue Mu {
            get {
                return (mu);
            }
        }

        /// <summary>
        /// Gets and estimate, with uncertainty, of the standard deviation.
        /// </summary>
        public UncertainValue Sigma {
            get {
                return (sigma);
            }
        }


    }
}
