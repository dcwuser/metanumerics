using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Contains the result of the fit of a sample to a Wald (Inverse Gaussian) distribution.
    /// </summary>
    public sealed class WaldFitResult : DistributionFitResult<WaldDistribution> {

        internal WaldFitResult(double mu, double lambda, double muVar, double lambdaVar, WaldDistribution distribution, TestResult goodnessOfFit) :
            base(distribution, goodnessOfFit) {
            this.mu = mu;
            this.lambda = lambda;
            this.muVar = muVar;
            this.lambdaVar = lambdaVar;
        }

        private readonly double mu, muVar;
        private readonly double lambda, lambdaVar;

        /// <summary>
        /// Gets an estimate, with uncertainty, of the mean of the underlying Wald distribution.
        /// </summary>
        public UncertainValue Mean {
            get {
                return (new UncertainValue(mu, Math.Sqrt(muVar)));
            }
        }

        /// <summary>
        /// Gets an estimate, with uncertainty, of the shape parameter of the underlying Wald distribution.
        /// </summary>
        public UncertainValue Shape {
            get {
                return (new UncertainValue(lambda, Math.Sqrt(lambdaVar)));
            }
        }

        internal override ParameterCollection CreateParameters () {
            return (new ParameterCollection(nameof(Mean), mu, muVar, nameof(Shape), lambda, lambdaVar, 0.0));
        }
    }
}
