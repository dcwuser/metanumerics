using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Statistics.Distributions
{
    /// <summary>
    /// Represents the result of a fit to a Bernoulli distribution.
    /// </summary>
    public sealed class BernoulliFitResult : DistributionFitResult<BernoulliDistribution>
    {
        internal BernoulliFitResult(UncertainValue p, TestResult goodnessOfFit) : base(new BernoulliDistribution(p.Value), goodnessOfFit)
        {
            this.p = p;
        }

        private readonly UncertainValue p;

        /// <summary>
        /// Gets an estimate, with uncertainty, of the probability of success parameter.
        /// </summary>
        public UncertainValue P
        {
            get
            {
                return (p);
            }
        }

        internal override ParameterCollection CreateParameters()
        {
            return (new ParameterCollection(nameof(P), p.Value, MoreMath.Sqr(p.Uncertainty)));
        }
    }
}
