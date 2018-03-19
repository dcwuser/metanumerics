using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Contains the result of a fit to a Rayleigh distribution.
    /// </summary>
    public sealed class RayleighFitResult : DistributionFitResult<RayleighDistribution> {

        internal RayleighFitResult (UncertainValue scale, RayleighDistribution distribution, TestResult goodnessOfFit) : base(distribution, goodnessOfFit)
        {
            this.scale = scale;
        }

        private readonly UncertainValue scale;

        internal override ParameterCollection CreateParameters () {
            return (new ParameterCollection(nameof(Scale), scale.Value, MoreMath.Sqr(scale.Uncertainty)));
        }

        /// <summary>
        /// Gets an estimate, with uncertainty, of the scale parameter.
        /// </summary>
        public UncertainValue Scale {
            get {
                return (scale);
            }
        }
    }
}
