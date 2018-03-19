using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents the result of a sample to a normal distribution.
    /// </summary>
    public sealed class WeibullFitResult : DistributionFitResult<WeibullDistribution> {

        internal WeibullFitResult (double scale, double shape, SymmetricMatrix covariance, WeibullDistribution distribution, TestResult goodnessOfFit) :
            base(distribution, goodnessOfFit) {
            this.scale = scale;
            this.shape = shape;
            this.covariance = covariance;
        }

        private readonly double scale;
        private readonly double shape;
        private readonly SymmetricMatrix covariance;

        internal override ParameterCollection CreateParameters () {
            ParameterCollection parameters = new ParameterCollection(
                new string[] { nameof(Scale), nameof(Shape) },
                new ColumnVector(scale, shape), covariance
            );
            return (parameters);
        }

        /// <summary>
        /// Gets an estimate, with uncertainty, of the scale parameter.
        /// </summary>
        public UncertainValue Scale {
            get {
                return (new UncertainValue(scale, Math.Sqrt(covariance[0, 0])));
            }
        }

        /// <summary>
        /// Gets and estimate, with uncertainty, of the shape parameter.
        /// </summary>
        public UncertainValue Shape {
            get {
                return (new UncertainValue(shape, Math.Sqrt(covariance[1,1])));
            }
        }

    }

}
