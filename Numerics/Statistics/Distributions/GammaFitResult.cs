using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Contains the result of a fit of a sample to a gamma distribution.
    /// </summary>
    public sealed class GammaFitResult : DistributionFitResult<GammaDistribution> {

        internal GammaFitResult (double shape, double scale, SymmetricMatrix covariance, GammaDistribution distribution, TestResult goodnessOfFit) : base(distribution, goodnessOfFit) {
            Debug.Assert(covariance != null);
            Debug.Assert(covariance.Dimension == 2);
            this.shape = shape;
            this.scale = scale;
            this.covariance = covariance;
        }

        private readonly double shape;
        private readonly double scale;
        private readonly SymmetricMatrix covariance;

        internal override ParameterCollection CreateParameters () {
            ColumnVector p = new ColumnVector(shape, scale);
            ParameterCollection parameters = new ParameterCollection(
                new string[] { nameof(Shape), nameof(Scale) }, p, covariance
            );
            return (parameters);
        }

        /// <summary>
        /// Gets an estimate, with uncertainty, of the mean.
        /// </summary>
        public UncertainValue Shape {
            get {
                return (new UncertainValue(shape, Math.Sqrt(covariance[0, 0])));
            }
        }

        /// <summary>
        /// Gets and estimate, with uncertainty, of the standard deviation.
        /// </summary>
        public UncertainValue Scale {
            get {
                return (new UncertainValue(scale, Math.Sqrt(covariance[1, 1])));
            }
        }

    }
}
