using System;
using System.Diagnostics;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Contains the result of a fit of a sample to a Beta distribution.
    /// </summary>
    public sealed class BetaFitResult : DistributionFitResult<BetaDistribution> {

        internal BetaFitResult(ColumnVector ab, SymmetricMatrix C, BetaDistribution distribution, TestResult goodnessOfFit) : base(distribution, goodnessOfFit) {
            Debug.Assert(ab != null);
            Debug.Assert(C != null);
            Debug.Assert(ab.Dimension == 2);
            Debug.Assert(C.Dimension == 2);
            this.ab = ab;
            this.C = C;
        }

        private readonly ColumnVector ab;
        private readonly SymmetricMatrix C;

        /// <summary>
        /// Gets the alpha parameter.
        /// </summary>
        public UncertainValue Alpha {
            get {
                return (new UncertainValue(ab[0], Math.Sqrt(C[0, 0])));
            }
        }

        /// <summary>
        /// Gets the beta parameter.
        /// </summary>
        public UncertainValue Beta {
            get {
                return (new UncertainValue(ab[1], Math.Sqrt(C[1, 1])));
            }
        }

        internal override ParameterCollection CreateParameters () {
            return (new ParameterCollection(new string[] { nameof(Alpha), nameof(Beta) }, ab, C));
        }

    }
}
