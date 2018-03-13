using System;

namespace Meta.Numerics.Statistics.Distributions {
    /// <summary>
    /// Represents the result of fitting sample data to a Gumbel distribution.
    /// </summary>
    public sealed class GumbelFitResult : DistributionFitResult<GumbelDistribution> {

        internal GumbelFitResult (double m, double s, double mVar, double sVar, double msCov, TestResult goodnessOfFit) : base(new GumbelDistribution(m, s), goodnessOfFit)
        {
            this.m = m;
            this.s = s;
            this.mVar = mVar;
            this.sVar = sVar;
            this.msCov = msCov;
        }

        private readonly double m, s, mVar, sVar, msCov;

        /// <summary>
        /// Gets the best fit value of the location parameter and its associated uncertainty.
        /// </summary>
        public UncertainValue Location {
            get {
                return (new UncertainValue(m, Math.Sqrt(mVar)));
            }
        }

        /// <summary>
        /// Gets the best fit value of the scale parameter and its associated uncertainty.
        /// </summary>
        public UncertainValue Scale {
            get {
                return (new UncertainValue(s, Math.Sqrt(sVar)));
            }
        }

        internal override ParameterCollection CreateParameters () {
            return (new ParameterCollection(nameof(Location), m, mVar, nameof(Scale), s, sVar, msCov));
        }

    }

}
