using System;

namespace Meta.Numerics.Statistics.Distributions {

    internal class TransformedDistribution : Distribution {

        public TransformedDistribution (Distribution baseDistribution, double shift, double scale) {
            if (baseDistribution == null) throw new ArgumentNullException("baseDistribution");
            this.baseDistribution = baseDistribution;
            this.shift = shift;
            this.scale = scale;
        }

        private readonly Distribution baseDistribution;
        private readonly double shift;
        private readonly double scale;

        private double TransformXtoY (double x) {
            return ((x - shift) / scale);
        }

        private double TransformYtoX (double y) {
            return (shift + y * scale);
        }

        public override Interval Support {
	        get {
		        Interval ySupport = baseDistribution.Support;
                return(Interval.FromEndpoints(TransformYtoX(ySupport.LeftEndpoint), TransformYtoX(ySupport.RightEndpoint)));
	        }
        }

        public override double ProbabilityDensity (double x) {
            return (baseDistribution.ProbabilityDensity(TransformXtoY(x)) / scale);
        }

        public override double LeftProbability (double x) {
            return (baseDistribution.LeftProbability(TransformXtoY(x)));
        }

        public override double RightProbability (double x) {
            return (baseDistribution.RightProbability(TransformXtoY(x)));
        }

        public override double Median {
            get {
                return (TransformYtoX(baseDistribution.Median));
            }
        }

        public override double Mean {
            get {
                return (TransformYtoX(baseDistribution.Mean));
            }
        }
        public override double StandardDeviation {
            get {
                return (baseDistribution.StandardDeviation * scale);
            }
        }

        public override double Variance {
            get {
                return (baseDistribution.Variance * scale * scale);
            }
        }

        public override double Skewness {
            get {
                return (baseDistribution.Skewness);
            }
        }

        public override double MomentAboutMean (int n) {
            return (baseDistribution.MomentAboutMean(n) * MoreMath.Pow(scale, n));
        }

        public override double InverseLeftProbability (double P) {
            return (TransformYtoX(baseDistribution.InverseLeftProbability(P)));
        }

        public override double InverseRightProbability (double P) {
            return (TransformYtoX(baseDistribution.InverseRightProbability(P)));
        }

        public override double GetRandomValue (Random rng) {
            return (TransformYtoX(baseDistribution.GetRandomValue(rng)));
        }

    }

}
