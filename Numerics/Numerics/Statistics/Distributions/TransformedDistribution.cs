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

        // <x> = \int dx \, x \frac{dP}{dx} = \int dx \, (a + b y) \frac{dP}{dx}
        //     = a \int dx \frac{dP}{dx} + b \int dx y \frac{dP}{dx} \frac{dy}{dy}
        //     = a + b \int dy y \frac{dP}{dy} = a + b <y> 

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

        public override double Moment (int r) {
            if (shift == 0.0) {
                return (baseDistribution.Moment(r) * MoreMath.Pow(scale, r));
            } else {
                // If shift is non-zero, we need to do a binomial expansion inovling multiple y-moments.
                // In that case, just do integral.
                return(base.Moment(r));
            }
        }

        public override double MomentAboutMean (int r) {
            return (baseDistribution.MomentAboutMean(r) * MoreMath.Pow(scale, r));
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
