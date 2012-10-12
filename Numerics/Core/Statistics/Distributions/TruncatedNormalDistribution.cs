using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    // Start with a standard normal distribution with p(x) ~ e^{-x^2 / 2}, and truncate it at \pm D.
    // You can determine its required normalization from
    //   \int_{-D}^{+D} e^{-x^2 / 2} = \sqrt{2\pi} \erf(D/\sqrt{2})
    // Moments are also straightforward integrals

    public class TruncatedNormalDistribution : Distribution {

        public TruncatedNormalDistribution (double D) {
            if (D <= 0.0) throw new ArgumentOutOfRangeException("D");
            this.D = D;
            this.ED = AdvancedMath.Erf(D / Global.SqrtTwo);
        }

        private double D, ED;

        public override Interval Support {
            get {
                return (Interval.FromEndpoints(-D, D));
            }
        }

        public override double ProbabilityDensity (double x) {
            if (Math.Abs(x) < D) {
                return (Math.Exp(-x * x / 2.0) /  Global.SqrtTwoPI / ED);
            } else {
                return (0.0);
            }
        }

        public override double Median {
            get {
                return (0.0);
            }
        }

        public override double  Mean {
	        get {
                return (0.0);
	        }
        }

        public override double Variance {
            get {
                return (1.0 - Math.Sqrt(2.0 / Math.PI) * D * Math.Exp(-D * D / 2.0) / ED);
            }
        }

        public override double Skewness {
            get {
                return (0.0);
            }
        }

        public override double Moment (int n) {
            return (MomentAboutMean(n));
        }

        public override double MomentAboutMean (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else if (n % 2 != 0) {
                return (0.0);
            } else {
                // Using
                //   \int_0^D x^{2m} e^{-x^2 / 2} = 2^{m-1/2} \left[ \Gamma(m+1/2) - \Gamma(m+1/2, D^2 / 2) \right]
                // and \Gamma(m + 1/2) = \frac{(2m-1)!!}{2^m} \sqrt{\pi} and normalization above, you get
                //   C_{2m} = \frac{(2m-1)!!}{\erf(D/\sqrt{2}) \left[ 1 - \frac{\Gamma(m+1/2, D^2 / 2)}{\Gamma(m + 1/2)} \right]
                // which reduced to standard normal (2m-1)!! in D \rightarrow \infty limit as expected.
                return (AdvancedIntegerMath.DoubleFactorial(n - 1) * (1.0 - AdvancedMath.RightRegularizedGamma(n / 2 + 0.5, D * D / 2.0)) / ED);
            }
        }

        private double CentralProbability (double x) {
            return (AdvancedMath.Erf(x / Global.SqrtTwo) / ED);
        }

        private double TailProbability (double x) {
            return (AdvancedMath.Erfc(x / Global.SqrtTwo) / ED);
        }

        public override double LeftProbability (double x) {
            if (x < -D) {
                return (0.0);
            } else if (x < 0.0) {
                return (TailProbability(-x) / 2.0);
            } else if (x < D) {
                return ((1.0 - CentralProbability(x)) / 2.0);
            } else {
                return (1.0);
            }
        }

        public override double RightProbability (double x) {
            if (x < -D) {
                return (1.0);
            } else if (x < 0.0) {
                return ((1.0 - CentralProbability(-x)) / 2.0);
            } else if (x < D) {
                return (TailProbability(x) / 2.0);
            } else {
                return (0.0);
            }
        }
    }

}
