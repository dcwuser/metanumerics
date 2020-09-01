using System;

using Meta.Numerics.Analysis;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a standard skew nomral distribution.
    /// </summary>
    /// <seealso href="https://en.wikipedia.org/wiki/Skew_normal_distribution"/>
    public sealed class SkewNormalDistribution : ContinuousDistribution {

        /// <summary>
        /// Initializes a new instance of the skew nomral distribution with the given shape parameter.
        /// </summary>
        /// <param name="alpha">The shape parameter.</param>
        public SkewNormalDistribution(double alpha) {
            this.alpha = alpha;
            this.gamma = 1.0 / MoreMath.Hypot(1.0, alpha);
            this.delta = alpha * gamma;
            this.mu = delta * Math.Sqrt(2.0 / Math.PI);
        }

        private readonly double alpha;

        private readonly double gamma, delta, mu;

        private readonly NormalDistribution n = new NormalDistribution();

        /// <inheritdoc />
        public override double ProbabilityDensity(double x) {
            return (2.0 * n.ProbabilityDensity(x) * n.LeftProbability(alpha * x));
        }

        /// <inheritdoc />
        public override double LeftProbability(double x) {
            return n.LeftProbability(x) - 2.0 * OwenT(x, alpha);
        }

        /// <inheritdoc />
        public override double RightProbability(double x) {
            return n.RightProbability(x) + 2.0 * OwenT(x, alpha);
        }

        private static double OwenT(double h, double a) {
            if (h < 0.0) return OwenT(-h, a);
            if (a < 0.0) return -OwenT(h, -a);
            double g = 0.5 * h * h;
            double I = FunctionMath.Integrate(
                x => {
                    double w = 1.0 + x * x;
                    return Math.Exp(-g * w) / w;
                }, 0.0, a
            );
            return I / (2.0 * Math.PI);
        }

        /// <inheritdoc />
        public override double GetRandomValue(Random rng) {
            double u = n.GetRandomValue(rng);
            double v = n.GetRandomValue(rng);
            return (alpha * Math.Abs(u) + v) * gamma;
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return mu;
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return 1.0 - mu * mu;
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                double muSquared = mu * mu;
                return (4.0 - Math.PI) / 2.0 * muSquared * mu / Math.Pow(1.0 - muSquared, 3.0 / 2.0);
            }
        }

        /// <inheritdoc />
        public override double RawMoment(int r) {
            if (r % 2 == 0) {
                // Even raw moments are the same as for standard normal.
                return n.RawMoment(r);
            } else {
                return base.RawMoment(r);
            }
        }

    }

}
