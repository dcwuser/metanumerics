using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Fréchet distribution.
    /// </summary>
    public sealed class FrechetDistribution : ContinuousDistribution {

        /// <summary>
        /// Initializes a new instance of a Fréchet distribution with the given shape parameter.
        /// </summary>
        /// <param name="shape"></param>
        public FrechetDistribution (double shape) : this(shape, 1.0, 0.0) {

        }

        internal FrechetDistribution (double shape, double scale, double location) {
            if (shape <= 0.0) throw new ArgumentOutOfRangeException(nameof(shape));
            if (scale <= 0.0) throw new ArgumentOutOfRangeException(nameof(scale));
            this.a = shape;
            this.m = location;
            this.s = scale;
        }

        private readonly double a, m, s;


        /// <inheritdoc/>
        public override Interval Support {
            get {
                return Interval.FromEndpoints(m, Double.PositiveInfinity);
            }
        }

        /// <inheritdoc/>
        public override double ProbabilityDensity (double x) {
            if (x < m) {
                return (0.0);
            } else {
                double z = (x - m) / s;
                double za = 1.0 / Math.Pow(z, a);
                return (a / (s * z) * za * Math.Exp(-za));
            }
        }

        /// <inheritdoc/>
        public override double LeftProbability (double x) {
            if (x <= m) {
                return (0.0);
            } else {
                double z = (x - m) / s;
                double za = 1.0 / Math.Pow(z, a);
                return (Math.Exp(-za));
            }
        }

        /// <inheritdoc/>
        public override double RightProbability (double x) {
            if (x <= m) {
                return (1.0);
            } else {
                double z = (x - m) / s;
                double za = 1.0 / Math.Pow(z, a);
                return (-MoreMath.ExpMinusOne(-za));
            }
        }

        /// <inheritdoc/>
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            return (m + s * Math.Pow(-Math.Log(P), -1.0 / a));
        }

        /// <inheritdoc/>
        public override double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException(nameof(Q));
            return (m + s * Math.Pow(-MoreMath.LogOnePlus(-Q), -1.0 / a));
        }

        /// <inheritdoc/>
        public override double Median {
            get {
                return (m + s * Math.Pow(Global.LogTwo, -1.0 / a));
            }
        }

        /// <inheritdoc/>
        public override double Mean {
            get {
                if (a <= 1.0) {
                    return (Double.PositiveInfinity);
                } else {
                    return (m + s * AdvancedMath.Gamma(1.0 - 1.0 / a));
                }
            }
        }

        /// <inheritdoc/>
        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r < a) {
                // This only works for m = 0
                return (MoreMath.Pow(s, r) * AdvancedMath.Gamma(1.0 - r / a));
            } else {
                return (Double.PositiveInfinity);
            }
        }

    }
}
