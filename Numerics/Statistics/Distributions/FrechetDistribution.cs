﻿using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Fréchet distribution.
    /// </summary>
    /// <remarks>
    /// <para>The Fréchet distribution is also called the inverse Weibull distribution.</para>
    /// </remarks>
    /// <seealso href="https://en.wikipedia.org/wiki/Fr%C3%A9chet_distribution"/>
    public sealed class FrechetDistribution : ContinuousDistribution {

        /// <summary>
        /// Initializes a new instance of a Fréchet distribution with the given shape parameter.
        /// </summary>
        /// <param name="shape">The shape parameter, which must be positive.</param>
        public FrechetDistribution (double shape) : this(shape, 1.0) {

        }

        /// <summary>
        /// Initializes a new instance of a Fréchet distribution with the given shape and scale parameters.
        /// </summary>
        /// <param name="shape">The shape parameter, which must be positive.</param>
        /// <param name="scale">The scale parameter, which must be positive.</param>
        public FrechetDistribution (double shape, double scale) {
            if (shape <= 0.0) throw new ArgumentOutOfRangeException(nameof(shape));
            if (scale <= 0.0) throw new ArgumentOutOfRangeException(nameof(scale));
            this.a = shape;
            this.s = scale;
            this.m = 0.0;
        }

        private readonly double a, s, m;

        /// <summary>
        /// Gets the shape parameter of the distribution.
        /// </summary>
        public double Shape {
            get {
                return(a);
            }
        }

        /// <summary>
        /// Gets the scale parameter of the distribution.
        /// </summary>
        public double Scale {
            get {
                return (s);
            }
        }

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
        public override double Variance {
            get {
                if (a <= 2.0) {
                    return (Double.PositiveInfinity);
                } else if (a <= 4.0) {
                    // If a gets even moderately big, this suffers from cancelation
                    // For large a, value goes to \frac{\pi^2}{6 a^2}, but higher
                    // terms fall off only slowly and are not trivial to compute
                    double ai = 1.0 / a;
                    double g1 = AdvancedMath.Gamma(1.0 - 2.0 * ai);
                    double g2 = MoreMath.Sqr(AdvancedMath.Gamma(1.0 - ai));
                    return (s * s * (g1 - g2));
                } else {
                    return (base.Variance);
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

        /// <inheritdoc/>
        public override double CentralMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (0.0);
            } else if (r < a) {
                return (base.CentralMoment(r));
            } else {
                return (Double.PositiveInfinity);
            }
        }

    }
}
