using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a &#x3C7; distribution.
    /// </summary>
    /// <remarks>
    /// <para>The &#x3C7; distribution is the distribution of the magnitude of a &#x3BD;-dimensional vector whose individual components are
    /// normally distribution.</para>
    /// <para>Since the distribution of the same sum of squares without taking the square root is called the &#x3BD;<sup>2</sup> distribution,
    /// (<see cref="ChiSquaredDistribution"/>), this distribution is called the &#x3BD; distribution.</para>
    /// </remarks>
    /// <seealso cref="ChiSquaredDistribution"/>
    /// <seealso href="https://en.wikipedia.org/wiki/Chi_distribution"/>
    public sealed class ChiDistribution : ContinuousDistribution {

        /// <summary>
        /// Initializes a new &#x3C7; distribution.
        /// </summary>
        /// <param name="nu">The number of degrees of freedom, which must be positive.</param>
        public ChiDistribution (int nu) {
            if (nu < 1) throw new ArgumentOutOfRangeException(nameof(nu));
            this.nu = nu;
            this.gamma = new GammaDistribution(0.5 * nu);
        }

        private readonly int nu;
        private readonly GammaDistribution gamma;

        /// <summary>
        /// Gets the number of degrees of freedom &#x3BD; of the distribution.
        /// </summary>
        public int DegreesOfFreedom {
            get {
                return (nu);
            }
        }

        /// <inheritdoc/>
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

        /// <inheritdoc/>
        public override double ProbabilityDensity (double x) {
            if (x < 0.0) {
                return (0.0);
            } else {
                return (x * gamma.ProbabilityDensity(0.5 * x * x));
            }
        }

        /// <inheritdoc/>
        public override double LeftProbability (double x) {
            if (x < 0.0) {
                return (0.0);
            } else {
                return (gamma.LeftProbability(0.5 * x * x));
            }
        }

        /// <inheritdoc/>
        public override double RightProbability (double x) {
            if (x <= 0.0) {
                return (1.0);
            } else {
                return (gamma.RightProbability(0.5 * x * x));
            }
        }

        /// <inheritdoc/>
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            return (Math.Sqrt(2.0 * gamma.InverseLeftProbability(P)));
        }

        /// <inheritdoc/>
        public override double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException(nameof(Q));
            return (Math.Sqrt(2.0 * gamma.InverseRightProbability(Q)));
        }

        /// <inheritdoc/>
        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else {
                return (Math.Pow(2.0, 0.5 * r) * AdvancedMath.Pochhammer(0.5 * nu, 0.5 * r));
            }
        }

        /// <inheritdoc/>
        public override double Variance {
            get {
                // C_2 = M_2 - M_1^2 = \nu - 2 [ (\nu/2)_{1/2} ]^2
                // Amazingly, these two terms differ by almost exactly 1/2 for all \nu.
                // For large \nu, this implies significant cancelation. We can use the asymptotic
                // expansion for the ratio of Gamma functions differing by 1/2 to derive
                //    C_2 = \frac{1}{2} - \frac{1}{8\nu} - \frac{1}{16\nu^2} + \frac{5}{128\nu^3} + \cdots
                // which we use for large \nu. 
                if (nu < 32) {
                    return (nu - 2.0 * MoreMath.Sqr(AdvancedMath.Pochhammer(0.5 * nu, 0.5)));
                } else {
                    double C = varianceCoefficients[0];
                    double nuk = 1.0;
                    for (int k = 1; k < varianceCoefficients.Length; k++) {
                        double C_old = C;
                        nuk /= nu;
                        C += varianceCoefficients[k] * nuk;
                        if (C == C_old) break;
                    }
                    return (C);
                }
            }
        }

        private static readonly double[] varianceCoefficients = new double[] {
            1.0 / 2.0, -1.0 / 8.0, - 1.0 / 16.0, 5.0 / 128.0, 23.0 / 256.0, -53.0 / 1024.0, -593.0 / 2048.0,
            5165.0 / 32768.0, 110123.0 / 65536.0, -231743.0 / 262144.0
        };

    }
}
