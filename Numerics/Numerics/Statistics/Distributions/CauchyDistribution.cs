using System;
using System.Collections.Generic;
using System.Text;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Cauchy distribution.
    /// </summary>
    /// <remarks>
    /// <para>In physical applications, the Cauchy distribution is usually called a Lorentz distribution. It models
    /// the shape of a spectral line.</para>
    /// <para>The Cauchy distribution has "fat tails". In fact, it falls off at the minimum possible rate consistent
    /// with having a convergent integral. For this same reason, none of its moments (above the zeroth) are defined.</para>
    /// </remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/Cauchy_distribution"/>
    public sealed class CauchyDistribution : ContinuousDistribution {

        /// <summary>
        /// Initializes a new standard Cauchy distribution.
        /// </summary>
        public CauchyDistribution () : this(0.0, 1.0) {
        }

        /// <summary>
        /// Initializes a new Cauchy distribution.
        /// </summary>
        /// <param name="mu">The centroid of the distribution.</param>
        /// <param name="gamma">The width parameter of the distribution.</param>
        public CauchyDistribution (double mu, double gamma) {
            if (gamma <= 0.0) throw new ArgumentNullException("gamma");
            this.mu = mu;
            this.gamma = gamma;
            this.cauchyRng = DeviateGeneratorFactory.GetCauchyGenerator();
        }

        private readonly double mu;
        private readonly double gamma;

        private readonly IDeviateGenerator cauchyRng;

        /// <summary>
        /// Gets the full width at half maximum (FWHM) of the Cauchy distribution.
        /// </summary>
        /// <remarks>
        /// <para>The full-width at half maximum (FWHM) is the width of the distribution peak at half its maximum value.</para>
        /// </remarks>
        public double FullWithAtHalfMaximum {
            get {
                return (2.0 * gamma);
            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            double z = (x - mu) / gamma;
            return (1.0 / (1.0 + z * z) / (Math.PI * gamma));
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (Double.NaN);
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (mu);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (Double.NaN);
            }
        }

        /// <inheritdoc />
        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentNullException("r");
            } else if (r == 0) {
                return (1.0);
            } else {
                return (Double.NaN);
            }
        }

        /// <inheritdoc />
        public override double CentralMoment (int r) {
            if (r < 0) {
                throw new ArgumentNullException("r");
            } else if (r == 0) {
                return (1.0);
            } else {
                return (Double.NaN);
            }
        }

        /// <inheritdoc />
        public override double Cumulant (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException("r");
            } else if (r == 0) {
                return (0.0);
            } else {
                return (Double.NaN);
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            double z = (x - mu) / gamma;
            if (z < 0.0) {
                // To avoid cancelation, use pi/2 - atan(x)  = atan(1/x).
                return (Math.Atan(-1.0 / z) / Math.PI);
            } else {
                return (0.5 + Math.Atan(z) / Math.PI);
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            double z = (x - mu) / gamma;
            if (z <= 0.0) {
                return (0.5 - Math.Atan(z) / Math.PI);
            } else {
                // To avoid cancelation, use pi/2 - atan(x)  = atan(1/x).
                return (Math.Atan(1.0 / z) / Math.PI);
            }
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            return (InverseProbability(P, 1.0 - P));
        }

        /// <inheritdoc />
        public override double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException("Q");
            return (InverseProbability(1.0 - Q, Q));
        }

        // Code an inverse CDF that accepts both P and Q so that we can handle the case of either very small.

        private double InverseProbability (double P, double Q) {
            double z;
            if (P < Q) {
                z = -1.0 / MoreMath.TanPi(P);
            } else {
                z = 1.0 / MoreMath.TanPi(Q);
            }
            return (mu + z * gamma);
        }

        /// <inheritdoc />
        public override double GetRandomValue (Random rng) {
            if (rng == null) throw new ArgumentNullException("rng");
            double z = cauchyRng.GetNext(rng);
            return (mu + z * gamma);
        }

    }

}
