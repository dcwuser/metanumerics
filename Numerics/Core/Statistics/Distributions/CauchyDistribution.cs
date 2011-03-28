using System;
using System.Collections.Generic;
using System.Text;

namespace Meta.Numerics.Statistics.Distributions {

#if FUTURE

    /// <summary>
    /// Represents a Cauchy distribution.
    /// </summary>
    public sealed class CauchyDistribution : Distribution {

        /// <summary>
        /// Initializes a new Cauchy distribution.
        /// </summary>
        /// <param name="mu"></param>
        /// <param name="gamma"></param>
        public CauchyDistribution (double mu, double gamma) {
            if (gamma <= 0.0) throw new ArgumentNullException("gamma");
            this.mu = mu;
            this.gamma = gamma;
        }

        private double mu;
        private double gamma;

        /// <summary>
        /// Gets the full width at half maximum of the Cauchy distribution.
        /// </summary>
        /// <remarks>
        /// <para>The full-width at half maximum is the width of the Cauchy peak at half its maximum value.</para>
        /// </remarks>
        public double FullWithAtHalfMaximum {
            get {
                return (gamma);
            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            return (1.0 / (1.0 + MoreMath.Pow2((x - mu) / gamma)) / Math.PI / gamma);
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (mu);
            }
        }

        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentNullException("n");
            } else if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (mu);
            } else {
                return (Double.NaN);
            }
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) {
                throw new ArgumentNullException("n");
            } else if (n == 0) {
                return (1.0);
            } else {
                if (n % 2 == 0) {
                    return (Double.PositiveInfinity);
                } else {
                    return (0.0);
                }
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            // expand for |x-mu| >> gamma
            return (0.5 + Math.Atan((x - mu) / gamma) / Math.PI);
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentNullException("P");
            return (mu + gamma * Math.Tan(Math.PI * (P - 0.5)));
        }

    }

#endif

}
