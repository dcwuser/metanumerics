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
    /// </remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/Cauchy_distribution"/>
    public sealed class CauchyDistribution : Distribution {

        /// <summary>
        /// Initializes a new standard Cauchy distribution.
        /// </summary>
        public CauchyDistribution () {
            this.mu = 0.0;
            this.gamma = 1.0;
        }

        /// <summary>
        /// Initializes a new Cauchy distribution.
        /// </summary>
        /// <param name="mu">The centroid of the distribution.</param>
        /// <param name="gamma">The full width at half maximum FWHM of the distribution.</param>
        public CauchyDistribution (double mu, double gamma) {
            if (gamma <= 0.0) throw new ArgumentNullException("gamma");
            this.mu = mu;
            this.gamma = gamma;
        }

        private double mu;
        private double gamma;

        /// <summary>
        /// Gets the full width at half maximum (FWHM) of the Cauchy distribution.
        /// </summary>
        /// <remarks>
        /// <para>The full-width at half maximum (FWHM) is the width of the distribution peak at half its maximum value.</para>
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
                return (Double.PositiveInfinity);
            }
        }

        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentNullException("n");
            } else if (n == 0) {
                return (1.0);
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
                    return (Double.NaN);
                }
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            double z = (x - mu) / gamma;
            if (z < -8.0) {
                // to get small values accurately without canelation, use the series for arctan(x) for large |x|
                return (ArcTanResidual(-z) / Math.PI);
            } else {
                return (0.5 + Math.Atan(z) / Math.PI);
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            double z = (x - mu) / gamma;
            if (z > 8.0) {
                // to get small values accurately without cancelation, use the series for arctan(x) for large |x|
                return (ArcTanResidual(z) / Math.PI);
            } else {
                return (0.5 - Math.Atan(z) / Math.PI);
            }
        }

        // For large argument, we have the series
        //   arctan(x) = \pm \pi / 2 - \sum_{k=0}^{\infty} (-1)^k / (2k+1) / z^{2k+1}
        // We implement ArcTanResidual to evaluate the sum, which gives the difference between arctan(x) and pi/2 for large x

        private static double ArcTanResidual (double z) {

            double f = 1.0;
            double z2 = -z * z;
            double zz = z2;
            for (int i = 1; i < Global.SeriesMax; i++) {
                double f_old = f;
                f += 1.0 / (2 * i + 1) / zz;
                if (f == f_old) return (f / z);
                zz *= z2;
            }
            throw new NonconvergenceException();

        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");

            double z;
            if (P < 0.0625) {
                // to preserve accuracy, use a series in the left tail
                z = TanNearHalfPi(Math.PI * P);
            } else if (P > 0.9375) {
                // and the right tail
                double Q = 1.0 - P;
                z = -TanNearHalfPi(Math.PI * Q);
            } else {
                // this is the simple inversion of the simple CDF formula
                z = Math.Tan(Math.PI * (P - 0.5));
            }
            return (mu + z * gamma);
        }

        /// <inheritdoc />
        public override double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException("Q");

            double z;
            if (Q < 0.0625) {
                // to preserve accuracy, use a series in the right tail
                z = -TanNearHalfPi(Math.PI * Q);
            } else if (Q > 0.9375) {
                // and the right tail
                double P = 1.0 - Q;
                z = TanNearHalfPi(Math.PI * P);
            } else {
                // this is the simple inversion of the simple CDF formula
                z = -Math.Tan(Math.PI * (Q - 0.5));
            }
            return (mu + z * gamma);
        }

        // tan(\pi/2 + e) = -(1/x) [ 1 + \sum_{k=1}^{\infty} (-1)^k B_{2k} / (2k)! (2 e)^{2k} ]
        // isn't the also cotangent expansion?

        private static double TanNearHalfPi (double x) {

            double f = 1.0;
            double x2 = MoreMath.Pow2(2.0 * x);
            double xx = 1.0;
            for (int i = 1; i < AdvancedIntegerMath.Bernoulli.Length; i++) {
                double f_old = f;
                xx *= -x2 / (2 * i) / (2 * i - 1);
                f += AdvancedIntegerMath.Bernoulli[i] * xx;
                if (f == f_old) return (-f / x);
            }
            throw new NonconvergenceException();

        }

    }

}
