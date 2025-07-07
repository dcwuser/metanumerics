using System;


using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents the distribution of Student's t statistic.
    /// </summary>
    /// <remarks><para>The mean of n independent standard-normal distributed variables, divided by their root mean square,
    /// is distributed according a Student distribution with n degrees of freedom. Since this is the form of the expression
    /// for the mean of a sample divided by its standard deviation, the Student distribution expresses the distribution of
    /// sample means arround the population mean, for a normally distributed population.</para>
    /// <para>The origin of the name Student distribution is a nice bit of statistical trivia. William Gosset published a
    /// paper describing the distribution and its statistical applications in 1908. His employer, the Guiness Brewery,
    /// did not want other brewers to be tipped off to its application to comparing small samples of beer, so they asked
    /// him to publish under a pseudonym. He published the paper under the name "Student" and the distribution has been
    /// known by that name since.
    /// </para>
    /// </remarks>
    /// <seealso cref="Univariate.StudentTTest(System.Collections.Generic.IReadOnlyCollection{double}, double)" />
    /// <seealso cref="Univariate.StudentTTest(System.Collections.Generic.IReadOnlyCollection{double}, System.Collections.Generic.IReadOnlyCollection{double})"/>
    /// <seealso href="https://en.wikipedia.org/wiki/Student%27s_t-distribution" />
    /// <seealso href="https://mathworld.wolfram.com/Studentst-Distribution.html"/>
    public sealed class StudentDistribution : ContinuousDistribution {

        private readonly double nu;

        /// <summary>
        /// Initializes a new Student distribution.
        /// </summary>
        /// <param name="nu">The number of degrees of freedom, which must be positive.</param>
        public StudentDistribution (double nu) {
            if (nu <= 0.0) throw new ArgumentOutOfRangeException(nameof(nu));
            this.nu = nu;
        }

        /// <summary>
        /// Gets the number of degrees of freedom.
        /// </summary>
        public double DegreesOfFreedom {
            get {
                return nu;
            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            return Math.Pow(1.0 + x * x / nu, -0.5 * (nu + 1.0)) / AdvancedMath.Beta(0.5, 0.5 * nu) / Math.Sqrt(nu);
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            double xSquared = x * x;
            if (x < 0.0) {
                return 0.5 * AdvancedMath.LeftRegularizedBeta(0.5 * nu, 0.5, nu / (nu + xSquared));
            } else {
                return 0.5 * (1.0 + AdvancedMath.Beta(0.5, 0.5 * nu, xSquared / (nu + xSquared)) / AdvancedMath.Beta(0.5, 0.5 * nu));
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            double xSquared = x * x;
            if (x < 0.0) {
                return 0.5 * (1.0 + AdvancedMath.LeftRegularizedBeta(0.5, 0.5 * nu, xSquared / (nu + xSquared)));
            } else {
                return 0.5 * AdvancedMath.LeftRegularizedBeta(0.5 * nu, 0.5, nu / (nu + xSquared));
            }
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));

            BetaDistribution b = new BetaDistribution(nu / 2.0, 1.0 / 2.0);
            if (P <= 0.5) {
                double x = b.InverseLeftProbability(2.0 * P);
                return (-Math.Sqrt(nu * (1.0 - x) / x));
            } else {
                double Q = 1.0 - P;
                double x = b.InverseLeftProbability(2.0 * Q);
                return (Math.Sqrt(nu * (1.0 - x) / x));
            }
        }

        /// <inheritdoc />
        public override double GetRandomValue (Random rng) {
            if (rng is null) throw new ArgumentNullException(nameof(rng));

            // this is a modified form of Box-Mueller due to Bailey

            double u, v, w;
            do {
                u = 2.0 * rng.NextDouble() - 1.0;
                v = 2.0 * rng.NextDouble() - 1.0;
                w = u * u + v * v;
            } while (w > 1.0);

            return u * Math.Sqrt(nu * (Math.Pow(w, -2.0 / nu) - 1.0) / w);

            // the corresponding v-deviate is not independent of the u-deviate,
            // so don't attempt store it and use it for the next call

        }

        /// <inheritdoc />
        public override double RawMoment (int r) {
            if (r < 0) throw new ArgumentOutOfRangeException(nameof(r));

            if ((r % 2) == 0) {
                if (r < nu) {
                    double m = 1.0;
                    for (int j = 1; j <= r / 2; j++) {
                        m *= (2 * j - 1) * nu / (nu - 2 * j);
                    }
                    return (m);
                } else {
                    return (Double.PositiveInfinity);
                }
            } else {
                if (r < nu) {
                    return (0.0);
                } else {
                    return (Double.NaN);
                }
            }

        }

        /// <inheritdoc />
        public override double CentralMoment (int r) {
            return RawMoment(r);
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                if (nu > 1.0) {
                    return 0.0;
                } else {
                    return Double.NaN;
                }
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                if (nu > 2.0) {
                    return nu / (nu - 2.0);
                } else {
                    return Double.PositiveInfinity;
                }
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                if (nu > 3.0) {
                    return 0.0;
                } else {
                    return Double.NaN;
                }
            }
        }

        /// <inheritdoc />
        public override double ExcessKurtosis {
            get {
                if (nu > 4.0) {
                    return 6.0 / (nu - 4.0);
                } else {
                    return Double.PositiveInfinity;
                }
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return 0.0;
            }
        }

    }

}