using System;


using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents the distribution of Student't t statistic.
    /// </summary>
    /// <remarks><para>The mean of n independent standard-normal distributed variables, divided by their root mean square,
    /// is distributed according a Student distribution with n degrees of freedom. Since this is the form of the expression
    /// for the mean of a sample divided by its standard deviation, the Student distribution expresses the distribution of
    /// sample means arround the population mean, for a normally distributed population.</para></remarks>
    /// <seealso cref="Sample.StudentTTest(Double)"/>
    /// <seealso href="http://en.wikipedia.org/wiki/Student_t_distribution" />
    public sealed class StudentDistribution : Distribution {

        private double nu;

        /// <summary>
        /// Initializes a new Student distribution.
        /// </summary>
        /// <param name="nu">The number of degrees of freedom, which must be positive.</param>
        public StudentDistribution (double nu) {
            if (nu <= 0.0) throw new ArgumentOutOfRangeException("nu");
            this.nu = nu;
        }

        /// <summary>
        /// Gets the number of degrees of freedom.
        /// </summary>
        public double DegreesOfFreedom {
            get {
                return (nu);
            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            return (Math.Pow(1.0 + x * x / nu, -0.5 * (nu + 1.0)) / AdvancedMath.Beta(0.5, 0.5 * nu) / Math.Sqrt(nu));
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            double t2 = x * x;
            if (x < 0.0) {
                return (0.5 * AdvancedMath.Beta(0.5 * nu, 0.5, nu / (nu + t2)) / AdvancedMath.Beta(0.5, 0.5 * nu));
            } else {
                return (0.5 * (1.0 + AdvancedMath.Beta(0.5, 0.5 * nu, t2 / (nu + t2)) / AdvancedMath.Beta(0.5, 0.5 * nu)));
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            double t2 = x * x;
            if (x < 0.0) {
                return ((1.0 + AdvancedMath.LeftRegularizedBeta(1.0 / 2.0, nu / 2.0, t2 / (nu + t2))) / 2.0);
                //return (0.5 * (1.0 + AdvancedMath.Beta(0.5, 0.5 * nu, t2 / (nu + t2)) / AdvancedMath.Beta(0.5, 0.5 * nu)));
            } else {
                return (AdvancedMath.LeftRegularizedBeta(nu / 2.0, 1.0 / 2.0, nu / (nu + t2)) / 2.0);
                //return (0.5 * AdvancedMath.Beta(0.5 * nu, 0.5, nu / (nu + t2)) / AdvancedMath.Beta(0.5, 0.5 * nu));
            }
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");

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
            if (rng == null) throw new ArgumentNullException("rng");

            // this is a modified form of Box-Mueller due to Bailey

            double u, v, w;
            do {
                u = 2.0 * rng.NextDouble() - 1.0;
                v = 2.0 * rng.NextDouble() - 1.0;
                w = u * u + v * v;
            } while (w > 1.0);

            return (u * Math.Sqrt(nu * (Math.Pow(w, -2.0 / nu) - 1.0) / w));

            // the corresponding v-deviate is not independent of the u-deviate,
            // so don't attempt store it and use it for the next call

        }

        /// <inheritdoc />
        public override double Moment (int r) {
            if (r < 0) throw new ArgumentOutOfRangeException("r");

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
        public override double MomentAboutMean (int r) {
            return (Moment(r));
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                if (nu > 1.0) {
                    return (0.0);
                } else {
                    return (Double.NaN);
                }
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                if (nu > 2.0) {
                    return (nu / (nu - 2.0));
                } else {
                    return (Double.PositiveInfinity);
                }
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                if (nu > 3) {
                    return (0.0);
                } else {
                    return (Double.NaN);
                }
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (0.0);
            }
        }

    }

}