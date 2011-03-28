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
                return (0.5 * (1.0 + AdvancedMath.Beta(0.5, 0.5 * nu, t2 / (nu + t2)) / AdvancedMath.Beta(0.5, 0.5 * nu)));
            } else {
                return (0.5 * AdvancedMath.Beta(0.5 * nu, 0.5, nu / (nu + t2)) / AdvancedMath.Beta(0.5, 0.5 * nu));
            }
        }

        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0) throw new ArgumentOutOfRangeException("n");

            if ((n % 2) == 0) {
                if (n < nu) {
                    double m = 1.0;
                    for (int j = 1; j <= n / 2; j++) {
                        m *= (2 * j - 1) * nu / (nu - 2 * j);
                    }
                    return (m);
                } else {
                    return (Double.PositiveInfinity);
                }
            } else {
                return (0.0);
            }

        }

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            return (Moment(n));
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (0.0);
            }
        }

        /// <inheritdoc />
        public override double StandardDeviation {
            get {
                if (nu <= 2.0) {
                    return (Double.PositiveInfinity);
                } else {
                    return (Math.Sqrt(nu / (nu - 2.0)));
                }
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (0.0);
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