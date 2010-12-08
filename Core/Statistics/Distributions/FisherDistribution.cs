using System;
using System.Collections.Generic;

using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents the distribution of Fisher's F-statistic.
    /// </summary>
    /// <remarks>
    /// <para>The ratio of the variances of two sets of normally distributed variables is distributed according to Fisher's F-distribution.</para>
    /// <img src="../images/FisherFromNormal.png" />
    /// <para>Many test statistics are ratios of variances and are therefore distributed according to the F-distribution. These include
    /// the F-test (<see cref="Sample.FisherFTest"/>),
    /// the goodness-of-fit test for a multi-linear regression (<see cref="MultivariateSample.LinearRegression(int)"/>),
    /// and ANOVA tests (<see cref="Sample.OneWayAnovaTest(IList{Sample})"/>).</para>
    /// </remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/F_distribution"/>
    public class FisherDistribution : Distribution {

        private double nu1;
        private double nu2;

        /// <summary>
        /// Gets the number of degrees of freedom in the numerator.
        /// </summary>
        public double NumeratorDegreesOfFreedom {
            get {
                return (nu1);
            }
        }

        /// <summary>
        /// Gets the number of degrees of freedom in the denominator.
        /// </summary>
        public double DenominatorDegreesOfFreedom {
            get {
                return (nu2);
            }
        }

    
        /// <summary>
        /// Instantiates a new Fisher distribution.
        /// </summary>
        /// <param name="nu1">The number of degrees of freedom in the numerator, which must be positive.</param>
        /// <param name="nu2">The number of degrees of freedom in the denominator, which must be positive.</param>
        /// <seealso cref="Sample.FisherFTest"/>
		public FisherDistribution (double nu1, double nu2) {
			if (nu1 <= 0.0) throw new ArgumentOutOfRangeException("nu1");
			if (nu2 <= 0.0) throw new ArgumentOutOfRangeException("nu2");
			this.nu1 = nu1;
			this.nu2 = nu2;
		}

        /// <inheritdoc />
        public override double ProbabilityDensity (double F) {
            if (F <= 0.0) {
                return (0.0);
            } else {
                double N = Math.Pow(nu1, 0.5 * nu1) * Math.Pow(nu2, 0.5 * nu2) /
                    AdvancedMath.Beta(0.5 * nu1, 0.5 * nu2);
                return (N * Math.Pow(F, 0.5 * nu1 - 1.0) * Math.Pow(nu2 + nu1 * F, -0.5 * (nu1 + nu2)));
            }
		}

        /// <inheritdoc />
        public override double LeftProbability (double F) {
            if (F <= 0.0) {
                return (0.0);
            } else {
                return (AdvancedMath.Beta(0.5 * nu1, 0.5 * nu2, nu1 * F / (nu2 + nu1 * F)) / AdvancedMath.Beta(0.5 * nu1, 0.5 * nu2));
            }
		}

        /// <inheritdoc />
        public override double RightProbability (double F) {
            if (F <= 0.0) {
                return (1.0);
            } else {
                return (AdvancedMath.Beta(0.5 * nu2, 0.5 * nu1, nu2 / (nu2 + nu1 * F)) / AdvancedMath.Beta(0.5 * nu2, 0.5 * nu1));
            }
		}

        /// <inheritdoc />
        public override double Mean {
			get {
                if (nu2 > 2.0) {
                    return (nu2 / (nu2 - 2.0)); 
                } else {
                    return (System.Double.PositiveInfinity);
                }
			}
		}

        /// <inheritdoc />
        public override double Variance {
            get {
                if (nu2 > 4.0) {
                    double m = Mean;
                    return (2.0 * m * m * (nu1 + nu2 - 2.0) / nu1 / (nu2 - 4.0));
                } else {
                    return (System.Double.PositiveInfinity);
                }
            }
        }

        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (Mean);
            } else {
                if (nu2 <= 2.0 * n) {
                    return (System.Double.PositiveInfinity);
                } else {

                    // (nu2)^n  Gamma(nu1/2 + n)  Gamma(nu2/2 - n)
                    // (---)    ----------------  ----------------
                    // (nu1)      Gamma(nu1/2)      Gamma(nu2/2)

                    // this can be computed using the recursion relation for the Gamma function

                    double r = nu2 / nu1;
                    double M = 1.0;
                    for (int k = 0; k < n; k++) {
                        M = M * r * (nu1 + 2*k) / (nu2 - 2*(k+1));
                    }

                    return (M);

                }
            }
		}

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (0.0);
            } else if (n == 2) {
                return (Variance);
            } else if (n == 3) {
                if (nu2 <= 6.0) {
                    return (Double.PositiveInfinity);
                } else {
                    double m = Mean;
                    return (8 * m * m * m * (2 * nu1 + nu2 - 2.0) * (nu1 + nu2 - 2.0) / (nu1 * nu1) / (nu2 - 4.0) / (nu2 - 6.0));
                }
            } else {
                if (nu2 < 2.0 * n) {
                    return (Double.PositiveInfinity);
                } else {
                    return (CentralMomentFromRawMoment(n));
                }
            }
		}
		
	}

}