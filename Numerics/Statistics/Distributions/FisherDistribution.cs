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
    /// <para>The Fisher distribution is related to the Beta distribution (<see cref="BetaDistribution"/>) by a simple
    /// variable transformation.</para>
    /// <img src="../images/FisherBeta.png" />
    /// </remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/F_distribution"/>
    public sealed class FisherDistribution : Distribution {

        /// <summary>
        /// Initializes a new Fisher distribution.
        /// </summary>
        /// <param name="nu1">The number of degrees of freedom in the numerator, which must be positive.</param>
        /// <param name="nu2">The number of degrees of freedom in the denominator, which must be positive.</param>
        /// <seealso cref="Sample.FisherFTest"/>
        public FisherDistribution (double nu1, double nu2) {
            if (nu1 <= 0.0) throw new ArgumentOutOfRangeException("nu1");
            if (nu2 <= 0.0) throw new ArgumentOutOfRangeException("nu2");
            this.nu1 = nu1;
            this.nu2 = nu2;
            this.beta = new BetaDistribution(nu1 / 2.0, nu2 / 2.0);
        }

        private double nu1;
        private double nu2;
        private BetaDistribution beta;

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

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else {
                double p = nu1 * x;
                double q = nu2 + p;
                double y = p / q;
                double u = nu1 * nu2 / (q * q);
                return (u * beta.ProbabilityDensity(y));
                /*
                double N = Math.Pow(nu1, 0.5 * nu1) * Math.Pow(nu2, 0.5 * nu2) /
                    AdvancedMath.Beta(0.5 * nu1, 0.5 * nu2);
                return (N * Math.Pow(x, 0.5 * nu1 - 1.0) * Math.Pow(nu2 + nu1 * x, -0.5 * (nu1 + nu2)));
                 */
            }
		}

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else {
                double p = nu1 * x;
                double q = nu2 + p;
                double y = p / q;
                return (beta.LeftProbability(y));
                //return (AdvancedMath.Beta(nu1 / 2.0, nu2 / 2.0, nu1 * x / (nu2 + nu1 * x)) / AdvancedMath.Beta(0.5 * nu1, 0.5 * nu2));
            }
		}

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x <= 0.0) {
                return (1.0);
            } else {
                double p = nu1 * x;
                double q = nu2 + p;
                double y = p / q;
                return (beta.RightProbability(y));
                //return (AdvancedMath.Beta(0.5 * nu2, 0.5 * nu1, nu2 / (nu2 + nu1 * x)) / AdvancedMath.Beta(0.5 * nu2, 0.5 * nu1));
            }
		}

        // moments are not directly calculable from relationship to Beta distribution
        // is there some way to do this that I don't know about?

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

        // The rth central moment is given by the hypergeometric function.
        //   C_r = \left( \frac{m}{2-m} \right)^r _2F_1 ( n/2, -r; 1 - m/2; (2-m)/n )
        // Using the known recursion on the 2nd argument of the hypergeometric function (A&S ##), we can derive the recursion
        //   C_{r+1} = \frac{2 r m}{(m-2)(m-2(r+1))n} \left[ \frac{m(n+m-2)}{(m-2)} C_{r-1} + (2n + m - 2) C_r \right]
        // We implement this recursion here.

        /// <inheritdoc />
        public override double MomentAboutMean (int r) {

            if (r < 0) {
                throw new ArgumentOutOfRangeException("r");
            } else if (r == 0) {
                return (1.0);
            } else {

                if (2 * r >= nu2) return (Double.NaN);

                double CM = 1.0;
                double C = 0.0;

                for (int k = 1; k < r; k++) {
                    double CP = 2 * k * nu2 / ((nu2 - 2.0) * (nu2 - 2.0 * (k + 1)) * nu1) *
                        (nu2 * (nu1 + nu2 - 2.0) / (nu2 - 2.0) * CM + (2.0 * nu1 + nu2 - 2.0) * C);
                    CM = C;
                    C = CP;
                }

                return (C);
            }
        }


        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            double y = beta.InverseLeftProbability(P);
            return(nu2 / nu1 * y / (1.0 - y));
        }

        /// <inheritdoc />
        public override double InverseRightProbability (double Q) {
            double y = beta.InverseRightProbability(Q);
            return (nu2 / nu1 * y / (1.0 - y));
        }
		
	}

}