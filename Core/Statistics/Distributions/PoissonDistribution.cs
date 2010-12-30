using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represented a Poisson distribution.
    /// </summary>
    public class PoissonDistribution : DiscreteDistribution {

        /// <summary>
        /// Creates a new instance of a Poisson distribution.
        /// </summary>
        /// <param name="mu">The mean, which must be positive.</param>
        public PoissonDistribution (double mu) {
            if (mu <= 0.0) throw new ArgumentOutOfRangeException("mu");
            this.mu = mu;
        }

        private double mu;

        /// <inheritdoc />
        public override double ProbabilityMass (int k) {
            if (k < 0) {
                return (0.0);
            } else {
                // these are the same expression, but the form for small arguments is faster and the form for large arguments avoids overflow
                if (k < 16) {
                    return (Math.Exp(-mu) * MoreMath.Pow(mu, k) / AdvancedIntegerMath.Factorial(k));
                } else {
                    return (Math.Exp(
                        k * Math.Log(mu) - AdvancedIntegerMath.LogFactorial(k) - mu
                    ));
                }
            }
        }

        /// <inheritdoc />
        public override DiscreteInterval Support {
            get {
                return (DiscreteInterval.FromEndpoints(0, Int32.MaxValue));
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (mu);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (mu);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (1.0 / Math.Sqrt(mu));
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (int k) {
            if (k < 0) {
                return (0.0);
            } else {
                // these are equivilent expressions, but the explicit sum is faster for if the number of terms is small
                if (k < 16) {
                    double ds = 1.0;
                    double s = ds;
                    for (int i = 1; i <= k; i++) {
                        ds *= mu / i;
                        s += ds;
                    }
                    return (Math.Exp(-mu) * s);
                } else {
                    return(AdvancedMath.RightRegularizedGamma(k + 1, mu));
                }
            }
        }

        /// <inheritdoc />
        public override double RightProbability (int k) {
            if (k < 0) {
                return (1.0);
            } else {
                return (AdvancedMath.LeftRegularizedGamma(k + 1, mu));
            }
        }

        internal override double Cumulant (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (0.0);
            } else {
                return (mu);
            }
        }

        /// <inheritdoc />
        public override int InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");

            //Console.WriteLine("P={0}", P);

            // expand interval until we bracket P
            int ka = 0;
            int kb = (int) Math.Truncate(mu);
            //Console.WriteLine("[{0} {1}] {2}", ka, kb, LeftProbability(kb));
            while (P > LeftProbability(kb)) {
                ka = kb;
                kb = 2 * kb;
                //Console.WriteLine("[{0} {1}] {2}", ka, kb, LeftProbability(kb));
            }

            // reduce interval until we have isolated P
            // this logic is copied from Binomial; we should factor it out
            while (ka != kb) {
                int k = (ka + kb) / 2;
                if (P > LeftProbability(k)) {
                    ka = k + 1;
                } else {
                    kb = k;
                }
            }
            return (ka);

        }

    }

}
