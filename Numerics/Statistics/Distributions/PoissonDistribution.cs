using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represented a Poisson distribution.
    /// </summary>
    public sealed class PoissonDistribution : DiscreteDistribution {

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

#if PAST
        /// <inheritdoc />
        public override DiscreteInterval Support {
            get {
                return (DiscreteInterval.FromEndpoints(0, Int32.MaxValue));
            }
        }
#endif

        /// <inheritdoc />
        public override int Minimum {
            get {
                return (0);
            }
        }

        /// <inheritdoc />
        public override int Maximum {
            get {
                return (Int32.MaxValue);
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
        public override double LeftExclusiveProbability (int k) {
            if (k < 0) {
                return (0.0);
            } else {
                // these are equivilent expressions, but the explicit sum is faster for if the number of terms is small
                if (k < 16) {
                    double ds = 1.0;
                    double s = ds;
                    for (int i = 1; i < k; i++) {
                        ds *= mu / i;
                        s += ds;
                    }
                    return (Math.Exp(-mu) * s);
                } else {
                    return(AdvancedMath.RightRegularizedGamma(k, mu));
                }
            }
        }

        /// <inheritdoc />
        public override double RightExclusiveProbability (int k) {
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
            //int count = 0;

            // expand interval until we bracket P
            int ka = 0;
            int kb = (int) Math.Ceiling(mu);
            //Console.WriteLine("[{0} {1}] {2}", ka, kb, LeftProbability(kb));
            while (P > LeftExclusiveProbability(kb)) {
                ka = kb;
                kb = 2 * kb;
                //Console.WriteLine("[{0} {1}] {2}", ka, kb, LeftProbability(kb));
                //count++;
                //if (count > 32) throw new NonconvergenceException();
            }

            //count = 0;

            // reduce interval until we have isolated P
            // this logic is copied from Binomial; we should factor it out
            while (ka != kb) {
                int k = (ka + kb) / 2;
                if (P > LeftExclusiveProbability(k)) {
                    ka = k + 1;
                } else {
                    kb = k;
                }
                //if (count > 32) throw new NonconvergenceException();
            }
            return (ka);

        }

        /// <inheritdoc />
        public override double Moment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException("r");
            } else if (r == 0) {
                return (1.0);
            } else {

                // Moments are given by Touchard polynomials M_r = T_r(\mu).
                // Touchard polynomials follow recurrence T_{r+1}(x) = x \sum_{k=0}^{r} {r choose k} T_{k}(x).
                // Basically, this is a hidden way of computing Stirling numbers of the 2nd kind.

                double[] M = new double[r + 1];
                M[0] = 1.0;
                M[1] = mu;

                for (int i = 1; i < r; i++) {
                    double s = 0.0;
                    for (int j = 0; j <= i; j++) {
                        s += AdvancedIntegerMath.BinomialCoefficient(i, j) * M[j];
                    }
                    M[i + 1] = mu * s;
                }

                return (M[r]);

                // If mu is small enough and r large enough, it will be faster to just sum terms.

            }

        }


        /// <inheritdoc />
        public override double MomentAboutMean (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException("r");
            } else if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (0.0);
            } else {

                // Moments are given by Touchard polynomials M_r = T_r(\mu).
                // Touchard polynomials follow recurrence T_{r+1}(x) = x \sum_{k=0}^{r} {r choose k} T_{k}(x).
                // Basically, this is a hidden way of computing Stirling numbers of the 2nd kind.

                double[] C = new double[r + 1];
                C[0] = 1.0;
                C[1] = 0.0;

                for (int i = 1; i < r; i++) {
                    double s = 0.0;
                    for (int j = 0; j < i; j++) {
                        s += AdvancedIntegerMath.BinomialCoefficient(i, j) * C[j];
                    }
                    C[i + 1] = mu * s;
                }

                return (C[r]);

                // If mu is small enough and r large enough, it will be faster to just sum terms.

            }

        }

    }

}
