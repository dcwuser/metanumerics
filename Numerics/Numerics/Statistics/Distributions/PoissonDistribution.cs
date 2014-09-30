using System;
using System.Collections.Generic;

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

        private readonly double mu;

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
        public override double ExcessKurtosis {
            get {
                return (1.0 / mu);
            }
        }

        /// <inheritdoc />
        public override double LeftInclusiveProbability (int k) {
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
        public override double LeftExclusiveProbability (int k) {
            return(LeftInclusiveProbability(k-1));
        }

        /// <inheritdoc />
        public override double RightExclusiveProbability (int k) {
            if (k < 0) {
                return (1.0);
            } else {
                return (AdvancedMath.LeftRegularizedGamma(k + 1, mu));
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
            while (P > LeftInclusiveProbability(kb)) {
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
                if (P > LeftInclusiveProbability(k)) {
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
                double[] M = new double[r + 1];
                M[0] = 1.0;
                M[1] = mu;
                ComputePoissonRawMoments(M);
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
            } else {
                double[] C = new double[r + 1];
                C[0] = 1.0;
                C[1] = 0.0;
                ComputePoissonCentralMoments(C);
                return (C[r]);
            }

        }

        // The raw moments of the Poisson distribution are the Touchard polynomials T_r(x) evaluated at x = \mu.
        // Since the Touchard polynomials satisfy the recurrence
        //   T_{r+1}(x) = x \sum_{s=0}^{r} { r \choose s } T_s(x)
        // The raw moments satisfy the same recurrence:
        //   M_{r+1) = \mu \sum_{s=0}^{r} { r \choose s } M_r
        // Notice that this computation is O(r^2).

        private void ComputePoissonRawMoments (double[] M) {
            for (int r = 2; r < M.Length; r++) {
                IEnumerator<double> binomial = AdvancedIntegerMath.BinomialCoefficients(r - 1).GetEnumerator();
                double t = 0.0;
                for (int s = 0; s < r; s++) {
                    binomial.MoveNext();
                    t += binomial.Current * M[s];
                }
                M[r] = mu * t;
            }
        }

        // Privault, "Generalized Bell polynomials and the combinatorics of Poisson central moments",
        // The Electronic Jounrnal of Combinatorics 2011, derives the recurrence
        //   C_{r+1} = \mu \sum{s=0}^{r-1} { r \choose s } C_s
        // for central moments of the Poisson distribution, directly from the definition.

        // This is almost, but not quite, the same recurrence as for the raw moments. For the raw moment recursion,
        // the bottom binomial argument runs over its full range. For the central moment recursion, the last value of
        // the bottom binomial argument is left out of the range.

        private void ComputePoissonCentralMoments (double[] C) {
            for (int r = 2; r < C.Length; r++) {
                IEnumerator<double> binomial = AdvancedIntegerMath.BinomialCoefficients(r - 1).GetEnumerator();
                double t = 0.0;
                for (int s = 0; s < r - 1; s++) {
                    binomial.MoveNext();
                    t += binomial.Current * C[s];
                }
                C[r] = mu * t;
            }
        }

        /// <inheritdoc />
        public override double Cumulant (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException("r");
            } else if (r == 0) {
                return (0.0);
            } else {
                return (mu);
            }
        }

    }

}
