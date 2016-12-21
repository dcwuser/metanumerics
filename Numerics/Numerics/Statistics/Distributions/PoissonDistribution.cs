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
            if (mu <= 0.0) throw new ArgumentOutOfRangeException(nameof(mu));
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
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));

            // We will need Q = 1 - P, and if it's zero we know we are at infinity.
            double Q = 1.0 - P;
            if (Q == 0.0) return (Int32.MaxValue);

            // Just adding up PMF values is fastest if we can we sure there won't be
            // too many of them, particularly since they obey a recursion relation
            // that requires only a few flops to iterate. But we don't want to get
            // into a situation where we might have to add up 100s or 1000s; that would
            // take too long and accumulate too many errors. We can use the Markov
            // limit on the tail to bound k.
            int kmax = (int) Math.Min(Math.Ceiling(mu / Q), Int32.MaxValue);

            if (kmax < 32) {
                int k = 0;
                double pmf = Math.Exp(-mu);
                double pmfSum = pmf;
                while (pmfSum < P) {
                    k++;
                    pmf *= mu / k;
                    pmfSum += pmf;
                }
                return (k);
            }

            // If kmax is too big to be sure direct sumation would be quick, we
            // will do binary search instead. We want to start with the best limits
            // we can quickly compute to avoid extra evaluations of the CDF.

            // We can use the median bounds
            //   mu - \log 2 \e \nu \le mu + 1/3
            // to quickly set a limit based on whether P is below or above 1/2.
            int kmin;
            if (P < 0.5) {
                kmin = 0;
                kmax = (int) Math.Ceiling(mu + 1.0 / 3.0);

                // Use the Chernov inequality to improve the lower bound, if possible
                int kChernov = (int) Math.Floor(mu - Math.Sqrt(-2.0 * mu * Math.Log(P)));
                if (kChernov > kmin) kmin = kChernov;
                // This isn't really useful if P > 1/2 because the median already gives a lower bound ~\mu
            } else {
                kmin = Math.Max((int) Math.Floor(mu - Global.LogTwo), 0);

                // Use the Cantelli inequality to improve the upper bound, if possible
                int kCantelli = (int) Math.Ceiling(mu + Math.Sqrt(mu * P / Q));
                if (kCantelli < kmax) kmax = kCantelli;
                // This isn't really useful if P < 1/2 because the median already gives an upper bound ~\mu
            }

            return (InverseLeftProbability(kmin, kmax, P));

        }

        /// <inheritdoc />
        public override double Moment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
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
                throw new ArgumentOutOfRangeException(nameof(r));
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

        // The (falling) factorial moments of the Poisson distribution are F_r = \mu^r
        // The raw moments in terms of the factorial moments are
        //   M_r = \sum_{s=0}^r \{ r \over s \} F_s = \sum_{s=0}^r \{ r \over \}s \mu^s = T_r(\mu)
        // These are the Touchard polynomials T_r(x) evaluated at x = \mu.

        // From the Sterling number recurrence
        //   \{ n + 1 \over k + 1 \} = \sum_{j = k}^{n} { n \choose j } \{ j \over k \}
        // follows the Touchard polynomial recurrence
        //   T_{r + 1}(x) = x \sum_{s = 0}^{r} { r \choose s } T_{s}(x)
        // We can use this to compute M_{r}, but note it is a O(r^2) operation and requires O(r)
        // storage, since to compute each new value requires all the previous values.

        // From another Sterling number recurrence
        //   \{ n + 1 \over k \} = k \{ n \over k \} + \{ n \over k - 1 \}
        // follows another Touchard polynomial recurrence

        // The raw moments of the Poisson distribution are the Touchard polynomials T_r(x) evaluated at x = \mu.
        // Since the Touchard polynomials satisfy the recurrence
        //   T_{r+1}(x) = x \sum_{s=0}^{r} { r \choose s } T_s(x)
        // The raw moments satisfy the same recurrence:
        //   M_{r+1} = \mu \sum_{s=0}^{r} { r \choose s } M_r
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
        // The Electronic Jounrnal of Combinatorics 2011
        // (http://www.ntu.edu.sg/home/nprivault/papers/central_moments.pdf),
        // derives the recurrence
        //   C_{r+1} = \mu \sum{s=0}^{r-1} { r \choose s } C_s
        // for central moments of the Poisson distribution, directly from the definition.

        // For the record, here is the derivation, which I have not found elsewhere. By
        // Poisson PMF, mean, and definition of central moment:
        //   C_{r+1} = e^{-\mu} \sum_{k=0}^{\infty} \frac{\mu^k}{k!} (k - \mu)^{r+1}
        // Expand one factor of (k - \mu) to get
        //   C_{r+1} = e^{-\mu} \sum_{k=1}^{\infty} \frac{\mu^k (k - \mu)^{r}}{(k-1)!} 
        //            -e^{-\mu} \sum_{k=0}^{\infty} \frac{\mu^{k+1} (k - \mu)^{r}}{k!}
        // Note k=0 does not contribute to first term because of multiplication by k.
        // Now in first term, redefine dummy summation variable k -> k + 1 to get
        //   C_{r+1} = e^{-\mu} \sum_{k=0}^{\infty} \frac{\mu^{k+1}}{k!}
        //             \left[ (k + 1 - \mu)^{r} - (k - \mu)^{r} \right]
        // Now use the binomial theorem to expand (k - \mu + 1) in factors of (k -\mu) and 1.
        //   C_{r+1} = e^{-\mu} \sum_{k=0}^{\infty} \frac{\mu^{k+1}}{k!}
        //             \sum_{s=0}^{r-1} { r \choose s } (k - \mu)^{s}
        // The s=r term is canceled by the subtracted (k - \mu)^{r}, so the binomial sum only goes to s=r-1.
        // Switch the order of the sums and notice \sum_{k} \frac{\mu^{k}}{k!} (k - \mu)^{s} = C_{s}
        // to obtain the result. Wow.  

        // This is almost, but not quite, the same recurrence as for the raw moments.
        // For the raw moment recursion, the bottom binomial argument runs over its full range.
        // For the central moment recursion, the final value of the bottom binomial argument is left out.
        // Also, for raw moments start the recursion with M_0 = 1, M_1 = \mu. For central moments,
        // start the recursion with C_0 = 1, C_1 = 0.

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
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (0.0);
            } else {
                return (mu);
            }
        }

    }

}
