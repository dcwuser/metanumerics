using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represented a Poisson distribution.
    /// </summary>
    /// <remarks>
    /// <para>If events occur independely at uniformly distributed random times,
    /// the Poisson distribution gives the probabability distribution of the number of
    /// events than occur in a given time interval. If on average &#x3BC; events occur
    /// in the time interval the number of actual events that occur is Poisson distributed
    /// with parameter &#x3BC;.</para>
    /// <para>The Poisson distribution has the unusual property that the distribution
    /// of a sum of Poisson distributed values itself has a poisson distribution,
    /// with a mean equal to the sum of the means of the contributing addends.</para>
    /// <para><img src="../images/SumOfPoisson.png" /></para>
    /// </remarks>
    /// <seealso href="https://en.wikipedia.org/wiki/Poisson_distribution"/>
    /// <seealso href="https://mathworld.wolfram.com/PoissonDistribution.html"/>
    public sealed class PoissonDistribution : DiscreteDistribution {

        /// <summary>
        /// Creates a new instance of a Poisson distribution.
        /// </summary>
        /// <param name="mu">The mean, which must be positive.</param>
        public PoissonDistribution (double mu) {
            if (mu <= 0.0) throw new ArgumentOutOfRangeException(nameof(mu));
            this.mu = mu;
            this.poissonGenerator = GetPoissonGenerator(mu);
        }

        internal static IDeviateGenerator<int> GetPoissonGenerator (double mu) {
            Debug.Assert(mu > 0.0);
            if (mu < 4.0) {
                return new PoissonGeneratorMultiplicative(mu);
            } else if (mu < 16.0) {
                return new PoissonGeneratorTabulated(mu);
            } else {
                return new PoissonGeneratorPTRS(mu);
            }
        }

        private readonly double mu;

        private readonly IDeviateGenerator<int> poissonGenerator;

        /// <inheritdoc />
        public override double ProbabilityMass (int k) {
            if (k < 0) {
                return 0.0;
            } else {
                return AdvancedMath.PoissonProbability(mu, k);
            }
        }

        /// <inheritdoc />
        public override DiscreteInterval Support {
            get {
                return DiscreteInterval.Semiinfinite;
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return mu;
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return mu;
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return 1.0 / Math.Sqrt(mu);
            }
        }

        /// <inheritdoc />
        public override double ExcessKurtosis {
            get {
                return 1.0 / mu;
            }
        }

        /// <inheritdoc />
        public override double LeftInclusiveProbability (int k) {
            if (k == Int32.MaxValue) {
                return 1.0;
            } else {
                return LeftExclusiveProbability(k + 1);
            }
        }

        /// <inheritdoc />
        public override double LeftExclusiveProbability (int k) {
            if (k < 1) {
                return 0.0;
            } else if (k < 16) {
                // For small k, it's fastest to sum probabilities directly.
                double ds = 1.0;
                double s = ds;
                for (int i = 1; i < k; i++) {
                    ds *= mu / i;
                    s += ds;
                }
                return Math.Exp(-mu) * s;
            } else {
                // Otherwise, fall back to regularized incomplete Gamma.
                return AdvancedMath.RightRegularizedGamma(k, mu);
            }
        }

        /// <inheritdoc />
        public override double RightExclusiveProbability (int k) {
            if (k < 0) {
                return 1.0;
            } else {
                return AdvancedMath.LeftRegularizedGamma(k + 1, mu);
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

            // If kmax is too big to be sure direct summation would be quick, we
            // will do binary search instead. We want to start with the best limits
            // we can quickly compute to avoid extra evaluations of the CDF.

            // We can use the median bounds
            //   mu - \log 2 \ge \nu \le mu + 1/3
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
        public override double RawMoment (int r) {
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
        public override double CentralMoment (int r) {
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
                return 0.0;
            } else {
                return mu;
            }
        }

        // There is a long and rich history of algorithms for Poisson deviate generation.
        // See Hoermann https://pdfs.semanticscholar.org/00f1/8468dd53e4e3f0c28a5d504f8cd10376a673.pdf for a summary.

        // The simplest generation method besides inversion is to make use the definition of a poisson process and
        // count the number of random numbers needed before their product is less than e^{-\lambda} (or their sum
        // exceeds \lambda). This is fast for very low \lambda (less than ~3), but time (and number of uniform
        // deviates consumed) clearly scales linearly with \lambda, so it rapidly becomes untentable for larger \lambda.

        // Ahrens and Dieter, "Computer generation of Poisson deviates from modified normal distributions",
        // ACM Transactions on Mathematical Software 8 (1982) 163-179 presented the first O(1) generator, PD.
        // The code is complicated and requires normal deviates as inputs.

        // We use Hermann's PTRS generator.

        // NR contains a ratio-of-uniforms variant with no citation and a squeeze I don't recognize. It isn't
        // as fast as PTRS, though.

        // For intermediate value of \lambda (say ~3 to ~30), the fastest method is actually an optimized
        // form of simple inversion: binary search, starting from the median, using a pre-computed table
        // of CDF values that covers most of the space.

        /// <inheritdoc />
        public override int GetRandomValue (Random rng) {
            if (rng is null) throw new ArgumentNullException(nameof(rng));
            return poissonGenerator.GetNext(rng);
        }

    }

    internal class PoissonGeneratorMultiplicative : IDeviateGenerator<int> {

        public PoissonGeneratorMultiplicative (double lambda) {
            Debug.Assert(lambda > 0.0);
            expMinusLambda = Math.Exp(-lambda);
        }

        private readonly double expMinusLambda;

        public int GetNext (Random rng) {
            int k = 0;
            double t = rng.NextDouble();
            while (t > expMinusLambda) {
                k++;
                t *= rng.NextDouble();
            }
            return k;
        }

    }

    internal class PoissonGeneratorTabulated : IDeviateGenerator<int> {

        public PoissonGeneratorTabulated (double lambda) {
            Debug.Assert(lambda > 0.0);

            cdf = new double[(int) Math.Ceiling(4.0 * lambda)];

            double P = Math.Exp(-lambda);
            cdf[0] = P;
            for (int k = 1; k < cdf.Length; k++) {
                P *= lambda / k;
                cdf[k] = cdf[k - 1] + P;
            }
        }

        private readonly double[] cdf;

        public int GetNext (Random rng) {

            double u = rng.NextDouble();

            int k = Array.BinarySearch<double>(cdf, u);
            if (k < 0) k = ~k;

            return k;
        }

    }

    internal class PoissonGeneratorPTRS : IDeviateGenerator<int> {

        public PoissonGeneratorPTRS (double lambda) {
            Debug.Assert(lambda > 0.0);
            this.lambda = lambda;
            this.lnmu = Math.Log(lambda);
            this.b = 0.931 + 2.53 * Math.Sqrt(lambda);
            this.a = -0.059 + 0.02483 * b;
            this.vr = 0.9277 - 3.6224 / (b - 2.0);
        }

        private readonly double lambda, lnmu, b, a, vr;

        public int GetNext (Random rng) {

            while (true) {
                double u = rng.NextDouble() - 0.5;
                double v = rng.NextDouble();

                double us = 0.5 - Math.Abs(u);
                int k = (int) Math.Floor((2.0 * a / us + b) * u + lambda + 0.43);
                if (us > 0.07 && v < vr) return (k);
                if (k < 0) continue;
                if (us < 0.013 && v > us) continue;
                double ai = 1.1239 + 1.1328 / (b - 3.4);
                if (Math.Log(v * ai / (a / (us * us) + b)) <= -lambda + k * lnmu - AdvancedIntegerMath.LogFactorial(k)) {
                    return (k);
                }
            }

        }
    }

}
