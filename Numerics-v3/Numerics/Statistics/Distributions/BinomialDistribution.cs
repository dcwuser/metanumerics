using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a discrete binomial distribution.
    /// </summary>
    /// <remarks>
    /// <para>The binomial distribution gives the probability of obtaining k successes
    /// in n independent Bernoulli trials in which the probability of success in each trial is p.</para>
    /// <para>For a single trial, the binomial distribution reduces to a Bernoulli distribution (<see cref="BernoulliDistribution"/>).</para>
    /// <para>The test statistic for a sign test (<see cref="Sample.SignTest"/>) is distributed according to the Bernoulli distribution.</para>
    /// </remarks>
    /// <seealso href="http://mathworld.wolfram.com/BinomialDistribution.html"/>
    /// <seealso href="http://en.wikipedia.org/wiki/Binomial_distribution"/>
    public sealed class BinomialDistribution : DiscreteDistribution {

        /// <summary>
        /// Initializes a new binomial distribution.
        /// </summary>
        /// <param name="p">The probability of the success in a single trial, which must lie between zero and one.</param>
        /// <param name="n">The number of trials, which must be positive.</param>
        public BinomialDistribution (double p, int n) {
            if ((p < 0.0) || (p > 1.0)) throw new ArgumentOutOfRangeException(nameof(p));
            if (n < 1) throw new ArgumentOutOfRangeException(nameof(n));
            this.p = p;
            this.q = 1.0 - p;
            this.n = n;
        }

        private readonly int n;
        private readonly double p, q;

        /// <inheritdoc />
        public override double ProbabilityMass (int k) {
            if ((k < 0) || (k > n)) {
                return (0.0);
            } else {
                // for small enough integers, use the exact binomial coefficient,
                // which can be evaluated quickly and exactly
                return (AdvancedIntegerMath.BinomialCoefficient(n, k) * MoreMath.Pow(p, k) * MoreMath.Pow(q, n - k));
                // this could fail if the binomial overflows; should we go do factorials or
                // use an expansion around the normal approximation?
            }
        }

#if PAST
        /// <inheritdoc />
        public override DiscreteInterval Support {
            get {
                return (DiscreteInterval.FromEndpoints(0, m));
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
                return (n);
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (n * p);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (n * p * q);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return ((q - p) / Math.Sqrt(n * p * q));
            }
        }


        /// <inheritdoc />
        public override double Moment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else {

                if (n < 2 * r) {

                    // If n is small, just do a weighted sum over the support.
                    return (ExpectationValue(delegate (int k) { return (MoreMath.Pow(k, r)); }));

                } else {

                    // If n >> r, convert analytically known factorial moments to raw moments
                    // using Sterling numbers of the 2nd kind. The (falling) factorial moment
                    //   F_r = (n)_r p^r
                    // where (n)_r is a falling factorial.

                    double[] s = AdvancedIntegerMath.StirlingNumbers2(r);
                    double t = 1.0;
                    double M = 0.0;
                    for (int k = 1; k <= r; k++) {
                        // s[0] = 0 (except for r=0, in which case we have already returned above)
                        t *= (n - (k - 1)) * p;
                        M += s[k] * t;
                    }
                    return (M);
                }
                
            }
        }

        // The definition of the central moment produces:
        //   C_{r} = \sum{k=0}^{n} (k - np)^r { n \choose k } p^k q^{n-k}
        // This requires a large number of terms when n is large.

        // There is a well-known recurrsion:
        //  C_{r+1} = p q [ n r C_{r-1} + \frac{d}{dp} C_{r} ]
        // but this requires symbolic differentiation and produces polynomials whose evaluation may be subject to cancelation.

        // A numerically useful recurrsion:
        //  C_{r} = \sum_{k=0}^{r-2} { r-1 \choose k } ( n p q C_{k} - p C_{k+1} )
        // is derived in Shenton, Bowman, & Sheehan, "Sampling Moments of Moments Associated with Univariate Distributions",
        // Journal of the Royal Statistical Society B, 33 (1971) 444

        /// <inheritdoc />
        public override double MomentAboutMean (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (0.0);
            } else {
                double[] C = CentralMoments(r);
                return (C[r]);
            }
        }

        internal override double[] CentralMoments (int rMax) {
 
            double[] C = new double[rMax + 1];
            C[0] = 1.0;

            if (rMax == 0) return (C);

            C[1] = 0.0;

            for (int r = 2; r <= rMax; r++) {
                double s = 0.0;
                IEnumerator<double> B = AdvancedIntegerMath.BinomialCoefficients(r - 1).GetEnumerator();
                for(int k = 0; k <= r - 2; k++) {
                    B.MoveNext();
                    s += B.Current * (n * q * C[k] - C[k + 1]);
                }
                C[r] = p * s;
            }

            return (C);

        }

        // for any m larger than a few, this is a slow way to compute explicit moment
        // Riordan gives a recurrence; we should probably use it, but its implementation is a little complicated

        /// <inheritdoc />
        public override double LeftExclusiveProbability (int k) {
            if (k <= 0) {
                return (0.0);
            } else if (k > n) {
                return (1.0);
            } else {
                // use direct summation in tails
                return (AdvancedMath.LeftRegularizedBeta(n - k + 1, k, q));
            }
        }

        /// <inheritdoc />
        public override double RightExclusiveProbability (int k) {
            if (k < 0) {
                return (1.0);
            } else if (k >= n) {
                return (0.0);
            } else {
                // use direct summation in tails
                return (AdvancedMath.Beta(k + 1, n - k, p) / AdvancedMath.Beta(k + 1, n - k));
            }
        }

        /// <inheritdoc />
        public override int InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            if (n < 16) {

                // for small m, just add up probabilities directly
                // here PP is p^k q^(m-k), PS is the sum of PPs up to k, and r = p / q is used to increment PP
                double r = p / q;
                double PP = MoreMath.Pow(q, n);
                double PS = 0.0;
                
                int k = -1;
                foreach (double B in AdvancedIntegerMath.BinomialCoefficients(n)) {
                    k++;
                    PS += B * PP;
                    if (PS >= P) return(k);
                    PP *= r;
                }
                // we "shouldn't" reach here, but if floating point jitter makes all the values add up to 0.999999 instead of 1,
                // and we have P = 0.99999999, it could happen, so in this case we return m
                return (n);

            } else {

                double Q = 1.0 - P;
                if (Q == 0.0) return (n);
                double mu = n * p;

                int kmin, kmax;
                if (P < 0.5) {
                    kmin = 0;
                    kmax = (int) Math.Ceiling(mu);
                    int kChernov = (int) Math.Floor(mu - Math.Sqrt(-2.0 * mu * Math.Log(P)));
                    if (kChernov > kmin) kmin = kChernov;
                } else {
                    kmin = (int) Math.Floor(mu);
                    kmax = n;
                    int kMarkov = (int) Math.Ceiling(mu / Q);
                    if (kMarkov < kmax) kmax = kMarkov;
                    int kCantelli = (int) Math.Ceiling(mu + Math.Sqrt(mu * q * P / Q));
                    if (kCantelli < kmax) kmax = kCantelli;
                    int kChernov = (int) Math.Ceiling(mu + Math.Sqrt(-2.0 * mu * Math.Log(Q)));
                    if (kChernov < kmax) kmax = kChernov;
                }

                return (InverseLeftProbability(kmin, kmax, P));

            }
        }


    }

}
