using System;
using System.Collections.Generic;

using Meta.Numerics;
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
            if ((p < 0.0) || (p > 1.0)) throw new ArgumentOutOfRangeException("p");
            if (n < 1) throw new ArgumentOutOfRangeException("m");
            this.p = p;
            this.q = 1.0 - p;
            this.n = n;
        }

        private int n;
        private double p, q;

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


        // There are multiple ways to get 

        /// <inheritdoc />
        public override double Moment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (r == 0) {
                return (1.0);
            } else {
                // if m is large, this is slow; find a better way
                return (ExpectationValue(delegate (int k) { return (MoreMath.Pow(k, r)); }));
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
                throw new ArgumentOutOfRangeException("n");
            } else if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (0.0);
            } else {

                // This requires computing r(r-1)/2 binomial coefficients. Summing over all values requires computing n binomial coefficients.

                double[] C = new double[r + 1];
                C[0] = 1.0;
                C[1] = 0.0;

                for (int i = 2; i <= r; i++) {
                    double s = 0.0;
                    for (int j = 0; j <= i - 2; j++) {
                        s += AdvancedIntegerMath.BinomialCoefficient(i - 1, j) * (n * q * C[j] - C[j + 1]);
                    }
                    C[i] = p * s;
                }

                return (C[r]);

            }
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
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            if (n < 16) {

                // for small m, just add up probabilities directly
                // here PP is p^k q^(m-k), PS is the sum of PPs up to k, and r = p / q is used to increment PP
                double r = p / q;
                double PP = MoreMath.Pow(q, n);
                double PS = 0.0;
                /*
                using (IEnumerator<double> B = AdvancedIntegerMath.BinomialCoefficients(m).GetEnumerator()) {
                    // by ending this loop at k = m - 1 instead of k = m, we not only avoid evaluating P(m), but we also
                    // ensure that if PS = 0.999999 due to floating point noise, and we have P = 0.99999999, we dont
                    // return m+1 when we should return m
                    for (int k = 0; k < m; k++) {
                        B.MoveNext();
                        PS += B.Current * PP;
                        if (PS >= P) return(k);
                        PP *= r;
                    }
                    return (m);
                }
                */
                
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
                /*
                // for small distributions, just add probabilities directly
                int k = 0;
                double P0 = MoreMath.Pow(q, m);
                double PP = P0;
                while (k < m) {
                    if (P <= PP) break;
                    k++;
                    P0 *= (m + 1 - k) / k * p / q;
                    PP += P0;
                }
                return (k);
                */
            } else {
                // for larger distributions, use bisection
                // this will require log_{2}(m) CDF evaluations, which is at most 31
                int ka = 0;
                int kb = n;
                while (ka != kb) {
                    int k = (ka + kb) / 2;
                    if (P > LeftExclusiveProbability(k)) {
                        ka = k + 1;
                    } else {
                        kb = k;
                    }
                }
                return (ka);
            }
        }


    }

}
