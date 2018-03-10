using System;

using Meta.Numerics;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Bernoulli distribution.
    /// </summary>
    /// <remarks>
    /// <para>A Bernoulli distribution describes a trial with two possible outcomes. These outcomes are usually called
    /// "success" and "failure", but the same framework is applicable to any binary outcome: right or left, true or false,
    /// male or female, dead or alive, etc.
    /// We represent the outcomes by 0 and 1, which are the only two integers for which the Bernoulli probability
    /// does not vanish.</para>
    /// <para>When multiple, independent Bernoulli trials are conducted, the binomial distribution (<see cref="BinomialDistribution"/>)
    /// describes the probablity of obtaining any particular number of successes.</para>
    /// </remarks>
    /// <seealso cref="BinomialDistribution"/>
    /// <seealso href="http://mathworld.wolfram.com/BinomialDistribution.html"/>
    /// <seealso href="http://en.wikipedia.org/wiki/Bernoulli_distribution"/>
    public sealed class BernoulliDistribution : DiscreteDistribution {

        /// <summary>
        /// Initializes a new Bernoulli distribution.
        /// </summary>
        /// <param name="p1">The probability of outcome 1, which must lie between zero and one.</param>
        public BernoulliDistribution (double p1) {
            if ((p1 < 0.0) || (p1 > 1.0)) throw new ArgumentOutOfRangeException(nameof(p1));
            this.p = p1;
            this.q = 1.0 - p1;
        }

        private readonly double p, q;

        /// <inheritdoc />
        public override double ProbabilityMass (int k) {
            if (k == 0) {
                return (q);
            } else if (k == 1) {
                return (p);
            } else {
                return (0.0);
            }
        }


        /// <inheritdoc />
        public override DiscreteInterval Support {
            get {
                return (new DiscreteInterval(0, 1));
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (p);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (p * q);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return ((q - p) / Math.Sqrt(p * q));
            }
        }

        /// <inheritdoc />
        public override double ExcessKurtosis {
            get {
                return (1.0 / (p * q) - 6.0);
            }
        }

        /// <inheritdoc />
        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else {
                return (p);
            }
        }

        /// <inheritdoc />
        public override double CentralMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else {
                // The mean is p, so C_r = q (0 - p)^r + p (1 - p)^r = q (-p)^r + p q^r.
                // For odd r, there is likely to be significant cancelation between the two opposite-signed terms.
                if (r % 2 == 0) {
                    // For even r, there is no cancelation. Just simplify a bit by pulling out shared factors p and q.
                    //   C_r = p q^r + q p^r = p q ( q^{r-1} + p^{r-1} )
                    return (p * q * (MoreMath.Pow(p, r - 1) + MoreMath.Pow(q, r - 1)));
                } else {
                    // For odd r, we have a difference of powers.
                    //   C_r = p q^r - q p^r = p q ( q^{r-1} - p^{r-1} )
                    // We put our tricks to more accurately evaluate a difference of powers in a seperate method.
                    return (p * q * DifferenceOfPowers(r - 1));
                }
            }
        }

        // Compute q^m - p^m
        private double DifferenceOfPowers (int m) {

            switch (m) {
                case 0:
                    return 0.0;
                case 1:
                    return (q - p);
                case 2:
                    return (q - p) * (q + p);
                case 3:
                    return (q - p) * (q * q + q * p + p * p);
                default:
                    if (m % 2 == 0) {
                        // For even m, use (q^m - p^m) = (q^n - p^n)(q^n + p^n) where m = 2n. 
                        int n = m / 2;
                        return DifferenceOfPowers(n) * (MoreMath.Pow(q, n) + MoreMath.Pow(p, n));
                    } else {
                        // For odd m, we could use (q^m - p^m) = (q - p) (q^{m-1} + p q^{m-2} + \cdots + p^{m-1}),
                        // and we might eventually need to do this to get the accuracy we need. But that is
                        // an awful lot of terms that takes a awful lot of flops and not entirely trivial code
                        // to evaluate. For now, we will accept the occasional difference of powers.
                        return (MoreMath.Pow(q, m) - MoreMath.Pow(p, m));
                        // The first C_r that hits this is C_11, which makes us evaluate q^5 - p^5.
                    }
            }

        }

        /// <inheritdoc />
        public override double LeftExclusiveProbability (int k) {
            if (k <= 0) {
                return (0.0);
            } else if (k == 1) {
                return (q);
            } else {
                return (1.0);
            }
        }

        /// <inheritdoc />
        public override double RightExclusiveProbability (int k) {
            if (k < 0) {
                return (1.0);
            } else if (k == 0) {
                return (p);
            } else {
                return (0.0);
            }
        }

        /// <inheritdoc />
        public override int InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            if (P < q) {
                return (0);
            } else {
                return (1);
            }
        }

        /// <inheritdoc />
        public override double ExpectationValue (Func<int, double> f) {
            if (f == null) throw new ArgumentNullException(nameof(f));
            return (q * f(0) + p * f(1));
        }

        /// <summary>
        /// Finds the Bernoulli distribution that best fits the given counts.
        /// </summary>
        /// <param name="n0">The count of zeros (failure outcomes), which must be non-negative.</param>
        /// <param name="n1">The count of ones (success outcomes), which must be non-negative.</param>
        /// <returns>The best-fit Bernoulli distribution parameters.</returns>
        public static BernoulliFitResult FitToSample (int n0, int n1) {

            if (n0 < 0) throw new ArgumentOutOfRangeException(nameof(n0));
            if (n1 < 0) throw new ArgumentOutOfRangeException(nameof(n1));

            int n = n0 + n1;
            if (n < 2) throw new InsufficientDataException();

            // The only parameter is p. Maximum likelyhood is straightforward.
            //    \ln L = \sum_i \ln P_i = n_0 \ln q + n_1 \ln p = n_0 \ln (1 - p) + n_1 \ln p
            //    \frac{\partial \ln L}{\partial p} = \frac{-n_0}{1-p} + \frac{n_1}{p}
            // Set equal to zero to obtain
            //    p = \frac{n_1}{n_0 + n_1} = \frac{n_1}{n}
            // as would be expected. Continue taking derivatives
            //    \frac{\partial^2 \ln L}{\partial p^2} = \frac{-n_0}{(1-p)^2} - \frac{n_1}{p^2}
            //        = - \frac{(n_0 + n_1)^3}{n_0 n_1} = - \frac{n}{p q}
            // so (\delta p)^2 = \frac{p q}{n}
            
            // When X ~ Bernoulli(p), \sum_{i=1}^{n} X_i ~ Binomial(n, p). We know that Binomial(n, p)
            // has mean n p and variance n p q, so the MLE results are exact and unbiased.

            double p = ((double) n1) / ((double) n);
            double q = ((double) n0) / ((double) n);
            double vp = p * q / n;

            double chi2 = MoreMath.Sqr(n * q - n0) + MoreMath.Sqr(n * p - n1);
            TestResult test = new TestResult("chi2", chi2, TestType.RightTailed, new ChiSquaredDistribution(1));

            return (new BernoulliFitResult(new UncertainValue(p, Math.Sqrt(vp)), test));
        }

    }

}
