using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a negative binomial distribution.
    /// </summary>
    /// <remarks>
    /// <para>Consider a series of independent, repeated Bernoulli trials, each of which results in success with probability p or failure
    /// with probability 1-p. If one repeats the trials until r failures occur, the negative binomial distribution gives the probability
    /// of having seen k successes before the rth failure.</para>
    /// <para>Keep in mind that there are several different conventions for the meaning of r, p, and k. The most common are that
    /// k denotes the number of successes before the rth failure (used here and in the referenced Wikipedia article), and that
    /// k denotes the number of failures before the rth success (used by Mathematica). Since these two conventions simply
    /// switch success and failure, they can be interconverted by interchanging p and 1-p.</para>
    /// </remarks>
    /// <seealso cref="BernoulliDistribution"/>
    /// <seealso href="https://en.wikipedia.org/wiki/Negative_binomial_distribution"/>
    /// <seealso href="http://mathworld.wolfram.com/NegativeBinomialDistribution.html"/>
    public sealed class NegativeBinomialDistribution : DiscreteDistribution {


        /// <summary>
        /// Initializes a new negative binomial distribution.
        /// </summary>
        /// <param name="r">The number of failures required, which must be positive.</param>
        /// <param name="p">The probability of success for each trial, which must lie between 0 and 1.</param>
        public NegativeBinomialDistribution (double r, double p) {
            if (r <= 0.0) throw new ArgumentOutOfRangeException(nameof(r));
            if ((p < 0.0) || (p > 1.0)) throw new ArgumentOutOfRangeException(nameof(p));
            this.r = r;
            this.p = p;
            this.q = 1.0 - p;
        }

        private readonly double r;
        private readonly double p, q;

        /// <inheritdoc />
        public override DiscreteInterval Support {
            get {
                return DiscreteInterval.Semiinfinite;
            }
        }

        /// <inheritdoc />
        public override double ProbabilityMass (int k) {
            if (k < 0) {
                return (0.0);
            } else {
                return (
                    //AdvancedIntegerMath.BinomialCoefficient(k + r - 1, k) *
                    MoreMath.Pow(p, k) * Math.Pow(q, r) /
                    AdvancedMath.Beta(r, k + 1) / (k + r)
                );
                // note: n \choose k = \frac{1}{(n+1) B(n - k + 1, k + 1)}
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return p * r / q;
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
            } else {
                return AdvancedMath.LeftRegularizedBeta(r, k, q);
            }
        }

        /// <inheritdoc />
        public override double RightExclusiveProbability (int k) {
            if (k < 0) {
                return 1.0;
            } else {
                return AdvancedMath.LeftRegularizedBeta(k + 1, r, p);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return p * r / (q * q);
            }
        }

        /// <inheritdoc />
        public override double StandardDeviation {
            get {
                return Math.Sqrt(p * r) / q;
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (1.0 + p) / Math.Sqrt(p * r);
            }
        }

        /// <inheritdoc />
        /// <param name="k">The order of moment to compute.</param>
        public override double RawMoment (int k) {

            if (k < 0) {
                throw new ArgumentOutOfRangeException(nameof(k));
            } else if (k == 0) {
                return 1.0;
            } else {

                // The (falling) factorial moments of the negative binomial distribution
                // are F_k = (p/q)^m (r)_m, where (r)_m is a rising factorial.

                // We can use this plus the Stirling numbers, which convert factorial
                // to raw moments.

                double[] s = AdvancedIntegerMath.StirlingNumbers2(k);
                double f = p / q;
                double t = 1.0;
                double M = 0.0;
                for (int i = 1; i <= k; i++) {
                    t *= f * (r + i - 1);
                    M += s[i] * t;
                }
                return M;

            }
        }

        /// <inheritdoc />
        public override int InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));

            double Q = 1.0 - P;
            if (Q == 0.0) return (Int32.MaxValue);

            // Find the Markov upper bound on k.
            int kmax = (int) Math.Min(Math.Ceiling(this.Mean / Q), Int32.MaxValue);

            // If kmax is small enough, we should use direct summation.

            // We should use the median bounds.
            // Payton, Young, and Young, "Bounds for the Difference between Median and Mean
            // of Beta and Negative Binomial Distributions", Metrika 36 (1989) 346-354
            // In our notation, if p > 1/2
            //   \mu - p / q \le m le \mu
            // and if p < 1/2, \mu - 1 \le m \le \mu.

            return (InverseLeftProbability(0, kmax, P));

        }

    }
}
