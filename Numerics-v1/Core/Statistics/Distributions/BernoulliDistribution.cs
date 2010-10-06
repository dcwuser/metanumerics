using System;

using Meta.Numerics;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a Bernoulli distribution.
    /// </summary>
    /// <remarks>
    /// <para>A Bernoulli distribution describes a trial with two possible outcomes. These outcomes are usually called
    /// "success" and "failure", but the same logic is applicable to any binary set of outcomes: right or left, true or false, etc.
    /// We represent the outcomes by 1 and 0, which are the only two integers for which the Bernoulli probability
    /// does not vanish. The parameter p is the probability of obtaining outcome 1.</para>
    /// <para>When multiple, independent Bernoulli trials are conducted, the binomial distribution (<see cref="BinomialDistribution"/>)
    /// describes the probablity of obtaining any particular number of successes.</para>
    /// </remarks>
    /// <seealso cref="BinomialDistribution"/>
    /// <seealso href="http://mathworld.wolfram.com/BinomialDistribution.html"/>
    /// <seealso href="http://en.wikipedia.org/wiki/Bernoulli_distribution"/>
    public class BernoulliDistribution : DiscreteDistribution {

        /// <summary>
        /// Initializes a new Bernoulli distribution.
        /// </summary>
        /// <param name="p">The probability of success, which must lie between zero and one.</param>
        public BernoulliDistribution (double p) {
            if ((p < 0.0) || (p > 1.0)) throw new ArgumentOutOfRangeException("p");
            this.p = p;
            this.q = 1.0 - p;
        }

        private double p, q;

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
                return (DiscreteInterval.FromEndpoints(0, 1));
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
        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else {
                return (p);
            }
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else {
                // these are simplifications of q (-p)^n + p q^n
                if (n % 2 == 0) {
                    return (p * q * (MoreMath.Pow(p, n - 1) + MoreMath.Pow(q, n - 1)));
                } else {
                    int m = (n - 1) / 2;
                    double pm = MoreMath.Pow(p, m);
                    double qm = MoreMath.Pow(q, m);
                    return (p * q * (qm - pm) * (qm + pm));
                }
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (int k) {
            if (k < 0) {
                return (0.0);
            } else if (k == 0) {
                return (q);
            } else {
                return (1.0);
            }
        }

        /// <inheritdoc />
        public override double RightProbability (int k) {
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
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            if (P < q) {
                return (0);
            } else {
                return (1);
            }
        }


    }

}
