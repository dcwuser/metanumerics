using System;

using Meta.Numerics;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Bernoulli distribution.
    /// </summary>
    /// <remarks>
    /// <para>A Bernoulli distribution describes a trial with two possible outcomes. These outcomes are usually called
    /// "success" and "failure", but the same framework is applicable to any binary set of outcomes: right or left, true or false, male or female, dead or alive, etc.
    /// We represent the outcomes by 1 and 0, which are the only two integers for which the Bernoulli probability
    /// does not vanish. The parameter p1 is the probability of obtaining outcome 1.</para>
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
            if ((p1 < 0.0) || (p1 > 1.0)) throw new ArgumentOutOfRangeException("p1");
            this.p = p1;
            this.q = 1.0 - p1;
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

#if PAST
        /// <inheritdoc />
        public override DiscreteInterval Support {
            get {
                return (DiscreteInterval.FromEndpoints(0, 1));
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
                return (1);
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
        public override double Moment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else {
                return (p);
            }
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else {
                // these are simplifications of q (-p)^n + p q^n
                if (r % 2 == 0) {
                    return (p * q * (MoreMath.Pow(p, r - 1) + MoreMath.Pow(q, r - 1)));
                } else {
                    int m = (r - 1) / 2;
                    double pm = MoreMath.Pow(p, m);
                    double qm = MoreMath.Pow(q, m);
                    return (p * q * (qm - pm) * (qm + pm));
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

    }

}
