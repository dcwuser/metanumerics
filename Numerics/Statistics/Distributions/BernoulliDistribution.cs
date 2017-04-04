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
                // Mean is p, so C_r = q (0 - p)^r + p (1 - p)^r = q (-p)^r + p q^r.
                // For odd r, there is likely to be significant cancelation between the two opposite-signed terms.
                // These are simplifications of q (-p)^r + p q^r
                if (r % 2 == 0) {
                    // For even r, there is no cancelation. Just simplify a bit by pulling out shared factors p and q.
                    return (p * q * (MoreMath.Pow(p, r - 1) + MoreMath.Pow(q, r - 1)));
                } else {
                    // For odd r, we have a
                    return (p * q * DifferenceOfPowers(r - 1));
                    /* 
                    int m = (r - 1) / 2;
                    double pm = MoreMath.Pow(p, m);
                    double qm = MoreMath.Pow(q, m);
                    return (p * q * (qm - pm) * (qm + pm));
                    */
                }
            }
        }

        // Compute q^m - p^m
        private double DifferenceOfPowers (int m) {
            if (m == 0) {
                return (0.0);
            } else if (m % 2 == 0) {
                // If m is even, write q^{2n} - p^{2n} = (q^n - p^n) (q^n + p^n).
                int n = m / 2;
                return (DifferenceOfPowers(n) * (MoreMath.Pow(q, n) + MoreMath.Pow(p, n)));
            } else {
                // If m is odd, write q^{n + 1} - p^{n + 1} = (p - q) (p^n + p^{n-1} q + p^{n-2} q^2 + \cdots + q^n)
                // = (p - q) p^n (1 + r + r^2 + \cdots + r^n) where r = q/p.
                // = (p - q) p^n (1 - r^{n + 1}) / (1 - r) = (p - q) p^{n+1} (1 - (q/p)^{n+1})/(p-q)
                double s = MoreMath.Pow(q, m) - MoreMath.Pow(p, m);
                return (s);  
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
