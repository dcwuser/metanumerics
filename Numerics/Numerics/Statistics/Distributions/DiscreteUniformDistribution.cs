using System;

using Meta.Numerics.Functions;


namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Describes a discrete uniform distribution.
    /// </summary>
    /// <remarks>
    /// <para>In a discrete uniform distribution, each integer in the allowed range is equally probable.</para>
    /// <para>For example, the distribution of results for one roll of a fair die is DiscreteUniformDistribution(1,6).</para>
    /// </remarks>
    public sealed class DiscreteUniformDistribution : DiscreteDistribution {

        /// <summary>
        /// Instantiates a new discrete uniform distribution with the given endpoints.
        /// </summary>
        /// <param name="a">One end-point.</param>
        /// <param name="b">The other end-point.</param>
        public DiscreteUniformDistribution (int a, int b) {
            if (b < a) Global.Swap(ref a, ref b);
            this.a = a;
            this.b = b;
            this.n = b - a + 1;
        }

        private int a, b, n;

#if PAST
        /// <inheritdoc />
        public override DiscreteInterval Support {
            get { return (DiscreteInterval.FromEndpoints(a, b)); }
        }
#endif

        /// <inheritdoc />
        public override int Minimum {
            get {
                return (a);
            }
        }

        /// <inheritdoc />
        public override int Maximum {
            get {
                return (b);
            }
        }

        /// <inheritdoc />
        public override double ProbabilityMass (int k) {
            if (k < a) {
                return (0.0);
            } else if (k > b) {
                return (0.0);
            } else {
                return (1.0 / n);
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return ((a + b) / 2.0);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return ((n + 1) * (n - 1) / 12.0);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (0.0);
            }
        }

        /// <inheritdoc />
        public override double LeftExclusiveProbability (int k) {
            if (k <= a) {
                return (0.0);
            } else if (k > b) {
                return (1.0);
            } else {
                return ((k - a) / ((double) n));
            }
        }

        /// <inheritdoc />
        public override double RightExclusiveProbability (int k) {
            if (k < a) {
                return (1.0);
            } else if (k >= b) {
                return (0.0);
            } else {
                return ((b - k) / ((double) n));
            }
        }

        /// <inheritdoc />
        public override int InverseLeftProbability (double P) {
            if ((P < 0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            return (a + (int) Math.Floor(n * P));
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else {
                if (r % 2 != 0) {
                    // Odd central moments are zero by symmetry, so no need to compute them. 
                    return (0.0);
                } else {
                    // Express using Hurwitz Zeta or Bernoulli polynomial
                    return (base.MomentAboutMean(r));
                }
            }
        }

        // Since every k has the same probability, the base class's
        // attempt to sample only ones near the average will fail,
        // and we might as well take advantage of the uniform PMF
        // to only multiply by it once.

        /// <inheritdoc />
        public override double ExpectationValue (Func<int, double> f) {
            if (f == null) throw new ArgumentNullException(nameof(f));
            double s = 0.0;
            for (int k = a; k <= b; k++) {
                s += f(k);
            }
            return (s / n);
        }

        /// <inheritdoc />
        public override double Cumulant (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (0.0);
            } else if (r == 1) {
                return (Mean);
            } else if (r % 2 != 0) {
                return (0.0);
            } else {
                return (AdvancedIntegerMath.BernoulliNumber(r) / r * (MoreMath.Pow(n, r) - 1.0));
            }
        }

    }

}
