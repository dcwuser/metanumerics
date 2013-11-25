using System;


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
                return ((n + 1) * (n - 1)/ 12.0);
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
            if ((P < 0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            return (a + (int) Math.Floor(n * P));
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else {
                if (n % 2 != 0) {
                    return (0.0);
                } else {
                    return (base.MomentAboutMean(n));
                }
            }
        }

    }

}
