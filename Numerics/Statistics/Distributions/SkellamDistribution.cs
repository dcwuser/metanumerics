using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Skellam distribution.
    /// </summary>
    /// <remarks>
    /// <para>The difference of two independent Poisson distributed variables has a Skellam distribution.</para>
    /// </remarks>
    /// <seealso href="https://en.wikipedia.org/wiki/Skellam_distribution"/>
    public class SkellamDistribution : DiscreteDistribution {

        /// <summary>
        /// Initializes a new Skellam distribution.
        /// </summary>
        /// <param name="mu1">The mean of the first, positively-contributing variate.</param>
        /// <param name="mu2">The mean of the second, negatively-contributing variate.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="mu1"/> or <paramref name="mu2"/> is non-positive.</exception>
        public SkellamDistribution(double mu1, double mu2) {
            if (mu1 <= 0.0) throw new ArgumentOutOfRangeException(nameof(mu1));
            if (mu2 <= 0.0) throw new ArgumentOutOfRangeException(nameof(mu2));
            this.mu1 = mu1;
            this.mu2 = mu2;
            this.generator1 = PoissonDistribution.GetPoissonGenerator(mu1);
            this.generator2 = PoissonDistribution.GetPoissonGenerator(mu2);
        }

        private readonly double mu1, mu2;

        private readonly IDeviateGenerator<int> generator1, generator2;

        /// <inheritdoc/>
        public override DiscreteInterval Support {
            get {
                return DiscreteInterval.Infinite;
            }
        }

        /// <inheritdoc/>
        public override double ProbabilityMass(int k) {
            return AdvancedMath.Skellam(mu1, mu2, k);
            //return Math.Pow(mu1 / mu2, 0.5 * k) * Math.Exp(-MoreMath.Sqr(Math.Sqrt(mu1) - Math.Sqrt(mu2))) * AdvancedMath.ScaledModifiedBesselI(Math.Abs((double) k), 2.0 * Math.Sqrt(mu1 * mu2));
        }


        /// <inheritdoc/>
        public override double LeftInclusiveProbability(int k) {
            double P;
            // Convert integer to double first to deal with endpoint effects (e.g. -Int32.MinValue is still negative)
            double dk = (double)k;
            if (k < 0) {
                AdvancedMath.Marcum(-dk, Math.Sqrt(2.0 * mu1), Math.Sqrt(2.0 * mu2), out P, out _);
            } else {
                AdvancedMath.Marcum(dk + 1.0, Math.Sqrt(2.0 * mu2), Math.Sqrt(2.0 * mu1), out _, out P);
            }
            return P;
        }

        /// <inheritdoc/>
        public override double RightExclusiveProbability(int k) {
            double Q;
            double dk = (double)k;
            if (k < 0) {
                AdvancedMath.Marcum(-dk, Math.Sqrt(2.0 * mu1), Math.Sqrt(2.0 * mu2), out _, out Q);
            } else {
                AdvancedMath.Marcum(dk + 1.0, Math.Sqrt(2.0 * mu2), Math.Sqrt(2.0 * mu1), out Q, out _);
            }
            return Q;
        }

        /// <inheritdoc/>
        public override double LeftExclusiveProbability(int k) {
            if (k == Int32.MinValue) {
                return 0.0;
            } else {
                return LeftInclusiveProbability(k - 1);
            }
        }

        /// <inheritdoc/>
        public override int InverseLeftProbability(double P) {
            if (P < 0.0 || P > 1.0) throw new ArgumentOutOfRangeException(nameof(P));
            if (P == 0.0) return Support.LeftEndpoint;
            if (P == 1.0) return Support.RightEndpoint;
            // Use Cantelli's inequality to set bounds on how far from the mean the value can be.
            double Q = 1.0 - P;
            int ka = (int)Math.Floor(Mean - Math.Sqrt(Q / P * Variance));
            int kb = (int)Math.Ceiling(Mean + Math.Sqrt(P / Q * Variance));
            return InverseLeftProbability(ka, kb, P);
        }

        /// <inheritdoc/>
        public override double Mean {
            get {
                return mu1 - mu2;
            }
        }

        /// <inheritdoc/>
        public override double Variance {
            get {
                return mu1 + mu2;
            }
        }

        /// <inheritdoc/>
        public override double Skewness {
            get {
                return (mu1 - mu2) / Math.Pow(mu1 + mu2, 3.0 / 2.0);
            }
        }

        /// <inheritdoc/>
        public override double ExcessKurtosis {
            get {
                return 1.0 / (mu1 + mu2);
            }
        }

        /// <inheritdoc/>
        public override double Cumulant(int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return 0.0;
            } else {
                // Even cumulants are sum of mus, odd cumulants are difference of mus.
                // This follows from cumulants adding under convolution and every cumulant of a Poisson distribution being mu.
                if (r % 2 == 0) {
                    return mu1 + mu2;
                } else {
                    return mu1 - mu2;
                }
            }
        }

        /// <inheritdoc/>
        public override double CentralMoment(int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return 1.0;
            } else {
                double[] K = new double[r + 1];
                K[0] = 0.0;
                for (int i = 1; i <= r; i++) {
                    K[i] = (i % 2 == 0) ? mu1 + mu2 : mu1 - mu2;
                }
                return MomentMath.CumulantToCentral(K, r);
            }
        }

        /// <inheritdoc/>
        public override double RawMoment(int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return 1.0;
            } else {
                double[] K = new double[r + 1];
                K[0] = 0.0;
                for (int i = 1; i <= r; i++) {
                    K[i] = (i % 2 == 0) ? mu1 + mu2 : mu1 - mu2;
                }
                return MomentMath.CumulantToRaw(K, r);
            }
        }

        /// <inheritdoc/>
        public override int GetRandomValue(Random rng) {
            if (rng is null) throw new ArgumentNullException(nameof(rng));
            int k1 = generator1.GetNext(rng);
            int k2 = generator2.GetNext(rng);
            return k1 - k2;
        }

    }


}
