using System;
using System.Diagnostics;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a hypergeometric distribution.
    /// </summary>
    /// <remarks>
    /// <para>A <see cref="BinomialDistribution"/> gives the probability of obtaining a given number of successes in a given
    /// number of drawns from an infinite population with a given fraction of successes an failures. (Equivilently, the infinite-population
    /// property is also sometimes expressed as making the draws "with replacement".) A hypergeometric distribution,
    /// by way of contrast, gives the probability of obtaining a given number of successes in a given number of draws
    /// from a finite population containing fixed numbers of successes and failures. Each of the draws in the binomial case is identical
    /// and independent, but in the hypergeometric case the draws are without replacement, so the outcome of previous draws affect the probability
    /// of obtaining success or failure in subsequent draws.</para>
    /// </remarks>
    /// <seealso href="https://en.wikipedia.org/wiki/Hypergeometric_distribution"/>
    /// <seealso href="http://mathworld.wolfram.com/HypergeometricDistribution.html"/>
    public sealed class HypergeometricDistribution : DiscreteDistribution {

        /// <summary>
        /// Initializes a new hypergeometric distribution with the given parameters.
        /// </summary>
        /// <param name="nPopulation">The total population, which must be non-negative.</param>
        /// <param name="nSuccessPopulation">The number of successes in the total population, which must be non-negative.</param>
        /// <param name="nDraws">The number of draws from the population, which must be non-negative.</param>
        public HypergeometricDistribution (int nPopulation, int nSuccessPopulation, int nDraws) {
            if (nPopulation < 0) throw new ArgumentOutOfRangeException(nameof(nPopulation));
            if ((nSuccessPopulation < 0) || (nSuccessPopulation > nPopulation)) throw new ArgumentOutOfRangeException(nameof(nSuccessPopulation));
            if ((nDraws < 0) || (nDraws > nPopulation)) throw new ArgumentOutOfRangeException(nameof(nDraws));
            this.nPopulation = nPopulation;
            this.nSuccessPopulation = nSuccessPopulation;
            this.nFailurePopulation = nPopulation - nSuccessPopulation;
            this.nDraws = nDraws;
        }

        private readonly int nPopulation;

        private readonly int nSuccessPopulation, nFailurePopulation;

        private readonly int nDraws;

        /// <summary>
        /// Gets the total population parameter of the distribution.
        /// </summary>
        public int Population {
            get {
                return (nPopulation);
            }
        }

        /// <summary>
        /// Gets the success population parameter of the distribution.
        /// </summary>
        public int SuccessPopulation {
            get {
                return (nSuccessPopulation);
            }
        }

        /// <summary>
        /// Gets the failure population parameter of the distribution.
        /// </summary>
        public int FailurePopulation {
            get {
                return (nFailurePopulation);
            }
        }

        /// <summary>
        /// Gets the number of draws parameter of the distribution.
        /// </summary>
        public int Draws {
            get {
                return (nDraws);
            }
        }

        /// <inheritdoc/>
        public override DiscreteInterval Support {
            get {
                int min = Math.Max(0, nDraws - nFailurePopulation);
                int max = Math.Min(nDraws, nSuccessPopulation);
                return (new DiscreteInterval(min, max));
            }
        }

        /// <inheritdoc/>
        public override double ProbabilityMass (int k) {
            if ((k < Support.LeftEndpoint) || (k > Support.RightEndpoint)) {
                return (0.0);
            } else {
                return (
                    AdvancedIntegerMath.BinomialCoefficient(nSuccessPopulation, k) /
                    AdvancedIntegerMath.BinomialCoefficient(nPopulation, nDraws) *
                    AdvancedIntegerMath.BinomialCoefficient(nFailurePopulation, nDraws - k)
                );
            }
        }

        /// <inheritdoc/>
        public override double Mean {
            get {
                double successFraction = ((double) nSuccessPopulation) / nPopulation;
                return (nDraws * successFraction);
            }
        }

        /// <inheritdoc/>
        public override double Variance {
            get {
                double successFraction = ((double) nSuccessPopulation) / nPopulation;
                double failureFraction = ((double) nFailurePopulation) / nPopulation;
                return (nDraws * successFraction * failureFraction * (nPopulation - nDraws) / (nPopulation - 1));
            }
        }

        /// <inheritdoc/>
        public override double Skewness {
            get {
                return (
                    Math.Sqrt((nPopulation - 1) / (1.0 * nSuccessPopulation * nFailurePopulation * nDraws * (nPopulation - nDraws))) *
                    (nPopulation - 2 * nSuccessPopulation) *(nPopulation - 2 * nDraws) / (nPopulation - 2)
                );
            }
        }

        /// <inheritdoc/>
        public override double LeftExclusiveProbability (int k) {
            if (k <= Support.LeftEndpoint) {
                return (0.0);
            } else if (k <= Mean) {
                return (LeftProbabilitySum(k - 1));
            } else if (k <= Support.RightEndpoint) {
                return (1.0 - RightProbabilitySum(k));
            } else {
                return (1.0);
            }
        }

        /// <inheritdoc/>
        public override double LeftInclusiveProbability (int k) {
            if (k < Support.LeftEndpoint) {
                return (0.0);
            } else if (k <= Mean) {
                return (LeftProbabilitySum(k));
            } else if (k < Support.RightEndpoint) {
                return (1.0 - RightProbabilitySum(k + 1));
            }  else {
                return (1.0);
            }
        }

        /// <inheritdoc/>
        public override double RightExclusiveProbability (int k) {
            if (k < Support.LeftEndpoint) {
                return (1.0);
            } else if (k < Mean) {
                return (1.0 - LeftProbabilitySum(k));
            } else if (k < Support.RightEndpoint) {
                return (RightProbabilitySum(k + 1));
            } else {
                return (0.0);
            }
        }

        private double LeftProbabilitySum (int k) {
            Debug.Assert((Support.LeftEndpoint <= k) && (k <= Support.RightEndpoint));
            double dP = ProbabilityMass(Support.LeftEndpoint);
            double P = dP;
            for (int j = Support.LeftEndpoint + 1; j <= k; j++) {
                dP *= ((double) (nSuccessPopulation - j + 1)) / j * (nDraws - j + 1) / (nFailurePopulation - (nDraws - j));
                P += dP;
            }
            return (P);
        }

        private double RightProbabilitySum (int k) {
            Debug.Assert((Support.LeftEndpoint <= k) && (k <= Support.RightEndpoint));
            double dP = ProbabilityMass(Support.RightEndpoint);
            double P = dP;
            for (int j = Support.RightEndpoint - 1; j >= k; j--) {
                dP *= ((double) (j + 1)) / (nSuccessPopulation - j) * (nFailurePopulation - (nDraws - j) + 1) / (nDraws - j);
                P += dP;
            }
            return (P);
        }

        // Factorial moments known in closed form
        
    }
}
