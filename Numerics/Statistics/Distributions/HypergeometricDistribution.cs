using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a hypergeometric distribution.
    /// </summary>
    /// <remarks>
    /// <para>A <see cref="BinomialDistribution"/> gives the probability of obtaining a given number of successes in a given
    /// number of drawns from an infinite population with a given fraction of successes an failures. A hypergeometric distribution,
    /// by way of contrast, gives the probability of obtaining a given number of successes in a given number of draws
    /// from a finite population with a given number of successes and failures. Each of the draws in the binomial case is identical
    /// and independent, but in the hypergeometric case the the outcome of previous draws affect the probability of obtaining
    /// success or failure in subsequent draws. This difference is often characterized by saying that while the binomial
    /// distribution describes draws with replacement (the obtained result is placed back into the population for the next draw),
    /// the hypergeometic distribution describes draws without replacement.</para>
    /// </remarks>
    public sealed class HypergeometricDistribution : DiscreteDistribution {

        /// <summary>
        /// Initializes a new hypergeometric distribution with the given parameters.
        /// </summary>
        /// <param name="nPopulation"></param>
        /// <param name="nSuccessPopulation"></param>
        /// <param name="nDraws"></param>
        public HypergeometricDistribution (int nPopulation, int nSuccessPopulation, int nDraws) {
            if (nPopulation < 0) throw new ArgumentOutOfRangeException(nameof(nPopulation));
            if ((nSuccessPopulation < 0) || (nSuccessPopulation > nPopulation)) throw new ArgumentOutOfRangeException(nameof(nSuccessPopulation));
            if ((nDraws < 0) || (nDraws > nPopulation)) throw new ArgumentOutOfRangeException(nameof(nDraws));
            this.nPopulation = nPopulation;
            this.nSuccessPopulation = nSuccessPopulation;
            this.nDraws = nDraws;
        }

        private readonly int nPopulation;

        private readonly int nSuccessPopulation;

        private readonly int nDraws;

        public override int Minimum {
            get {
                // Generally, but if all failures have been drawn, we have to start drawing some successes
                int nFailurePopulation = nPopulation - nSuccessPopulation;
                return (Math.Max(0, nDraws - nFailurePopulation));
            }
        }

        public override int Maximum {
            get {
                return (Math.Min(nSuccessPopulation, nSuccessPopulation));
            }
        }

        /// <inheritdoc/>
        public override double ProbabilityMass (int k) {
            if ((k < this.Minimum) || (k > this.Maximum)) {
                return (0.0);
            } else {
                int nFailurePopulation = nPopulation - nSuccessPopulation;
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
                int nFailurePopulation = nPopulation - nSuccessPopulation;
                double successFraction = ((double) nSuccessPopulation) / nPopulation;
                double failureFraction = ((double) nFailurePopulation) / nPopulation;
                return (nDraws * successFraction * failureFraction * (nPopulation - nDraws) / (nPopulation - 1));
            }
        }

        // Factorial moments known in closed form
        // Raw moments from derivatives of 2F1?
        
    }
}
