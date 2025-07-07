using System;
using System.Collections.Generic;
using System.Text;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Benford distribution.
    /// </summary>
    /// <remarks>
    /// <para>The Benford distribution is the distribution of leading digits
    /// of values with a scale-invariant distribution.</para>
    /// </remarks>
    public sealed class BenfordDistribution : DiscreteDistribution {

        /// <summary>
        /// Initializes a new Benford distribution for base-10 numbers.
        /// </summary>
        public BenfordDistribution () : this (10) { }

        /// <summary>
        /// Initializes a new Benford distribution.
        /// </summary>
        /// <param name="b">The base of the number system.</param>
        public BenfordDistribution (int b) {
            if (b < 2) throw new ArgumentOutOfRangeException(nameof(b));
            this.b = b;
            this.lnb = Math.Log(b);
        }

        private readonly int b;

        private readonly double lnb;

        /// <inheritdoc/>
        public override DiscreteInterval Support {
            get {
                return new DiscreteInterval(1, b - 1);
            }
        }

        /// <inheritdoc/>
        public override double ProbabilityMass (int k) {
            if (k < 1 || k >= b) {
                return 0.0;
            } else {
                return MoreMath.LogOnePlus(1.0 / k) / lnb;
            }
        }

        // We can compute the CDF in closed form:
        // p_1 + p_2 + p_3 + \cdots + p_k
        //   = \log (1 + 1) + \log(1 + 1/2) + \log(1 + 1/3) + \cdots + \log(1 + 1/k) =
        //   = \log ( 2 \times 3/2 \times 4/3 \times \cdots \times (k+1)/k )
        //   = \log (k + 1)

        /// <inheritdoc/>
        public override double LeftExclusiveProbability(int k) {
            if (k < 2) {
                return 0.0;
            } else if (k < b) {
                return Math.Log(k) / lnb;
            } else {
                return 1.0;
            }
        }

        // Random generation d <- floor(b^u)

    }
}
