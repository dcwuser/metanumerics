using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents the distribution of Pearsons's r statistic.
    /// </summary>
    public sealed class PearsonRDistribution : Distribution {

        private int n;

        /// <summary>
        /// Initializes a new instance of the Pearson r distribution for the given number of pairs.
        /// </summary>
        /// <param name="n">The number of pairs, which must be three or more.</param>
        public PearsonRDistribution (int n) {
            if (n < 3) throw new ArgumentOutOfRangeException("n");
            this.n = n;
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if (Math.Abs(x) > 1.0) {
                return (0.0);
            } else {
                return (Math.Pow(1.0 - x * x, (n - 4) / 2.0) / AdvancedMath.Beta(0.5, (n - 2) / 2.0));
            }
        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(-1.0, 1.0));
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return(0.0);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (1.0 / (n-1));
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (0.0);
            }
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n % 2 != 0) {
                return (0.0);
            } else {
                double M = 1.0;
                for (int i = 0; i < n; i+= 2) {
                    M = M * (i + 1) / (this.n + i - 1);
                }
                return (M);
            }
        }

        /// <inheritdoc />
        public override double Moment (int n) {
            return (MomentAboutMean(n));
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= -1.0) {
                return (0.0);
            } else if (x < 0.0) {
                return (AdvancedMath.Beta((n - 2) / 2.0, 0.5, 1.0 - x * x) / AdvancedMath.Beta((n-2) / 2.0, 0.5) / 2.0);
            } else if (x < 1.0) {
                return ((1.0 + AdvancedMath.Beta(0.5, (n - 2) / 2.0, x * x) / AdvancedMath.Beta(0.5, (n-2) / 2.0)) / 2.0);
            } else {
                return (1.0);
            }
        }

        // central probability between -x and x is I_{\sqrt{x}}(1/2, a+1)

    }
}
