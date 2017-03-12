using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Laplace distribution.
    /// </summary>
    public sealed class LaplaceDistribution : Distribution {

        /// <summary>
        /// Initializes a new standard Laplace distribution.
        /// </summary>
        public LaplaceDistribution () : this(0.0, 1.0) { }

        /// <summary>
        /// Initializes a new Laplace distribution with the given location and scale parameters.
        /// </summary>
        /// <param name="location">The location parameter.</param>
        /// <param name="scale">The scale parameter, which must be positive.</param>
        public LaplaceDistribution (double location, double scale) {
            if (scale <= 0.0) throw new ArgumentOutOfRangeException(nameof(scale));
            this.a = location;
            this.b = scale;
        }

        // Location parameter
        private readonly double a;

        // Scale paramter
        private readonly double b;

        /// <inheritdoc/>
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(Double.NegativeInfinity, Double.PositiveInfinity));
            }
        }

        /// <inheritdoc/>
        public override double ProbabilityDensity (double x) {
            return (Math.Exp(-Math.Abs(x - a) / b) / (2.0 * b));
        }

        /// <inheritdoc/>
        public override double LeftProbability (double x) {
            double z = (x - a) / b;
            if (z < 0.0) {
                return (Math.Exp(z) / 2.0);
            } else {
                return (1.0 - Math.Exp(-z) / 2.0);
            }
        }

        /// <inheritdoc/>
        public override double RightProbability (double x) {
            double z = (x - a) / b;
            if (z < 0.0) {
                return (1.0 - Math.Exp(z) / 2.0);
            } else {
                return (Math.Exp(-z) / 2.0);
            }
        }

        /// <inheritdoc/>
        public override double Median {
            get {
                return (a);
            }
        }

        /// <inheritdoc/>
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            if (P <= 0.5) {
                return (a + b * Math.Log(2.0 * P));
            } else {
                throw new NotImplementedException();
            }
        }

        /// <inheritdoc/>
        public override double Mean {
            get {
                return (a);
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation {
            get {
                return (Global.SqrtTwo * b);
            }
        }

        /// <inheritdoc/>
        public override double Variance {
            get {
                return (2.0 * b * b);
            }
        }

        /// <inheritdoc/>
        public override double Skewness {
            get {
                return (0.0);
            }
        }

        /// <inheritdoc/>
        public override double ExcessKurtosis {
            get {
                return (3.0);
            }
        }

        /// <inheritdoc/>
        public override double MomentAboutMean (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r % 2 == 0) {
                return (AdvancedIntegerMath.Factorial(r) / Math.Pow(b, r));
            } else {
                return (0.0);
            }
        }

    }

}
