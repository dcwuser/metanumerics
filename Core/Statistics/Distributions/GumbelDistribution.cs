using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Gumbel distribution.
    /// </summary>
    /// <seealso href="http://en.wikipedia.org/wiki/Gumbel_distribution"/>
    public class GumbelDistribution : Distribution {

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(Double.NegativeInfinity, Double.PositiveInfinity));
            }
        }

        /// <summary>
        /// Initializes a new standard Gumbel distribution.
        /// </summary>
        public GumbelDistribution () : this(0.0, 1.0) {
        }

        /// <summary>
        /// Initializes a new Gumbel distribution with the given parameters.
        /// </summary>
        /// <param name="location">The location parameter.</param>
        /// <param name="width">The width parameter, which must be positive.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="width"/> is negative or zero.</exception>
        public GumbelDistribution (double location, double width) {
            if (width <= 0.0) throw new ArgumentOutOfRangeException("width");
            this.m = location;
            this.s = width;
        }

        private double m, s;

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            double z = (x - m) / s;
            double e = Math.Exp(-z);
            // check for infinite e; when e = PositiveInfinity, Exp(-e) = 0, but the computer doesn't recognize that the
            // latter is even smaller than the former is big, and gives NaN rather than zero when they are multiplied
            // for some reason e ~ Infinity rather than PositiveInfinity, so we check for that; i actually thought that
            // infinities were always either PositiveInfinity or NegativeInfinity, but that appears not to be the case
            if (Double.IsInfinity(e)) {
                return (0.0);
            } else {
                return (e * Math.Exp(-e) / s);
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            double z = (x - m) / s;
            double e = Math.Exp(-z);
            return (Math.Exp(-e));
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            double z = (x - m) / s;
            double e = Math.Exp(-z);
            return (-MoreMath.ExpMinusOne(-e));
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            return (m - s * Math.Log(-Math.Log(P)));
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (m - s * Math.Log(Global.LogTwo));
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (m + s * AdvancedMath.EulerGamma);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (MoreMath.Pow2(Math.PI * s) / 6.0);
            }
        }

        /// <inheritdoc />
        public override double StandardDeviation {
            get {
                return (Math.PI * s / Math.Sqrt(3.0));
            }
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 1) {
                return (0.0);
            } else {
                return (AdvancedIntegerMath.Factorial(n - 1) * AdvancedMath.RiemannZeta(n) * MoreMath.Pow(s, n));
            }
        }

        internal override double Cumulant (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 1) {
                return (0.0);
            } else {
                return (AdvancedIntegerMath.Factorial(n - 1) * AdvancedMath.RiemannZeta(n) * MoreMath.Pow(s, n));
            }
        }

    }

}
