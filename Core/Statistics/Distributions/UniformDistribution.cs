using System;

using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {


    /// <summary>
    /// Represents a uniform distribution over an interval.
    /// </summary>
    /// <seealso href="http://en.wikipedia.org/wiki/Uniform_distribution_(continuous)"/>
    public class UniformDistribution : Distribution {

        private Interval range;

        /// <summary>
        /// Gets the range of the uniform distribution.
        /// </summary>
        public Interval Range {
            get {
                return (range);
            }
        }

        /// <summary>
        /// Initializes a new uniform distribution on the given interval.
        /// </summary>
        /// <param name="range">The range of the distribution.</param>
        public UniformDistribution (Interval range) {
            this.range = range;
        }


        /// <summary>
        /// Initializes a new standard uniform distribution.
        /// </summary>
        /// <remarks>
        /// <para>A standard uniform distribution is uniform on the interval [0,1].</para>
        /// </remarks>
        public UniformDistribution () {
            this.range = Interval.FromEndpoints(0.0, 1.0);
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if (range.ClosedContains(x)) {
                return (1.0 / (range.Width));
            } else {
                return (0.0);
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x < range.LeftEndpoint) {
                return (0.0);
            } else if (x > range.RightEndpoint) {
                return (1.0);
            } else {
                return ((x - range.LeftEndpoint) / range.Width);
            }
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            return (range.LeftEndpoint * (1.0 - P) + range.RightEndpoint * P);
        }

        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (Mean);
            } else {

                double m = range.Midpoint;
                double w = range.Width;

                if (Math.Abs(m) > w) {

                    // the midpoint is greater than the width
                    // start from the approximate value m^n and compute a correction factor in terms of (w/m)

                    double f = Math.Pow(m, n);

                    double r = w / m;
                    double rr = 4.0 * r * r;

                    double dg = 1.0;
                    double g = dg;
                    for (int k = 2; k <= n; k += 2) {
                        dg = dg * rr;
                        g += g * AdvancedIntegerMath.BinomialCoefficient(n, k) / (k + 1.0);
                    }

                    return (f * g);

                } else {
                    // the width is large compared to the midpoint
                    // it should be safe to do a simple subtraction of endpoint powers

                    return ((Math.Pow(range.RightEndpoint, n + 1) - Math.Pow(range.LeftEndpoint, n + 1)) / range.Width / (n + 1));
                }
            }
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) throw new ArgumentOutOfRangeException("n");
            if ((n % 2) != 0) {
                return (0.0);
            } else {
                return (Math.Pow(range.Width / 2.0, n) / (n + 1));
            }
        }

        internal override double Cumulant (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (0.0);
            } else if (n == 1) {
                return (range.Midpoint);
            } else {
                // B_n / n where B_n is nth Bernoulli number
                throw new NotImplementedException();
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (range.Midpoint);
            }
        }

        /// <inheritdoc />
        public override double StandardDeviation {
            get {
                return (range.Width / Math.Sqrt(12.0));
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (0.0);
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (Mean);
            }
        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (range);
            }
        }

    }

}