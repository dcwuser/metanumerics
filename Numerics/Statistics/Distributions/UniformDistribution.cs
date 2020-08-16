using System;

using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {


    /// <summary>
    /// Represents a uniform distribution over an interval.
    /// </summary>
    /// <seealso href="http://en.wikipedia.org/wiki/Uniform_distribution_(continuous)"/>
    /// <seealso href="https://mathworld.wolfram.com/UniformDistribution.html"/>
    public sealed class UniformDistribution : ContinuousDistribution {

        private readonly Interval range;

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
                return (1.0 / range.Width);
            } else {
                return (0.0);
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= range.LeftEndpoint) {
                return (0.0);
            } else if (x >= range.RightEndpoint) {
                return (1.0);
            } else {
                return ((x - range.LeftEndpoint) / range.Width);
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x <= range.LeftEndpoint) {
                return (1.0);
            } else if (x >= range.RightEndpoint) {
                return (0.0);
            } else {
                return ((range.RightEndpoint - x) / range.Width);
            }
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            return (InverseProbability(P, 1.0 - P));
        }

        /// <inheritdoc />
        public override double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException(nameof(Q));
            return(InverseProbability(1.0 - Q, Q));
        }

        private double InverseProbability (double P, double Q) {
            return (range.LeftEndpoint * Q + range.RightEndpoint * P);
        }

        /// <inheritdoc />
        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else {
                // By simple integration, M_r = (b^{r+1} - a^{r + 1}) / (b - a) / (r + 1).
                // We have some relatively complex machinery for comuting this without cancelation error.
                return (DifferenceOfPowers(range.RightEndpoint, range.LeftEndpoint, r + 1) / (r + 1));
            }
        }

        // returns (a^n - b^n), modulo a factor (a - b), computed so as to avoid cancelation errors

        private static double DifferenceOfPowers (double a, double b, int n) {
            if (n % 2 == 0) {
                // for even powers, cancelation can be a problem
                // use a^{2m} - b^{2m} = ( a^m - b^m ) (a^m + b^m)
                int m = n / 2;
                return (DifferenceOfPowers(a, b, m) * SumOfPowers(a, b, m));
            } else {
                // for odd powers, can be a problem only if a and b have the same sign
                // use a^n - b^n = (a - b) (a^{n-1} + a^{n-2} b + a^{n-3} b^2 + \cdots + a b^{n-2} + b^{n-1})
                // All terms in second factor have coefficient 1.
                double[] aPowers = Powers(a, n);
                double[] bPowers = Powers(b, n);
                double s = 0.0;
                for (int m = 0; m <= (n - 1); m++) {
                    s += aPowers[m] * bPowers[(n - 1) - m];
                }
                // We leave out the factor (a - b). Since this block will execute once and only once for any n,
                // this will give us (a^n - b^n) / (a - b).
                return (s);
            }
        }

        // return (a^n + b^n), computed so as to minimize cancelation errors

        private static double SumOfPowers (double a, double b, int n) {
            if (n % 2 == 0) {
                // for even powers, cancelation is never a problem
                return (MoreMath.Pow(a, n) + MoreMath.Pow(b, n));
            } else {
                // for odd powers, cancelation can be a problem if a and b have opposite signs
                if (Math.Sign(a) != Math.Sign(b)) {
                    // use (a^{2m+1} + b^{2m+1}) = (a + b) (a^{2m} + ... )
                    return (MoreMath.Pow(a, n) + MoreMath.Pow(b, n));
                } else {
                    return (MoreMath.Pow(a, n) + MoreMath.Pow(b, n));
                }
            }
        }

        // returns all powers of a from a^n up to a^{n-1}

        private static double[] Powers (double a, int n) {
            double[] powers = new double[n];
            powers[0] = 1.0;
            for (int m = 1; m < powers.Length; m++) {
                powers[m] = powers[m - 1] * a;
            }
            return (powers);
        }

        /// <inheritdoc />
        public override double CentralMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if ((r % 2) != 0) {
                return (0.0);
            } else {
                return (MoreMath.Pow(range.Width / 2.0, r) / (r + 1));
            }
        }

        /// <inheritdoc />
        public override double Cumulant (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (0.0);
            } else if (r == 1) {
                return (range.Midpoint);
            } else if (r % 2 != 0) {
                // This isn't strictly necessary since BernoulliNumber(r) will return 0, but there is no reason to compute (b-a)^r
                return (0.0);
            } else {
                return (MoreMath.Pow(range.Width, r) * AdvancedIntegerMath.BernoulliNumber(r) / r);
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
        public override double ExcessKurtosis {
            get {
                return (-6.0 / 5.0);
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