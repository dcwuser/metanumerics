using System;
using System.Collections.Generic;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Rayleigh distribution.
    /// </summary>
    /// <remarks>
    /// <para>A Rayleigh distribution is the distribution of the magnitude of a two-dimensional vector whose
    /// components are normally distributed.</para>
    /// <para>A standard Rayleigh distribution is equivalent to a <see cref="ChiDistribution"/> with &#x3BD; = 2.</para>
    /// </remarks>
    /// <seealso href="https://en.wikipedia.org/wiki/Rayleigh_distribution"/>
    /// <seealso href="https://mathworld.wolfram.com/RayleighDistribution.html"/>
    public sealed class RayleighDistribution : ContinuousDistribution {

        /// <summary>
        /// Initializes a new Rayleigh distribution with the given scale parameter.
        /// </summary>
        /// <param name="scale">The scale parameter, which must be positive.</param>
        public RayleighDistribution (double scale) {
            if (scale <= 0.0) throw new ArgumentOutOfRangeException(nameof(scale));
            this.s = scale;
        }

        private readonly double s;

        /// <summary>
        /// Gets the scale parameter of the distribution.
        /// </summary>
        public double Scale {
            get {
                return (s);
            }
        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else {
                double z = x / s;
                return (z * Math.Exp(-0.5 * z * z) / s);
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else {
                double z = x / s;
                return (-MoreMath.ExpMinusOne(-0.5 * z * z));
            }
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            return (Math.Sqrt(-2.0 * MoreMath.LogOnePlus(-P)) * s);
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x <= 0.0) {
                return (1.0);
            } else {
                double z = x / s;
                return (Math.Exp(-0.5 * z * z));
            }
        }

        /// <inheritdoc />
        public override double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException(nameof(Q));
            return (Math.Sqrt(-2.0 * Math.Log(Q)) * s);
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (Math.Sqrt(2.0 * Global.LogTwo) * s);
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (Global.SqrtHalfPI * s);
            }
        }

        /// <inheritdoc />
        public override double StandardDeviation {
            get {
                return (Math.Sqrt((4.0 - Math.PI) / 2.0) * s);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return ((4.0 - Math.PI) / 2.0 * s * s);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (2.0 * Global.SqrtPI * (Math.PI - 3.0) / Math.Pow(4.0 - Math.PI, 3.0 / 2.0));
            }
        }

        /// <inheritdoc />
        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else {
                return (Math.Pow(2.0, r / 2.0) * AdvancedMath.Gamma(1.0 + r / 2.0) * MoreMath.Pow(s, r));
            }
        }

        /// <inheritdoc />
        public override double Hazard (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else {
                return (x / (s * s));
            }
        }

        /// <summary>
        /// Fits a Rayleigh distribution to a sample.
        /// </summary>
        /// <param name="sample">The sample to fit, which must have at least 2 values.</param>
        /// <returns>The fit result. The only parameter is the scale parameter.</returns>
        public static RayleighFitResult FitToSample (Sample sample) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            return (Univariate.FitToRayleigh(sample.data));
        }
    }
}
