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
        public static RayleighFitResult FitToSample (IReadOnlyList<double> sample) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 2) throw new InsufficientDataException();

            // It's easy to follow maximum likelihood prescription, because there is only one
            // parameter and the functional form is fairly simple.

            // \ln L = \sum_i p_i = = \sum_i \left[ \ln x_i - 2 \ln \sigma - \frac{1}{2} \frac{x_i^2}{\sigma^2 \right]
            // \frac{\partial \ln L}{\partial \sigma} = \sum_i \left[ \frac{x_i^2}{\sigma^3} - \frac{2}{\sigma} \right]
            // \frac{\partial^2 \ln L}{\partial \sigma^2} = \sum_i \left[ \frac{2}{\sigma^2} - \frac{3 x_i^2}{\sigma^4} \right]

            // Set the first derivative to zero to obtain
            //   \hat{\sigma}^2 = \frac{1}{2n} \sum_i x_i^2
            //   \hat{\sigma} = \sqrt{\frac{1}{2} \left< x^2 \right>}
            // and plug this value into the second derivative to obtain
            //   \frac{\partial^2 \ln L}{\partial \sigma^2} = - \frac{4n}{\sigma^2}
            // at the minimum, so
            //   \delta \sigma = \frac{\sigma}{\sqrt{4n}}

            // Next consider bias. We know from the moments of the distribution that < x^2 > = 2 \sigma^2, so \hat{\sigma}^2 is
            // unbiased. But its square root \hat{\sigma} is not. To get exact distribution of \hat{\sigma},
            //    x_i \sim Rayleigh(\sigma)
            //    ( \frac{x_i}{\sigma} )^2 \sim \chi^2(2)
            //    \sum_i ( \frac{x_i}{\sigma} )^2 \sim \chi^2(2n)
            //    \left[ \sum_i ( \frac{x_i}{\sigma} )^2 \right]^{1/2} \sim \chi(2n)
            // And if z \sim \chi(k) then
            //    E(z) = \sqrt{2} \frac{\Gamma((k + 1)/2)}{\Gamma(k / 2)}
            //    V(z) = k - [ E(z) ]^2
            // Here k = 2n and our estimator \hat{\sigma} = z / sqrt{2n}, so
            //    E(\hat{\sigma}) = \frac{\Gamma(n + 1/2)}{\sqrt{n} \Gamma(n)} \sigma = \frac{(n)_{1/2}}{\sqrt{n}} \sigma
            //    V(\hat{\sigma}) = \left[ 1 - \frac{(n)_{1/2}^2}{n} \right] \sigma^2
            // We can use series expansion to verify that
            //    E(\hat{\sigma}) = \left[ 1 - \frac{1}{8n} + \cdots \right] \sigma
            //    V(\hat{\sigma}) = \left[ \frac{1}{4n} - \frac{1}{32n^2} + \cdots \right] \sigma^2
            // our estimator is asymptotically unbiased and its asymptotic variance agrees with our assessment.

            // We correct for the bias of \hat{\sigma} by multiplying by \frac{\sqrt{n}}{(n)_{1/2}}. We could
            // correct our variance estimate too, but in order to evaluate the correction at high n,
            // we need to do a series expansion, and I regard it as less important that the error estimate
            // be exact.

            int n = sample.Count;
            double s = Math.Sqrt(sample.RawMoment(2) / 2.0) * (Math.Sqrt(n) / AdvancedMath.Pochhammer(n, 1.0 / 2.0));
            double ds = s / Math.Sqrt(4.0 * n);

            RayleighDistribution distribution = new RayleighDistribution(s);
            TestResult test = sample.KolmogorovSmirnovTest(distribution);

            return (new RayleighFitResult(new UncertainValue(s, ds), distribution, test));
        }
    }
}
