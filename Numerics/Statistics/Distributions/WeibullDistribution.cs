using System;
using System.Collections.Generic;
using System.Linq;

using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Weibull distribution.
    /// </summary>
    /// <remarks>
    /// <para>The Weibull distribution is a generalized form of the exponential distribution,
    /// for which the decay probability is not constant, but instead increases or decreases
    /// with time. When the shape parameter is one, the Weibull distribution reduces to the
    /// exponential distribution.</para>
    /// <para>The Weibull distribution is commonly used in engineering applications to
    /// model the time-to-failure of industrial components.</para>
    /// </remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/Weibull_distribution" />
    public sealed class WeibullDistribution : ContinuousDistribution {

        /// <summary>
        /// Initializes a new Weibull distribution.
        /// </summary>
        /// <param name="scale">The scale parameter, which must be positive.</param>
        /// <param name="shape">The shape parameter, which must be positive.</param>
        public WeibullDistribution (double scale, double shape) {
            if (scale <= 0.0) throw new ArgumentOutOfRangeException(nameof(scale));
            if (shape <= 0.0) throw new ArgumentOutOfRangeException(nameof(shape));
            this.scale = scale;
            this.shape = shape;
        }

        private readonly double scale;

        private readonly double shape;

        /// <summary>
        /// Gets the scale parameter of the distribution.
        /// </summary>
        public double Scale {
            get {
                return (scale);
            }
        }

        /// <summary>
        /// Gets the shape parameter of the distribution.
        /// </summary>
        public double Shape {
            get {
                return (shape);
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
            if (x < 0.0) {
                return (0.0);
            } else if (x == 0.0) {
                if (shape < 1.0) {
                    return (Double.PositiveInfinity);
                } else if (shape > 1.0) {
                    return (0.0);
                } else {
                    return (1.0);
                }
            } else {
                double z = x / scale;
                double zz = Math.Pow(z, shape);
                return (Math.Exp(-zz) * zz * shape / x);
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else {
                double z = x / scale;
                double zz = Math.Pow(z, shape);
                return (-MoreMath.ExpMinusOne(-zz));
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x <= 0.0) {
                return (1.0);
            } else {
                double z = x / scale;
                double zz = Math.Pow(z, shape);
                return (Math.Exp(-zz));
            }
        }

        /// <inheritdoc />
        public override double Hazard (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else {
                return (Math.Pow(x / scale, shape) * shape / x);
            }
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            return (scale * Math.Pow(-MoreMath.LogOnePlus(-P), 1.0 / shape));
        }

        /// <inheritdoc />
        public override double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException(nameof(Q));
            return (scale * Math.Pow(-Math.Log(Q), 1.0 / shape)); 
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (scale * AdvancedMath.Gamma(1.0 + 1.0 / shape));
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (scale * Math.Pow(Global.LogTwo, 1.0 / shape));
            }
        }

        // standard deviation inherits from base; nothing special to say about it

        /// <inheritdoc />
        public override double RawMoment (int r) {
            if (r < 0.0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else {
                return (Math.Pow(scale, r) * AdvancedMath.Gamma(1.0 + r / shape));
            }
        }

        /// <inheritdoc />
        public override double CentralMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (0.0);
            } else {
                // for large shape parameters, central moments involve strong cancelations, so integrate to get them
                if (shape < 2.0) {
                    return (CentralMomentFromRawMoment(r));
                } else {
                    return(base.CentralMoment(r));
                }
            }
        }

        // Central moments are difficult to compute because cancelations between raw moments are particularly extreme
        // and I know of no expression that avoids them. For example
        //   C_2 = M_2 - M_1^2 = \lambda * (\Gamma(1 + 2/k) - \Gamma^2(1 + 1/k))
        // For k large, so (1/k) small, both terms are nearly \Gamma(1) ~ 1 and in fact the (1/k) terms cancel too
        // so the leading term is (1/k)^2. Worse, the series is (a) egregiously difficult to calculate and (b)
        // converges poorly unless k is quite large.

        /// <inheritdoc />
        public override double Variance {
            get {
                if (shape < 10.0) {
                    return (MoreMath.Sqr(scale) * (AdvancedMath.Gamma(1.0 + 2.0 / shape) - MoreMath.Sqr(AdvancedMath.Gamma(1.0 + 1.0 / shape))));
                } else {
                    return (base.Variance);
                }
            }
        }

        /// <inheritdoc />
        public override double GetRandomValue (Random rng) {
            return(base.GetRandomValue(rng));
            //if (rng == null) throw new ArgumentNullException("rng");
            //return (scale * Math.Pow(-Math.Log(rng.NextDouble()), 1.0 / shape));
            // Something is wrong with this formula, I need to figure out what it is.
        }
        

        /// <summary>
        /// Computes the Weibull distribution that best fits the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The best fit parameters.</returns>
        /// <remarks>
        /// <para>The returned fit parameters are the <see cref="Shape"/> and <see cref="Scale"/>, in that order.
        /// These are the same parameters, in the same order, that are required by the <see cref="WeibullDistribution(double,double)"/> constructor to
        /// specify a new Weibull distribution.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <exception cref="InvalidOperationException"><paramref name="sample"/> contains non-positive values.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than three values.</exception>
        public static WeibullFitResult FitToSample (Sample sample) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            return (Univariate.FitToWeibull(sample.data));
        }

    }

}