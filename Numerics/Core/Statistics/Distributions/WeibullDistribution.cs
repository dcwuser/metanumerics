using System;
using System.Collections.Generic;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Weibull distribution.
    /// </summary>
    /// <remarks>
    /// <para>The Weibull distribution is a generalized form of the exponential distriubtion
    /// for which the decay probability is not constant, but instead increases with time
    /// (for shape parameters greater than one) or, less commonly, decreases with time (for
    /// shape parameters less than one). When the shape parameter is one, the Weibull
    /// distribution reduces to the exponential distribution.</para>
    /// <para>The Weibull distribution is commonly used in engineering applications to
    /// model the time-to-failure of industrial componets.</para>
    /// </remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/Weibull_distribution" />
    public sealed class WeibullDistribution : Distribution, IParameterizedDistribution {

        /// <summary>
        /// Initializes a new Weibull distribution.
        /// </summary>
        /// <param name="scale">The scale parameter, which must be positive.</param>
        /// <param name="shape">The shape parameter, which must be positive.</param>
        public WeibullDistribution (double scale, double shape) {
            if (scale <= 0.0) throw new ArgumentOutOfRangeException("scale");
            if (shape <= 0.0) throw new ArgumentOutOfRangeException("shape");
            this.scale = scale;
            this.shape = shape;
        }

        private double scale;

        private double shape;

        /// <summary>
        /// Gets the scale parameter of the distribution.
        /// </summary>
        public double ScaleParameter {
            get {
                return (scale);
            }
        }

        /// <summary>
        /// Gets the shape parameter of the distribution.
        /// </summary>
        public double ShapeParameter {
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
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            return (scale * Math.Pow(-Math.Log(1.0 - P), 1.0 / shape));
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
        public override double Moment (int n) {
            if (n < 0.0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else {
                return (Math.Pow(scale, n) * AdvancedMath.Gamma(1.0 + n / shape));
            }
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (0.0);
            } else {
                // for large shape parameters, central moments involve strong calculations, so integrate to get them
                if (shape < 2.0) {
                    return (CentralMomentFromRawMoment(n));
                } else {
                    return(base.MomentAboutMean(n));
                }
            }
        }

        /// <summary>
        /// Computes the Weibull distribution that best fits the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The result of the fit.</returns>
        /// <exception cref="ArgumentNullException">The <paramref name="sample"/> is null.</exception>
        /// <exception cref="InsufficientDataException">The <paramref name="sample"/> contains less than three data points.</exception>
        /// <exception cref="InvalidOperationException">The <paramref name="sample"/> contains non-positive values.</exception>
        public static FitResult FitToSample (Sample sample) {

            if (sample == null) throw new ArgumentNullException("sample");
            if (sample.Count < 3) throw new InsufficientDataException();
            if (sample.Minimum <= 0.0) throw new InvalidOperationException();

            // guess parameters, using formulas for 1/3 and 2/3 quantile points and mean
            double x1 = sample.InverseLeftProbability(1.0 / 3.0);
            double x2 = sample.InverseLeftProbability(2.0 / 3.0);
            double k0 = 0.996768 / Math.Log(x2 / x1);
            double s0 = sample.Mean / AdvancedMath.Gamma(1.0 + 1.0 / k0);

            // construct the log likeyhood function
            Func<double[], double> f = delegate(double[] p) {
                double lnL = 0.0;
                foreach (double x in sample) {
                    double r = x / p[0];
                    lnL -= (p[1] - 1.0) * Math.Log(r) - Math.Pow(r, p[1]);
                }
                lnL -= sample.Count * Math.Log(p[1] / p[0]);
                return (lnL);
            };

            // mamimize it
            SpaceExtremum fm = FunctionMath.FindMinimum(f, new double[] { s0, k0 });

            double[] v = fm.Location();
            Distribution distribution = new WeibullDistribution(v[0], v[1]);
            TestResult test = sample.KolmogorovSmirnovTest(distribution);

            // return the result
            return (new FitResult(fm.Location(), fm.Curvature().CholeskyDecomposition().Inverse(), test));

        }


        double[] IParameterizedDistribution.GetParameters () {
            return (new double[] { scale, shape });
        }

        void IParameterizedDistribution.SetParameters (IList<double> parameters) {
            if (parameters == null) throw new ArgumentNullException("parameters");
            if (parameters.Count != 2) throw new DimensionMismatchException();
            if (parameters[0] <= 0.0) throw new ArgumentOutOfRangeException("parameters");
            if (parameters[1] <= 0.0) throw new ArgumentOutOfRangeException("parameters");
            scale = parameters[0];
            shape = parameters[1];
        }

        double IParameterizedDistribution.Likelihood (double x) {
            return (ProbabilityDensity(x));
        }

    }

}