using System;
using System.Collections.Generic;

using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

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
            return (scale * Math.Pow(-MoreMath.LogOnePlus(-P), 1.0 / shape));
        }

        /// <inheritdoc />
        public override double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException("Q");
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
        public override double Moment (int r) {
            if (r < 0.0) {
                throw new ArgumentOutOfRangeException("r");
            } else if (r == 0) {
                return (1.0);
            } else {
                return (Math.Pow(scale, r) * AdvancedMath.Gamma(1.0 + r / shape));
            }
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException("r");
            } else if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (0.0);
            } else {
                // for large shape parameters, central moments involve strong calculations, so integrate to get them
                if (shape < 2.0) {
                    return (CentralMomentFromRawMoment(r));
                } else {
                    return(base.MomentAboutMean(r));
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
                    return (MoreMath.Pow2(scale) * (AdvancedMath.Gamma(1.0 + 2.0 / shape) - MoreMath.Pow2(AdvancedMath.Gamma(1.0 + 1.0 / shape))));
                } else {
                    return (base.Variance);
                }
            }
        }

        /// <summary>
        /// Computes the Weibull distribution that best fits the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The best fit parameters.</returns>
        /// <remarks>
        /// <para>The returned fit parameters are the <see cref="ShapeParameter"/> and <see cref="ScaleParameter"/>, in that order.
        /// These are the same parameters, in the same order, that are required by the <see cref="WeibullDistribution(double,double)"/> constructor to
        /// specify a new Weibull distribution.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <exception cref="InvalidOperationException"><paramref name="sample"/> contains non-positive values.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than three values.</exception>
        public static FitResult FitToSample (Sample sample) {

            if (sample == null) throw new ArgumentNullException("sample");
            if (sample.Count < 3) throw new InsufficientDataException();
            if (sample.Minimum <= 0.0) throw new InvalidOperationException();

            // The log likelyhood function is
            //   \log L = N \log k + (k-1) \sum_i \log x_i - N K \log \lambda - \sum_i \left(\frac{x_i}{\lambda}\right)^k
            // Taking derivatives, we get
            //   \frac{\partial \log L}{\partial \lambda} = - \frac{N k}{\lambda} + \sum_i \frac{k x_i^k}{\lambda^{k+1}}
            //   \frac{\partial \log L}{\partial k} =
            // Setting the first expression to zero and solving for \lambda gives
            //   \lambda = \left( N^{-1} \sum_i x_i^k \right)^{1/k} = ( < x^k > )^{1/k}
            // which allows us to reduce the problem from 2D to 1D.
            // By the way, using the expression for the moment < x^k > of the Weibull distribution, you can show there is
            // no bias to this result even for finite samples.
            // Setting the second expression to zero gives
            //   \frac{1}{k} = \frac{1}{N} \sum_i \left[ \left( \frac{x_i}{\lambda} \right)^k - 1 \right] \log \left(\frac{x_i}{\lambda}\right)
            // which, given the equation for \lambda as a function of k derived from the first expression, is an implicit equation for k.
            // It cannot be solved in closed form, but we have now reduced our problem to finding a root in one-dimension.

            // We need a starting guess for k.
            // The method of moments equations are not solvable for the parameters in closed form
            // but the scale parameter drops out of the ratio of the 1/3 and 2/3 quantile points
            // and the result is easily solved for the shape parameter
            //   k = \frac{\log 2}{\log\left(\frac{x_{2/3}}{x_{1/3}}\right)}
            double x1 = sample.InverseLeftProbability(1.0 / 3.0);
            double x2 = sample.InverseLeftProbability(2.0 / 3.0);
            double k0 = Global.LogTwo / Math.Log(x2 / x1);
            // Given the shape paramter, we could invert the expression for the mean to get
            // the scale parameter, but since we have an expression for \lambda from k, we
            // dont' need it.
            //double s0 = sample.Mean / AdvancedMath.Gamma(1.0 + 1.0 / k0);

            // Construct the \lambda function
            Func<double,double> lambdaFunction = delegate (double k) {
                double s = 0.0;
                foreach (double x in sample) {
                    s += Math.Pow(x, k);
                }
                return (Math.Pow(s / sample.Count, 1.0 / k));
            };

            // Construct the k function whoose root we want
            Func<double, double> kFunction = delegate (double k) {
                double lambda = lambdaFunction(k);
                double s = 0.0;
                foreach (double x in sample) {
                    double r = x / lambda;
                    s += (Math.Pow(r, k) - 1.0) * Math.Log(r);
                }
                return (s / sample.Count - 1.0 / k);
            };

            // Find its root, starting from our guess
            double k1 = FunctionMath.FindZero(kFunction, k0);

            // Get the lambda value corresponing to our best k value
            double lambda1 = lambdaFunction(k1);

            // We need the curvature matrix at the minimum of our log likelyhood function
            // to determine the covariance matrix. Taking more derivatives...
            //    \frac{\partial^2 \log L} = \frac{N k}{\lambda^2} - \sum_i \frac{k(k+1) x_i^k}{\lambda^{k+2}}
            //    = - \frac{N k^2}{\lambda^2}
            // The second expression follows by inserting the first-derivative-equal-zero relation into the first.
            // For k=1, this agrees with the variance formula for the mean of the best-fit exponential.

            // Derivatives involving k are less simple.

            // So we will need the means < (x/lambda) log(x/lambda) > and < (x/lambda) log^2(x/lambda) >
            double mpl = 0.0; double mpl2 = 0.0;
            foreach (double x in sample) {
                double r = x / lambda1;
                double p = Math.Pow(r, k1);
                double l = Math.Log(r);
                double pl = p * l;
                double pl2 = pl * l;
                mpl += pl;
                mpl2 += pl2;
            }
            mpl = mpl / sample.Count;
            mpl2 = mpl2 / sample.Count;

            // Construct the curvature matrix and invert it.
            SymmetricMatrix B = new SymmetricMatrix(2);
            B[0, 0] = sample.Count * MoreMath.Pow2(k1 / lambda1);
            B[0, 1] = -sample.Count * k1 / lambda1 * mpl;
            B[1, 1] = sample.Count * (1.0 / MoreMath.Pow2(k1) + mpl2);
            SymmetricMatrix C = B.CholeskyDecomposition().Inverse();

            // Do a KS test to compare sample to best-fit distribution
            Distribution distribution = new WeibullDistribution(lambda1, k1);
            TestResult test = sample.KolmogorovSmirnovTest(distribution);

            // return the result
            return (new FitResult(new double[] {lambda1, k1}, C, test));

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