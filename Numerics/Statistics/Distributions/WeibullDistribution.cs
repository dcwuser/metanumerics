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
    /// <para>The Weibull distribution is a generalized form of the exponential distriubtion
    /// for which the decay probability is not constant, but instead increases with time
    /// (for shape parameters greater than one) or, less commonly, decreases with time (for
    /// shape parameters less than one). When the shape parameter is one, the Weibull
    /// distribution reduces to the exponential distribution.</para>
    /// <para>The Weibull distribution is commonly used in engineering applications to
    /// model the time-to-failure of industrial componets.</para>
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
                return (Math.Pow(x / scale, shape - 1.0));
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
                    return (MoreMath.Pow2(scale) * (AdvancedMath.Gamma(1.0 + 2.0 / shape) - MoreMath.Pow2(AdvancedMath.Gamma(1.0 + 1.0 / shape))));
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
        public static WeibullFitResult FitToSample (IReadOnlyList<double> sample) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 3) throw new InsufficientDataException();
            foreach (double value in sample) {
                if (value <= 0.0) throw new InvalidOperationException();
            }

            // The log likelyhood function is
            //   \log L = N \log k + (k-1) \sum_i \log x_i - N K \log \lambda - \sum_i \left(\frac{x_i}{\lambda}\right)^k
            // Taking derivatives, we get
            //   \frac{\partial \log L}{\partial \lambda} = - \frac{N k}{\lambda} + \sum_i \frac{k}{\lambda} \left(\frac{x_i}{\lambda}\right)^k
            //   \frac{\partial \log L}{\partial k} =\frac{N}{k} + \sum_i \left[ 1 - \left(\frac{x_i}{\lambda}\right)^k \right] \log \left(\frac{x_i}{\lambda}\right)
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

            // Simply handing our 1D function to a root-finder works fine until we start to encounter large k. For large k,
            // even just computing \lambda goes wrong because we are taking x_i^k which overflows. Horst Rinne, "The Weibull
            // Distribution: A Handbook" describes a way out. Basically, we first move to variables z_i = \log(x_i) and
            // then w_i = z_i - \bar{z}. Then lots of factors of e^{k \bar{z}} cancel out and, even though we still do
            // have some e^{k w_i}, the w_i are small and centered around 0 instead of large and centered around \lambda.

            //Sample transformedSample = sample.Copy();
            //transformedSample.Transform(x => Math.Log(x));
            double[] transformedSample = new double[sample.Count];
            for (int j = 0; j < sample.Count; j++) transformedSample[j] = Math.Log(sample[j]);
            double zbar = transformedSample.Mean();
            for (int j = 0; j < transformedSample.Length; j++) transformedSample[j] -= zbar;

            // After this change of variable the 1D function to zero becomes
            //   g(k) = \sum_i ( 1 - k w_i ) e^{k w_i}
            // It's easy to show that g(0) = n and g(\infinity) = -\infinity, so it must cross zero. It's also easy to take
            // a derivative
            //   g'(k) = - k \sum_i w_i^2 e^{k w_i}
            // so we can apply Newton's method.

            int i = 0;
            double k1 = k0;
            while (true) {
                i++;
                double g = 0.0;
                double gp = 0.0;
                foreach (double w in transformedSample) {
                    double e = Math.Exp(k1 * w);
                    g += (1.0 - k1 * w) * e;
                    gp -= k1 * w * w * e;
                }
                double dk = -g / gp;
                k1 += dk;
                if (Math.Abs(dk) <= Global.Accuracy * Math.Abs(k1)) break;
                if (i >= Global.SeriesMax) throw new NonconvergenceException();
            }
            
            // The corresponding lambda can also be expressed in terms of zbar and w's.

            double t = 0.0;
            foreach (double w in transformedSample) {
                t += Math.Exp(k1 * w);
            }
            t /= transformedSample.Length;
            double lambda1 = Math.Exp(zbar) * Math.Pow(t, 1.0 / k1);
            
            // We need the curvature matrix at the minimum of our log likelyhood function
            // to determine the covariance matrix. Taking more derivatives...
            //    \frac{\partial^2 \log L} = \frac{N k}{\lambda^2} - \sum_i \frac{k(k+1) x_i^k}{\lambda^{k+2}}
            //    = - \frac{N k^2}{\lambda^2}
            // The second expression follows by inserting the first-derivative-equal-zero relation into the first.
            // For k=1, this agrees with the variance formula for the mean of the best-fit exponential.

            // Derivatives involving k are less simple.

            // We end up needing the means < (x/lambda)^k log(x/lambda) > and < (x/lambda)^k log^2(x/lambda) >

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

            // See if we can't do any better here. Transforming to zbar and w's looked ugly, but perhaps it
            // can be simplified? One interesting observation: if we take expectation values (which gives
            // the Fisher information matrix) the entries become simple:
            //   B_{\lambda \lambda} = \frac{N k^2}{\lambda^2}
            //   B_{\lambda k} = -\Gamma'(2) \frac{N}{\lambda}
            //   B_{k k } = [1 + \Gamma''(2)] \frac{N}{k^2}
            // Would it be bad to just use these directly?
            
            // Construct the curvature matrix and invert it.
            SymmetricMatrix B = new SymmetricMatrix(2);
            B[0, 0] = sample.Count * MoreMath.Sqr(k1 / lambda1);
            B[0, 1] = -sample.Count * k1 / lambda1 * mpl;
            B[1, 1] = sample.Count * (1.0 / MoreMath.Pow2(k1) + mpl2);
            SymmetricMatrix C = B.CholeskyDecomposition().Inverse();

            // Do a KS test to compare sample to best-fit distribution
            ContinuousDistribution distribution = new WeibullDistribution(lambda1, k1);
            TestResult test = sample.KolmogorovSmirnovTest(distribution);

            // return the result
            return (new WeibullFitResult(lambda1, k1, C, test));
            //return (new FitResult(new double[] {lambda1, k1}, C, test));

        }

    }

}