using System;

using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;


namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Gamma distribution.
    /// </summary>
    /// <remarks>
    /// <para>The sum of n exponentially distributed variates is a Gamma distributed variate.</para>
    /// <img src="../images/GammaFromExponential.png" />
    /// <para>When the shape parameter is an integer, the Gamma distribution is also called the Erlang distribution. When
    /// the shape parameter is one, the Gamma distribution reduces to the exponential distribution.</para>
    /// </remarks>
    /// <seealso cref="ExponentialDistribution"/>
    public sealed class GammaDistribution : ContinuousDistribution {

        /// <summary>
        /// Initializes a new instance of a Gamma distribution with the given parameters.
        /// </summary>
        /// <param name="shape">The shape parameter, which must be positive.</param>
        /// <param name="scale">The scale parameter, which must be positive.</param>
        public GammaDistribution (double shape, double scale) {
            if (shape <= 0.0) throw new ArgumentOutOfRangeException(nameof(shape));
            if (scale <= 0.0) throw new ArgumentOutOfRangeException(nameof(scale));
            a = shape;
            s = scale;

            // depending on size of a, store Gamma(a) or -Ln(Gamma(a)) for future use
            if (a < 64.0) {
                ga = AdvancedMath.Gamma(a);
            } else {
                ga = -AdvancedMath.LogGamma(a);
            }

            gammaRng = DeviateGeneratorFactory.GetGammaGenerator(a);

        }

        /// <summary>
        /// Initializes a new instance of the standard Gamma distribution.
        /// </summary>
        /// <param name="shape">The shape parameter, which must be positive.</param>
        public GammaDistribution (double shape) : this(shape, 1.0) {
        }

        // the shape and scale parameters
        private readonly double a, s;

        // Gamma(a), or -Ln(Gamma(a)) for large a
        private readonly double ga;

        // a gamma deviate generator
        // note this generator has the appropriate shape parameter, but is not scaled
        private readonly IDeviateGenerator gammaRng;

        /// <summary>
        /// Gets the shape parameter for the distribution.
        /// </summary>
        public double Shape {
            get {
                return (a);
            }
        }

        /// <summary>
        /// Gets the scale parameter for the distribution.
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
            if (x < 0.0) {
                return (0.0);
            } else {
                double z = x / s;
                // use Gamma(a) or Ln(Gamma(a)) form depending on size of a
                if (ga > 0.0) {
                    return (Math.Pow(z, a - 1.0) * Math.Exp(-z) / ga / s);
                } else {
                    return (Math.Exp((a - 1.0) * Math.Log(z) - z + ga) / s);
                }
            }

        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (a * s);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (a * s * s);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (2.0 / Math.Sqrt(a));
            }
        }

        /// <inheritdoc />
        public override double ExcessKurtosis {
            get {
                return (6.0 / a);
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x < 0.0) {
                return (0.0);
            } else {
                return (AdvancedMath.LeftRegularizedGamma(a, x / s));
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x < 0.0) {
                return (1.0);
            } else {
                return (AdvancedMath.RightRegularizedGamma(a, x / s));
            }
        }

        /// <inheritdoc />
        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else {
                // M_r = s^r Gamma(a+r)/Gamma(c), implying M_{r} = s (a+r) M_{r-1}
                double M = 1.0;
                for (int i = 0; i < r; i++) {
                    M = M * s * (a + i);
                }
                return (M);
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

                // for higher moments, use the recurrence C_{n+1} = n ( C_n + a C_{n-1} )
                // which can be derived from C = U(-n, 1-n-a, -a), where U is the irregular confluent hypergeometric function,
                // and the recursion recursion U(a-1,b-1,z) = (1-b+z) U(a,b,z) + z a U(a+1,b+1,z)

                double C1 = 0.0;
                double C2 = s * s * a;
                for (int i = 2; i < r; i++) {
                    double C3 = s * i * (C2 + a * s * C1);
                    C1 = C2;
                    C2 = C3;
                }
                return (C2);
            }
        }

        /// <inheritdoc />
        public override double Cumulant (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (0.0);
            } else {
                return (a * AdvancedIntegerMath.Factorial(r - 1) * MoreMath.Pow(s, r));
            }
        }

        //  Now comes a whole lot of code just to invert the Gamma CDF, which is a major pain!

        // The probit function is the inverse CDF of the standard normal distribution
        // Since Gamma becomes approximately normal for large shape parameters, this is useful
        // in that regime

        private static double ApproximateProbit (double P) {

            if (P < 0.1) {
                return (-Global.SqrtTwo * AdvancedMath.ApproximateInverseErfc(2.0 * P));
            } else if (P < 0.9) {
                return (Global.SqrtTwo * AdvancedMath.ApproximateInverseErf(2.0 * P - 1.0));
            } else {
                return (Global.SqrtTwo * AdvancedMath.ApproximateInverseErfc(2.0 * (1.0 - P)));
            }

        }
       
        private double ApproximateInverseStandardGamma (double P, double Q) {

            if (P == 0.0) return (0.0);
            if (Q == 0.0) return (Double.PositiveInfinity);

            // compute  (Gamma(a+1) * P)^(1/a)
            double x0;
            if (ga > 0.0) {
                x0 = Math.Pow(a * ga * P, 1.0 / a);
            } else {
                x0 = Math.Exp((Math.Log(a * P) - ga) / a);
            }
            double x1 = x0 / (a + 1.0);
            if (x1 < 0.33) {
                // Do left-tail expansion if possible, because it is quite fast.
                // This is just a term-by-term inversion of the power series
                // P(a,x) = gamma(a,x) / Gamma(a) = e^{-x} x^a / Gamma(a+1) [ 1 + x / (a+1) + ... ]
                double S = 1.0 + x1 + (3.0 * a + 5.0) / (a + 2.0) * x1 * x1 / 2.0 +
                    (8.0 * a * a + 33.0 * a + 31.0) / (a + 2.0) / (a + 3.0) * x1 * x1 * x1 / 3.0;
                return (x0 * S);
            } else {
                if (a > 1.0) {
                    // It is well known that Gamma(a) -> Normal(a,sqrt(a)) for large a,
                    // but it does so exceedingly slowly (skew decreases ~1/sqrt(a) and kurtosis ~1/a).
                    // In the days of yore when people looked for normalizing transforms, Wilson and Hilferty derived that
                    // y = (x/a)^(1/3) ~ Normal(1-1/9a, 1/sqrt(9a)) with much faster decreasing cumulants.
                    // We use this here, having verified that it gives much more accurate values than the naive normal approximation
                    double na = 9.0 * a;
                    double y = 1.0 - 1.0 / na + ApproximateProbit(P) / Math.Sqrt(na);
                    return (a * MoreMath.Pow(y, 3));
                } else {
                    // For small a, use a very crude right-tail approximation
                    // this will fail if Gamma(a) * Q > 1, but given our range for the left-tail
                    // approximation above, we have ~0.5 < Gamma(a) * Q < ~0.8 for 0 < a < 1.
                    double log;
                    if (ga > 0.0) {
                        log = - Math.Log(ga * Q);
                    } else {
                        log = ga - Math.Log(Q);
                    }
                    return (log - (1.0 - a) * Math.Log(log));
                }
            }

        }

        private double InverseLeftStandardGamma (double P) {

            // start with an approximation
            double x = ApproximateInverseStandardGamma(P, 1.0 - P);

            // zero x or infinite x will cause problems with the iterative method below
            // near zero (and at infinity) the approximation becomes good to double precision anyway
            if ((x < Global.Accuracy) || Double.IsPositiveInfinity(x)) return (x);

            // refine with Halley's method, i.e. second order Newton method
            double y = AdvancedMath.LeftRegularizedGamma(a, x) - P;
            for (int i = 0; i < 16; i++) {
                
                double r;
                if (ga > 0.0) {
                    r = y * ga * Math.Exp(x) / Math.Pow(x, a);
                } else {
                    r = y * Math.Exp(x - ga - a * Math.Log(x));
                }
                double dx = -r * x / (1.0 - r * (a - 1.0 - x) / 2.0);
                x += dx;

                if (Math.Abs(dx) <= Global.Accuracy * Math.Abs(x)) return (x);

                double y_new = AdvancedMath.LeftRegularizedGamma(a, x) - P;
                if (Math.Abs(y_new) >= Math.Abs(y)) return (x);
                y = y_new;

            }
            
            throw new NonconvergenceException();
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            return (s * InverseLeftStandardGamma(P));
        }

        /// <inheritdoc />
        public override double GetRandomValue (Random rng) {
            if (rng == null) throw new ArgumentNullException(nameof(rng));
            return (s * gammaRng.GetNext(rng));
        }

        /// <summary>
        /// Computes the Gamma distribution that best fits the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The best fit parameters.</returns>
        /// <remarks>
        /// <para>The returned fit parameters are the <see cref="Shape"/> and <see cref="Scale"/>, in that order.
        /// These are the same parameters, in the same order, that are required by the <see cref="GammaDistribution(double,double)"/> constructor to
        /// specify a new Gamma distribution.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <exception cref="InvalidOperationException"><paramref name="sample"/> contains non-positive values.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than three values.</exception>
        public static GammaFitResult FitToSample (Sample sample) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            return (Univariate.FitToGamma(sample.data));
        }

#if FUTURE
        private void PsiDeficitAndDerivative (double x, out double f, out double fp) {

            // define the psi deficit function f = log(x) - psi(x)
            // for large x, psi(x) ~ log(x), so f gets small
            // f > 0 for all x > 0

            f = 0.0;
            fp = 0.0;

            // for small x, use f(x) = f(x+1) - 1/x, and therefore f'(x) = f'(x+1) + 1/x^2, to advance x
            while (x < 16.0) {
                double xi = 1.0 / x;
                f -= xi;
                fp += xi * xi;
                x += 1.0;
            }

            // once x is large enough, use asymptotic expansion
            //   f = \frac{1}{2x} + \sum_{n=1}^{\infty} \frac{B_{2n}}{2n x^{2n}}
            // and therefore
            //   f' = \frac{-1}{2x^2} \left[ \frac{1}{2x} + \sum_{n=1}^{\infty} \frac{B_{2n}}{x^{2n}}
            // to compute f and f'

            f += 1.0 / (2.0 * x);
            fp -= 1.0 / (2.0 * x * x);

            double xx = x2;
            for (int n = 1; n < AdvancedIntegerMath.Bernoulli.Length; n++) {
                double f_old = f; double fp_old = fp;
                f += AdvancedIntegerMath.Bernoulli[n] / (2 * n) * xx;
                fp -= AdvancedIntegerMath.Bernoulli[n] * xx / x;
                if (f == f_old && fp == fp_old) {
                    return;
                }
                xx *= x2;
            }

            throw new NonconvergenceException();
        }
#endif
    }
}
