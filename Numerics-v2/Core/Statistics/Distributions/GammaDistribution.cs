using System;

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
    public sealed class GammaDistribution : Distribution {

        /// <summary>
        /// Initializes a new instance of a Gamma distribution with the given parameters.
        /// </summary>
        /// <param name="shape">The shape parameter, which must be positive.</param>
        /// <param name="scale">The scale parameter, which must be positive.</param>
        public GammaDistribution (double shape, double scale) {
            if (shape <= 0.0) throw new ArgumentOutOfRangeException("shape");
            if (scale <= 0.0) throw new ArgumentOutOfRangeException("scale");
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
        public double ShapeParameter {
            get {
                return (a);
            }
        }

        /// <summary>
        /// Gets the scale parameter for the distribution.
        /// </summary>
        public double ScaleParameter {
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
        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else {
                // M_n = s^n Gamma(a+n)/Gamma(c), implying M_{n} = s (a+n) M_{n-1}
                double M = 1.0;
                for (int i = 0; i < n; i++) {
                    M = M * s * (a + i);
                }
                return (M);
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

                // for higher moments, use the recurrence C_{n+1} = n ( C_n + a C_{n-1} )
                // which can be derived from C = U(-n, 1-n-a, -a), where U is the irregular confluent hypergeometric function,
                // and the recurrsion recursion U(a-1,b-1,z) = (1-b+z) U(a,b,z) + z a U(a+1,b+1,z)

                double C1 = 0.0;
                double C2 = s * s * a;
                for (int i = 2; i < n; i++) {
                    double C3 = s * i * (C2 + a * s * C1);
                    C1 = C2;
                    C2 = C3;
                }
                return (C2);
            }
        }

        internal override double Cumulant (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (0.0);
            } else {
                return (a * AdvancedIntegerMath.Factorial(n - 1));
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
                // do left-tail expansion if possible, because it is quite fast
                // this is just a term-by-term inversion of the power series
                // P(a,x) = gamma(a,x) / Gamma(a) = e^{-x} x^a / Gamma(a+1) [ 1 + x / (a+1) + ... ]
                double S = 1.0 + x1 + (3.0 * a + 5.0) / (a + 2.0) * x1 * x1 / 2.0 +
                    (8.0 * a * a + 33.0 * a + 31.0) / (a + 2.0) / (a + 3.0) * x1 * x1 * x1 / 3.0;
                return (x0 * S);
            } else {
                if (a > 1.0) {
                    // it is well known that Gamma(a) -> Normal(a,sqrt(a)) for large a
                    // but it does so exceedingly slowly (skew decreases ~1/sqrt(a) and kurtosis ~1/a)
                    // in the days of yore when people looked for normalizing transforms, Wilson and Hilferty derived that
                    // y = (x/a)^(1/3) ~ Normal(1-1/9a, 1/sqrt(9a)) with much faster decreasing cumulants
                    // we use this here, having verified that it gives much more accurate values than the naive normal approximation
                    double na = 9.0 * a;
                    double y = 1.0 - 1.0 / na + ApproximateProbit(P) / Math.Sqrt(na);
                    return (a * MoreMath.Pow(y, 3));
                } else {
                    // for small a, use a very crude right-tail approximation
                    // this will fail if Gamma(a) * Q > 1, but given our range for the left-tail
                    // approximation above, we have ~0.5 < Gamma(a) * Q < ~0.8 for 0 < a < 1
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
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            return (s * InverseLeftStandardGamma(P));
        }

        /// <inheritdoc />
        public override double GetRandomValue (Random rng) {
            if (rng == null) throw new ArgumentNullException("rng");
            return (s * gammaRng.GetNext(rng));
        }

        // routines for maximum likelyhood fitting

        /// <summary>
        /// Computes the Gamma distribution that best fits the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The best fit parameters.</returns>
        /// <remarks>
        /// <para>The returned fit parameters are the <see cref="ShapeParameter"/> and <see cref="ScaleParameter"/>, in that order.
        /// These are the same parameters, in the same order, that are required by the <see cref="GammaDistribution(double,double)"/> constructor to
        /// specify a new Gamma distribution.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <exception cref="InvalidOperationException"><paramref name="sample"/> contains non-positive values.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than three values.</exception>
        public static FitResult FitToSample (Sample sample) {

            if (sample == null) throw new ArgumentNullException("sample");
            if (sample.Count < 3) throw new InsufficientDataException();

            // The log likelyhood of a sample given k and s is
            //   \log L = (k-1) \sum_i \log x_i - \frac{1}{s} \sum_i x_i - N \log \Gamma(k) - N k \log s
            // Differentiating,
            //   \frac{\partial \log L}{\partial s} = \frac{1}{s^2} \sum_i x_i - \frac{Nk}{s}
            //   \frac{\partial \log L}{\partial k} = \sum_i \log x_i - N \psi(k) - N \log s
            // Setting the first equal to zero gives
            //   k s = N^{-1} \sum_i x_i = <x>
            //   \psi(k) + \log s = N^{-1} \sum_i \log x_i = <log x>
            // Inserting the first into the second gives a single equation for k
            //   \log k - \psi(k) = \log <x> - <\log x>
            // Note the RHS need only be computed once.
            // \log k > \psi(k) for all k, so the RHS had better be positive. They get
            // closer for large k, so smaller RHS will produce a larger k.

            double s = 0.0;
            foreach (double x in sample) {
                if (x <= 0.0) throw new InvalidOperationException();
                s += Math.Log(x);
            }
            s = Math.Log(sample.Mean) - s / sample.Count;

            // We can get an initial guess for k from the method of moments
            //   \frac{\mu^2}{\sigma^2} = k

            double k0 = MoreMath.Pow2(sample.Mean) / sample.Variance;

            // Since 1/(2k) < \log(k) - \psi(k) < 1/k, we could get a bound; that
            // might be better to avoid the solver running into k < 0 territory

            double k1 = FunctionMath.FindZero(k => (Math.Log(k) - AdvancedMath.Psi(k) - s), k0);

            double s1 = sample.Mean / k1;

            // Curvature of the log likelyhood is straightforward
            //   \frac{\partial^2 \log L}{\partial s^2} = -\frac{2}{s^3} \sum_i x_i + \frac{Nk}{s^2} = - \frac{Nk}{s^2}
            //   \frac{\partial^2 \log L}{\partial k \partial s} = - \frac{N}{s}
            //   \frac{\partial^2 \log L}{\partial k^2} = - N \psi'(k)
            // This gives the curvature matrix and thus via inversion the covariance matrix.

            SymmetricMatrix B = new SymmetricMatrix(2);
            B[0, 0] = sample.Count * AdvancedMath.Psi(1, k1);
            B[0, 1] = sample.Count / s1;
            B[1, 1] = sample.Count * k1 / MoreMath.Pow2(s1);
            SymmetricMatrix C = B.CholeskyDecomposition().Inverse();

            // Do a KS test for goodness-of-fit
            TestResult test = sample.KolmogorovSmirnovTest(new GammaDistribution(k1, s1));

            return (new FitResult(new double[] { k1, s1 }, C, test));
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
