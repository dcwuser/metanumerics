using System;

using Meta.Numerics.Functions;

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
    public class GammaDistribution : Distribution {

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

        }

        /// <summary>
        /// Initializes a new instance of the standard Gamma distribution.
        /// </summary>
        /// <param name="shape">The shape parameter, which must be positive.</param>
        public GammaDistribution (double shape) : this(shape, 1.0) {
        }

        double a, s;
        double ga;

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

        // the remainder of this code is all this is just to invert, which is a major pain!

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

            double x = ApproximateInverseStandardGamma(P, 1.0 - P);
            // zero x or infinite x will cause problems with the iterative method below
            // near zero (and at infinity) the approximation becomes good to double precision anyway
            if ((x < Global.Accuracy) || Double.IsPositiveInfinity(x)) return (x);

            double y = AdvancedMath.LeftRegularizedGamma(a, x) - P;
            for (int i = 0; i < 16; i++) {

                //Console.WriteLine("  x={0:g16} y={1}", x, y);
                
                double r;
                if (ga > 0.0) {
                    r = y * ga * Math.Exp(x) / Math.Pow(x, a);
                } else {
                    r = y * Math.Exp(x - ga - a * Math.Log(x));
                }
                double dx = -r * x / (1.0 - r * (a - 1.0 - x) / 2.0);
                x += dx;
                //Console.WriteLine("  dx={0} x={1}", dx, x);

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

    }
}
