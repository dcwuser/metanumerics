using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Gamma distribution.
    /// </summary>
    /// <remarks>
    /// <para>The sum of n exponentially distributed variates is a Gamma distributed variate.</para>
    /// <img src="../images/GammaFromExponential.png" />
    /// <para>When the shape parameter is an integer, the Gamma distribution is also called the Erlang distribution.</para>
    /// </remarks>
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

        public double ShapeParameter {
            get {
                return (a);
            }
        }

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
            double z = x / s;
            // use Gamma(a) or Ln(Gamma(a)) form depending on size of a
            if (ga > 0.0) {
                return (Math.Pow(z, a - 1.0) * Math.Exp(-z) / ga / s);
            } else {
                return (Math.Exp((a - 1.0) * Math.Log(z) - z + ga) / s);
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
            if (x < 0) {
                return (0.0);
            } else {
                return (AdvancedMath.LeftRegularizedGamma(a, x / s));
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x < 0) {
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

        // the remainder of this code is all this is just to invert, which is a major pain!

        private static double ApproximateProbit (double P) {

            if (P < 0.1) {
                return (-ApproximateInverseTailProbability(2.0 * P));
            } else if (P < 0.9) {
                return (ApproximateInverseCentralProbability(2.0 * P - 1.0));
            } else {
                return (ApproximateInverseTailProbability(2.0 * (1.0 - P)));
            }

        }

        private static double ApproximateInverseCentralProbability (double P) {

            // this is just a term-by-term inversion of the series for erf(z)

            double z = Math.Sqrt(Math.PI / 2.0) * P;
            double zz = z * z;
            double S = 1.0 + zz / 6.0 + 7.0 / 120.0 * zz * zz + 127.0 / 5040.0 * zz * zz * zz;
            return (z * S);

        }

        private static double ApproximateInverseTailProbability (double P) {
            double zz = P * P;
            double log = Math.Log(2.0 / Math.PI / zz);
            double S = log - Math.Log(log);
            return (Math.Sqrt(S));
        }

        public double ApproximateInverseGamma (double P) {

            // compute  (Gamma(a+1) * P)^(1/a)
            double x0;
            if (ga > 0.0) {
                x0 = Math.Pow(a * ga * P, 1.0 / a);
            } else {
                x0 = Math.Exp((Math.Log(a * P) - ga) / a);
            }
            double x1 = x0 / (a + 1.0);
            if (x1 < 0.33) {
                // do left-tail if possible, because it is quite fast
                // this is just a term-by-term inversion of the power series
                // P(a,x) = gamma(a,x) / Gamma(a) = e^{-x} x^a / Gamma(a+1) [ 1 + x / (a+1) + ... ]
                return (x0 * (1.0 + x1 + (3.0 * a + 5.0) * x1 * x1 / (a + 2.0) / 2.0));
            } else {
                if (a >= 1.0) {
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
                    double Q = 1.0 - P;
                    if (Q == 0.0) {
                        return (Double.PositiveInfinity);
                    } else {
                        // compute -log(Gamma(a) * Q)
                        double log;
                        if (ga > 0.0) {
                            log = - Math.Log(ga * Q);
                        } else {
                            log = ga - Math.Log(Q);
                        }
                        return (log + (a - 1.0) * Math.Log(log));
                    }
                }
            }

        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");

            double x_min = 0.0;
            double x = ApproximateInverseGamma(P);
            double x_max = 10.0 * x;

            double dx_old = Double.MaxValue;

            Console.WriteLine("x0={0}", x);

            for (int i = 0; i < Global.SeriesMax; i++) {

                double F = LeftProbability(x) - P;


                if (F < 0.0) {
                    x_min = x;
                } else {
                    x_max = x;
                }

                double FP = ProbabilityDensity(x);

                double dx = -F / FP;

                double x_new = x + dx;

                if (x_new < x_min) {
                    x_new = (x_min + x) / 2.0;
                } else if (x_new > x_max) {
                    x_new = (x + x_max) / 2.0;
                } else if (2.0 * Math.Abs(dx) > Math.Abs(dx_old)) {
                    if (F < 0.0) {
                        x_new = (x_min + x) / 2.0;
                    } else {
                        x_new = (x + x_max) / 2.0;
                    }
                }

                Console.WriteLine("x_min={0} x_new={1} x_max={2}", x_min, x_new, x_max);

                if (x == x_new) return (x);

                dx_old = x_new - x;
                x = x_new;

            }

            throw new NonconvergenceException();
        }

    }
}
