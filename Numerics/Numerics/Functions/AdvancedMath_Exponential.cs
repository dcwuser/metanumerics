using System;
using System.Diagnostics;


namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {

        // this series is technically convergent everywhere, but it won't start converging until k ~ x, so it's best to use with small x
        // takes about 100 terms at x~40; it would be nice to move to the asymptotic expansion for lower x, but it fails to converge in that region
        private static double IntegralEi_Series (double x) {
            double dy = x;
            double y = EulerGamma + Math.Log(x) + dy;
            for (int k = 2; k < Global.SeriesMax; k++) {
                double y_old = y;
                dy = dy * x / k;
                y = y_old + dy / k;
                if (y == y_old) return (y);
            }
            throw new NonconvergenceException();
        }

        // An asymptotic series for Ei(x) useful for large x is
        //   Ei(x) = \frac{e^x}{x} \sum_{k=0}^{\infty} \frac{k!}{x^k} =
        //         = \frac{e^x}{x} \left( 1 + \frac{1}{x} + \frac{2!}{x^2} + \cdots \right) 
        // x needs to be pretty large for this to be useful, though. It takes about 40 terms to achieve full precision
        // at x ~ 40, and fails to converge to full precision for x <~ 40.

        private static double IntegralEi_Asymptotic (double x) {
            double dy = 1.0;
            double y = dy;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double y_old = y;
                dy = dy * k / x;
                y = y_old + dy;
                if (y == y_old) return (Math.Exp(x) / x * y);
            }
            throw new NonconvergenceException();
        }

        /// <summary>
        /// Computes the principal value of the exponential integral.
        /// </summary>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of Ei(x).</returns>
        /// <remarks>
        /// <para>The function Ei(x) appears in the evaluation of some indefinite integrals involving exponents and in
        /// number theory in the approximation li(x) = Ei(ln x) to the cumulative distribution of primes.</para>
        /// <para>It is related to the real part of the exponential integral for negative arguments by
        /// E<sub>1</sub>(-x &#177; i&#x3B5;) = -Ei(x) &#x2213; i&#x3C0;. </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        public static double IntegralEi (double x) {
            if (x < 0) throw new ArgumentOutOfRangeException("x");
            if (x < 40.0) {
                return (IntegralEi_Series(x));
            } else {
                return (IntegralEi_Asymptotic(x));
            }
        }

        /// <summary>
        /// Computes the exponential integral.
        /// </summary>
        /// <param name="n">The order parameter.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of E<sub>n</sub>(x).</returns>
        /// <remarks>
        /// <para>The exponential integral is defined as:</para>
        /// <img src="../images/EIntegral.png" />
        /// <para>It is related to the incomplete Gamma function for negative, integer shape parameters by &#x393;(-k, x) = Ei<sub>k+1</sub>(x) / x<sup>k</sup>.</para>
        /// <para>In hydrology, E<sub>1</sub>(x) is sometimes called the Well function.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        public static double IntegralE (int n, double x) {

            if (x < 0.0) throw new ArgumentOutOfRangeException("x");

            // special case x = 0
            if (x == 0.0) {
                if (n <= 1) {
                    return (Double.PositiveInfinity);
                } else {
                    return (1.0 / (n - 1));
                }
            }

            if (n < 0) {
                // negative n is expressible using incomplete Gamma
                return (AdvancedMath.Gamma(1 - n, x) / MoreMath.Pow(x, 1 - n));
            } else if (n == 0) {
                // special case n=0
                return (Math.Exp(-x) / x);
            } else if (x < 2.0) {
                // use series for x < 1
                return (IntegralE_Series(n, x));
            } else {
                // use continued fraction for x > 1
                return (IntegralE_ContinuedFraction(n, x));
            }
        }

        // a series for E_{n}(x) from NR
        // the series is a little weird because, when n > 1, there is a finite sum of terms before the infinite series begins
        // converges after about 20 terms for x~1, faster for smaller x
        private static double IntegralE_Series (int n, double x) {
            double f, df;
            if (n == 1) {
                f = 0.0;
                df = 1.0;
            } else {
                f = -1.0 / (1 - n);
                df = 1.0;
                for (int k = 1; k < (n - 1); k++) {
                    double f_old = f;
                    df = df * (-x / k);
                    f -= df / (k + 1 - n);
                    if (f == f_old) return (f);
                }
                df = df * (-x / (n - 1));
            }
            f += df * (Psi(n) - Math.Log(x));
            for (int k = n; k < Global.SeriesMax; k++) {
                double f_old = f;
                df = df * (-x / k);
                f -= df / (k - n + 1);
                if (f == f_old) return (f);
            }
            throw new NonconvergenceException();
        }

        // converges after about 100 terms for n~1, x~1; faster for larger n and x
        // we may need to deal with vanishing denominators
        private static double IntegralE_ContinuedFraction (int n, double x) {
            int a = 1;
            double b = x + n;
            double D = 1.0 / b;
            double Df = a / b;
            double f = 0.0 + Df;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double f_old = f;
                a = -k * (n + k - 1);
                b += 2.0;
                D = 1.0 / (b + a * D);
                Df = (b * D - 1.0) * Df;
                //Debug.WriteLine(String.Format("D={0} Df={1}", D, Df));
                f += Df;
                if (f == f_old) {
                    //Debug.WriteLine(String.Format("f={0} e={1} ef={2}", f, Math.Exp(-x), Math.Exp(-x) * f));
                    return (Math.Exp(-x) * f);
                }
            }
            throw new NonconvergenceException();
        }

        /// <summary>
        /// Computes the cosine integral.
        /// </summary>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of Ci(x).</returns>
        /// <remarks>
        /// <para>The cosine integral is defined as:</para>
        /// <img src="../images/CiIntegral.png" />
        /// <para>The cosine integral diverges logrithmically to negative inifity at the origin and executes a damped oscilation arround zero as its argument increases.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        public static double IntegralCi (double x) {
            if (x < 0.0) throw new ArgumentOutOfRangeException("x");
            if (x < 4.0) {
                return (IntegralCi_Series(x));
            } else {
                return (-IntegralE1_Imaginary_ContinuedFraction(x).Re);
            }
        }

        /// <summary>
        /// Computes the sine integral.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of Si(x).</returns>
        /// <remarks>
        /// <para>The sine integral is defined as:</para>
        /// <img src="../images/SiIntegral.png" />
        /// <para>The sine integral is zero at the origin and executes a damped oscilation arround &#x3C0;/2 as its argument increases.</para>
        /// </remarks>
        public static double IntegralSi (double x) {
            if (x < 0.0) {
                return (-IntegralSi(-x));
            } else if (x < 4.0) {
                return (IntegralSi_Series(x));
            } else {
                return (IntegralE1_Imaginary_ContinuedFraction(x).Im + Math.PI / 2.0);
            }
        }

        // converges to full accuracy in about 30 terms for x~4, less for smaller x
        private static double IntegralSi_Series (double x) {
            double xx = x * x;
            double dy = x;
            double y = dy;
            for (int k=3; k<Global.SeriesMax; k+= 2) {
                double y_old = y;
                dy = - dy * xx / k / (k - 1);
                y = y_old + dy / k;
                if (y == y_old) {
                    return (y);
                }
            }
            throw new NonconvergenceException();
        }

        // converges to full accuracy in about 30 terms for x~4, less for smaller x
        private static double IntegralCi_Series (double x) {
            double xx = x * x;
            double y = EulerGamma + Math.Log(x);
            double dy = 1.0;
            for (int k = 2; k < Global.SeriesMax; k += 2) {
                double y_old = y;
                dy = - dy * xx / k / (k - 1);
                y = y_old + dy / k;
                if (y == y_old) {
                    return (y);
                }
            }
            throw new NonconvergenceException();
        }

        // converges in about 30 terms for x~8, 50 terms for x~4, and 90 terms for x~2, always to full accuracy 
        private static Complex IntegralE1_Imaginary_ContinuedFraction (double x) {
            double a = 1.0;
            Complex b = new Complex(1.0,x);
            Complex D = 1.0 / b;
            Complex Df = a / b;
            Complex f = 0.0 + Df;
            for (int k = 1; k < Global.SeriesMax; k++) {
                Complex f_old = f;
                a = - k * k;
                b += 2.0;
                D = 1.0 / (b + a * D);
                Df = (b * D - 1.0) * Df;
                f += Df;
                if (f == f_old) {
                    return (new Complex(Math.Cos(x), -Math.Sin(x)) * f);
                }
            }
            throw new NonconvergenceException();

        }

        /// <summary>
        /// Computes the inverse tangent integral.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of Ti<sub>2</sub>(x).</returns>
        /// <remarks>
        /// <para>The inverse tanget integral is defined by:</para>
        /// <img src="..\images\IntegralTi.png" />
        /// </remarks>
        /// <seealso href="http://mathworld.wolfram.com/InverseTangentIntegral.html"/>
        public static double IntegralTi (double x) {
            if (x < 0.0) {
                return (-IntegralTi(-x));
            } else if (x < 0.75) {
                return (IntegralTi_Series(x));
            } else if (x < 4.0 / 3.0) {
                return (IntegralTi_LogSeries(x));
            } else {
                return (IntegralTi_Series(1.0 / x) + Math.PI / 2.0 * Math.Log(x));
            }
        }

        // We can get a series for Ti(x) by taking the series for \tan(x), reducing each power by one, and integrating term-by-term.
        // At x = 1, this series becomes simply the alternating sum of inverse squares, which we know converges to the Catalan
        // constant but does so very slowly. This series does not converge fast enough for computational purposes near x ~ 1.

        private static double IntegralTi_Series (double x) {

            double mx2 = -x * x;
            double df = x;
            double f = df;

            for (int k = 3; k < Global.SeriesMax; k += 2) {
                double f_old = f;
                df *= mx2;
                f += df / (k * k);
                if (f == f_old) return (f);
            }

            throw new NonconvergenceException();
        }

        private static double IntegralTi_LogSeries (double x) {
            x = Math.Log(x);
            double x2 = x * x;
            return (
                Catalan
                + Math.PI / 4.0 * x
                + 1.0 / 4.0 * x2
                - 1.0 / 48.0 * x2 * x2
                + 1.0 / 288.0 * x2 * x2 * x2
                - 61.0 / 80640.0 * x2 * x2 * x2 * x2
                + 277.0 / 1451520.0 * x2 * x2 * x2 * x2 * x2
                - 50521.0 / 958003200.0 * x2 * x2 * x2 * x2 * x2 * x2
                + 41581.0 / 2682408960.0 * x2 * x2 * x2 * x2 * x2 * x2 * x2
            );

        }

    }

    public static partial class AdvancedComplexMath {

        /// <summary>
        /// Computes the entire complex exponential integral.
        /// </summary>
        /// <param name="z">The complex argument.</param>
        /// <returns>The value of Ein(z).</returns>
        /// <remarks>
        /// <para>The entire exponential integral function can be defined by an integral or an equivalent series.</para>
        /// <para>Both Ei(z) and E<sub>1</sub>(x) and be obtained from Ein(z).</para>
        /// <para>Unlike either Ei(z) or E<sub>1</sub>(z), Ein(z) is entire, that is, it has no poles or cuts anywhere
        /// in the complex plane.</para>
        /// </remarks>
        public static Complex Ein (Complex z) {

            if (z.Re < -40.0) {
                // For sufficiently negative z, use the asymptotic expansion for Ei
                return (AdvancedMath.EulerGamma + ComplexMath.Log(-z) - IntegralEi_AsymptoticSeries(-z));
            } else {
                // Ideally, we would like to use the series within some simple radius of the origin and the continued fraction
                // outside it. That works okay for the right half-plane, but for z.Re < 0 we have to contend with the fact that the
                // continued fraction fails on or near the negative real axis. That makes us use the series in an oddly shaped
                // region that extends much further from the origin than we would like into the left half-plane. It would be
                // great if we could find an alternative approach in that region.
                if (IsEinSeriesPrefered(z)) {
                    return (Ein_Series(z));
                } else {
                    return (AdvancedMath.EulerGamma + ComplexMath.Log(z) + IntegeralE1_ContinuedFraction(z));
                }
            }

        }

        // We use it for |z| < 4 on the positive real axis, moving to |z| < 40 on the negative real axis.

        private static Complex Ein_Series (Complex z) {

            Complex zp = z;
            Complex f = zp;

            for (int k = 2; k < Global.SeriesMax; k++) {
                Complex f_old = f;
                zp *= -z / k;
                f += zp / k;
                if (f == f_old) return (f);
            }

            throw new NonconvergenceException();
        }

        // This continued fraction is valid for |z| >> 1, except close to the negative real axis.

        private static Complex IntegeralE1_ContinuedFraction (Complex z) {
            int a = 1;                  // a_1
            Complex b = z + 1.0;        // b_1
            Complex D = 1.0 / b;        // D_1 = 1 / b_1 (denominator of a_1 / b_1)
            Complex Df = a / b;         // Df_1 = f_1 - f_0 = a_1 / b_1
            Complex f = 0.0 + Df;       // f_1 = f_0 + Df_1 = b_0 + a_1 / b_1 (here b_0 = 0)
            for (int k = 1; k < Global.SeriesMax; k++) {
                Complex f_old = f;
                a = -k * k;
                b += 2.0;
                D = 1.0 / (b + a * D);
                Df = (b * D - 1.0) * Df;
                f += Df;
                if (f == f_old) {
                    return (ComplexMath.Exp(-z) * f);
                }
            }
            throw new NonconvergenceException();
        }

        private static Complex IntegralEi_AsymptoticSeries (Complex z) {

            Complex dy = 1.0;
            Complex y = dy;
            for (int k = 1; k < Global.SeriesMax; k++) {
                Complex y_old = y;
                dy = dy * k / z;
                y = y_old + dy;
                if (y == y_old) return (ComplexMath.Exp(z) / z * y);
            }
            throw new NonconvergenceException();

        }

        // This defines the region where, empirically and approximately, the number of terms required by the Ein series is
        // less than the number of terms required by the E1 continued fraction. We use piecewise linear and quadratic terms
        // to define it. Basically, it says stay away from the negative real axis, particularly for small negative real parts.

        private static bool IsEinSeriesPrefered (Complex z) {

            if (z.Re > 4.0) {
                return(false);
            } else if (z.Re > 0.0) {
                return (Math.Abs(z.Im) < 6.0 * (1.0 - MoreMath.Pow2(z.Re / 4.0)));
            } else if (z.Re > -20.0) {
                return (Math.Abs(z.Im) < 7.0 - Math.Abs(z.Re + 10.0) / 10.0);
            } else if (z.Re > -50.0) {
                return (Math.Abs(z.Im) < 10.0 + z.Re / 5.0);
            } else {
                return (false);
            }

        }

    }

}
