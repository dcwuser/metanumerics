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

        // an asymptotic series useful for large x
        // takes about 40 terms at x ~ 40 and fails to converge in double precision for x < 40
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

        // to do: what about finding the root in Ei(x) near x~0.37; very close to this root we loose accuracy because of
        // the finite number of digits we have for EulerGamma.

        /// <summary>
        /// Computes the principal value of the exponential integral.
        /// </summary>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of Ei(<paramref name="x"/>).</returns>
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
        /// <param name="n">The order parameter, which must be non-negative.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of Ei<sub><paramref name="n"/></sub>(<paramref name="x"/>).</returns>
        /// <remarks>
        /// <para>The exponential integral is defined as Ei<sub>n</sub>(x) = <sub>1</sub>&#x222B;<sup>&#x221E;</sup>dt e<sup>-xt</sup>/t<sup>n</sup>.</para>
        /// <para>It is related to the incomplete Gamma function for negative, integer shape parameters by &#x393;<sub>Q</sub>(-k, x) = Ei<sub>k+1</sub>(x) / x<sup>k</sup>.</para></remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        public static double IntegralE (int n, double x) {
            // should be possible to do n<0; this is just alpha(n,z)
            if (n < 0) throw new ArgumentOutOfRangeException("n");
            if (x < 0.0) throw new ArgumentOutOfRangeException("x");

            if (n == 0) {
                // special case n=0
                return (Math.Exp(-x) / x);
            } else if (x == 0.0) {
                // special case x=0
                if (n == 1) {
                    return (Double.PositiveInfinity);
                } else {
                    return (1.0 / (n - 1));
                }
            } else if (x < 1.0) {
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
        /// <para>The sine integral can be defined via Si(x) = <sub>0</sub>&#x222B;<sup>x</sup>dt sin(t)/t.</para>
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


        // polylog
        // Li_{n}(z) = \sum_{k=1}^{\infty} \frac{z^{k}}{k^{n}} = z 

    }
}
