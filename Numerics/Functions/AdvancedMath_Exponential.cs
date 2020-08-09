using System;
using System.Diagnostics;


namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {

        // Our methods for Si, Ci, Ei, Shi, and Chi are all from Stegun & Zucker, "Automatic Computing Methods for Special Functions. Part III.
        // The Sine, Cosine, Exponential Integrals, and Related Functions", Journal of Research of the National Bureau of Standard 80 B (1976) 291
        // https://nvlpubs.nist.gov/nistpubs/jres/80B/jresv80Bn2p291_A1b.pdf

        // It's remarkable that I can't find any better algorithms. I'm reasonably pleased with these, except for Ei, Shi, and Chi for x ~ 16-48,
        // for which we are near the limits of numerical convergence for the series used.

        // This series is technically convergent everywhere, but it won't start converging until k ~ x, so it's best to use with small x
        // takes about 100 terms at x~40; it would be nice to move to the asymptotic expansion for lower x, but it fails to converge in that region.

        private static double IntegralEi_Series (double x) {
            double dy = x;
            double y = dy;
            for (int k = 2; k < Global.SeriesMax; k++) {
                double y_old = y;
                dy *= x / k;
                y += dy / k;
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
                dy *= k / x;
                y += dy;
                if (y == y_old) {
                    return (Math.Exp(x) / x * y);
                }
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
        /// <para>Ei(x) is related to the real part of E<sub>1</sub>(x) for negative arguments by:</para>
        /// <img src="../images/E1EiRelation.png" />
        /// <para>To compute Ei(z) in the entire complex plane, use <see cref="AdvancedComplexMath.Ein(Complex)"/>.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        /// <seealso cref="IntegralE(int, double)"/>
        /// <seealso href="https://en.wikipedia.org/wiki/Exponential_integral"/>
        /// <seealso href="http://mathworld.wolfram.com/ExponentialIntegral.html"/>
        /// <seealso href="https://dlmf.nist.gov/6"/>
        public static double IntegralEi (double x) {
            if (x < 0.0) {
                throw new ArgumentOutOfRangeException(nameof(x));
            } else if (x < 40.0) {
                return (EulerGamma + Math.Log(x) + IntegralEi_Series(x));
            } else if (x < Double.PositiveInfinity) {
                return (IntegralEi_Asymptotic(x));
            } else if (Double.IsPositiveInfinity(x)) {
                return (Double.PositiveInfinity);
            } else {
                Debug.Assert(Double.IsNaN(x));
                return (x);
            }
        }

        /// <summary>
        /// Computes the generalized exponential integral.
        /// </summary>
        /// <param name="n">The order parameter.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of E<sub>n</sub>(x).</returns>
        /// <remarks>
        /// <para>The generalized exponential integral is defined as:</para>
        /// <img src="../images/EIntegral.png" />
        /// <para>It is related to the incomplete Gamma function (<see cref="Gamma(double, double)"/>)
        /// for negative, integer shape parameters.</para>
        /// <img src="../images/EnGammaRelation.png" />
        /// <para>For n=1, it expressible as a simple power series.</para>
        /// <img src="../images/E1Series.png" />
        /// <para>For negative x, E<sub>1</sub>(x) develops an imaginary part, but its real part is given by the Ei(x) function
        /// (<see cref="IntegralEi(double)"/>).</para>
        /// <img src="../images/E1EiRelation.png" />
        /// <para>To compute E<sub>1</sub>(z) in the entire complex plane, use <see cref="AdvancedComplexMath.Ein(Complex)"/>.</para>
        /// <para>Sometimes the function E<sub>1</sub>(z) is called the exponential integral, and sometimes that name is used
        /// for Ei(x). In hydrology, E<sub>1</sub>(x) is sometimes called the Well function.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        /// <seealso href="http://mathworld.wolfram.com/En-Function.html"/>
        public static double IntegralE (int n, double x) {

            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));

            // Special case x = 0.
            if (x == 0.0) {
                if (n <= 1) {
                    return (Double.PositiveInfinity);
                } else {
                    return (1.0 / (n - 1));
                }
            }

            // Special case n negative and zero.
            if (n < 0) {
                // negative n is expressible using incomplete Gamma
                return (AdvancedMath.Gamma(1 - n, x) / MoreMath.Pow(x, 1 - n));
            } else if (n == 0) {
                // special case n=0
                return (Math.Exp(-x) / x);
            }

            // Now we are sure x > 0 and n > 0.
            if (x < 2.0) {
                return (IntegralE_Series(n, x));
            } else if (x < expLimit) {
                // Since E_n(x) < e^{-x}, we can short-cut to zero if x is big enough.
                // This nicely avoids our continued fraction's bad behavior for infinite x.
                return (IntegralE_ContinuedFraction(n, x));
            } else if (x <= Double.PositiveInfinity) {
                return (0.0);
            } else {
                Debug.Assert(Double.IsNaN(x));
                return (x);
            }
        }

        // The value of x at which e^{-x} drops below the smallest representable double.
        private static readonly double expLimit = -Math.Log(Double.Epsilon);

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
                f += Df;
                if (f == f_old) {
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
        /// <para>The cosine integral diverges logarithmically to negative infinity at the origin and executes a damped oscillation around zero as its argument increases.</para>
        /// <para>To obtain the non-divergent part of the consine integral near the origin, use <see cref="IntegralCin"/>.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Trigonometric_integral"/>
        /// <seealso href="http://mathworld.wolfram.com/CosineIntegral.html"/>
        /// <seealso href="https://dlmf.nist.gov/6"/>
        public static double IntegralCi (double x) {
            if (x < 0.0) {
                throw new ArgumentOutOfRangeException(nameof(x));
            } else if (x < 4.0) {
                return (EulerGamma + Math.Log(x) + IntegralCi_Series(-x * x));
            } else if (x < Double.PositiveInfinity) {
                return (-IntegralE1_Imaginary_ContinuedFraction(x).Re);
            } else {
                return (Double.NaN);
            }
        }

        /// <summary>
        /// Computes the entire cosine intgral.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of Cin(x).</returns>
        /// <remarks>
        /// <para>The entire cosine integral Cin(x) is related to the conventional cosine integral Ci(x) by:</para>
        /// <img src="../images/IntegralCin.png"/>
        /// <para>Unlike Ci(x), Cin(x) is regular at the origin, and may therefore be more useful in applications
        /// that need the non-divergent part of a cosine integral. In fact, Cin(x) is entire, meaning it has no poles
        /// or cuts anywhere in the complex plane. But, unlike Ci(x), Cin(x) does diverge (logarithmicaly) for large x.</para>
        /// </remarks>
        /// <seealso href="https://en.wikipedia.org/wiki/Trigonometric_integral"/>
        /// <seealso href="https://dlmf.nist.gov/6"/>
        /// <seealso cref="IntegralCi"/>
        public static double IntegralCin (double x) {
            if (x < 0.0) {
                return (IntegralCin(-x));
            } else if (x < 4.0) {
                return (-IntegralCi_Series(-x * x));
            } else if (x < Double.PositiveInfinity) {
                return (EulerGamma + Math.Log(x) + IntegralE1_Imaginary_ContinuedFraction(x).Re);
            } else {
                return (Double.NaN);
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
        /// <para>The sine integral is zero at the origin and executes a damped oscillation around &#x3C0;/2 as its argument increases.</para>
        /// </remarks>
        /// <seealso href="https://en.wikipedia.org/wiki/Trigonometric_integral"/>
        /// <seealso href="http://mathworld.wolfram.com/SineIntegral.html"/>
        /// <seealso href="https://dlmf.nist.gov/6"/>
        public static double IntegralSi (double x) {
            if (x < 0.0) {
                return (-IntegralSi(-x));
            } else if (x < 4.0) {
                return (x * IntegralSi_Series(-x * x));
            } else if (x < Double.PositiveInfinity) {
                return (IntegralE1_Imaginary_ContinuedFraction(x).Im + Math.PI / 2.0);
            } else if (Double.IsPositiveInfinity(x)) {
                return (Math.PI / 2.0);
            } else {
                Debug.Assert(Double.IsNaN(x));
                return (x);
            }
        }

        /// <summary>
        /// Computes the hyperbolic sine integral.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of Shi(x).</returns>
        /// <seealso href="https://en.wikipedia.org/wiki/Trigonometric_integral#Hyperbolic_sine_integral"/>
        /// <seealso href="https://mathworld.wolfram.com/Shi.html"/>
        public static double IntegralShi (double x) {
            if (x < 0.0) {
                return (-IntegralShi(-x));
            } else if (x < 40.0) {
                return (x * IntegralSi_Series(x * x));
            } else if (x <= Double.PositiveInfinity) {
                // For x > 40, Shi(x) differs from Ei(x)/2 only at the ~ 1.0E-34
                // level, far below double precision.
                return (0.5 * IntegralEi_Asymptotic(x));
            } else {
                Debug.Assert(Double.IsNaN(x));
                return (x);
            }
            // I would really like a better solution in the x ~ 16-48 range.
            // x ~ 40 is the very limit of convergence for asymptotic series,
            // and power series needs more than 100 terms to converge at that
            // point.
        }

        // This is written so as to be usable for both Si and Shi.
        // Converges to full accuracy in about 30 terms for x~4, less for smaller x
        private static double IntegralSi_Series (double xSquared) {
            double dy = 1.0;
            double y = dy;
            for (int k=3; k<Global.SeriesMax; k+= 2) {
                double y_old = y;
                dy *= xSquared / (k * (k - 1));
                y += dy / k;
                if (y == y_old) {
                    return (y);
                }
            }
            throw new NonconvergenceException();
        }

        // converges to full accuracy in about 30 terms for x~4, less for smaller x
        private static double IntegralCi_Series (double xSquared) {
            double dy = 1.0;
            double y = 0.0;
            for (int k = 2; k < Global.SeriesMax; k += 2) {
                double y_old = y;
                dy *= xSquared / (k * (k - 1));
                y += dy / k;
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
                    return (new Complex(MoreMath.Cos(x), -MoreMath.Sin(x)) * f);
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
        /// <para>The inverse tangent integral is defined by:</para>
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
        /// <img src="..\images\EinIntegralSeries.png" />
        /// <para>Both Ei(x) and E<sub>1</sub>(z) can be obtained from Ein(z).</para>
        /// <img src="..\images\E1EiEinRelation.png" />
        /// <para>Unlike either Ei(x) or E<sub>1</sub>(z), Ein(z) is entire, that is, it has no poles or cuts anywhere
        /// in the complex plane.</para>
        /// </remarks>
        /// <seealso cref="AdvancedMath.IntegralEi(double)"/>
        /// <seealso cref="AdvancedMath.IntegralE(int, double)"/>
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
                return (Math.Abs(z.Im) < 6.0 * (1.0 - MoreMath.Sqr(z.Re / 4.0)));
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
