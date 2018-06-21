using System;
using System.Diagnostics;


namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {


        /// <summary>
        /// Computes the error function.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of erf(x).</returns>
        /// <remarks>
        /// <para>The error function can be defined as a Gaussian integral.</para>
        /// <img src="../images/ErfIntegral.png" />
        /// <para>The area under a bell curve (<see cref="Meta.Numerics.Statistics.Distributions.NormalDistribution"/>) within &#x2213;z
        /// standard deviations of the mean is given by erf(z/&#x221A;2).</para>
        /// <para>For large values of x, erf(x) &#x2248; 1 to within floating-point accuracy. To obtain accurate values of erfc(x) = 1 - erf(x)
        /// in this range, use the <see cref="Erfc" /> function.</para>
        /// <para>The inverse of the error function is available as <see cref="AdvancedMath.InverseErf"/>.</para>
        /// <para>Values of the error function for complex arguments can be obtained using <see cref="AdvancedComplexMath.Erf"/> method, or using the
        /// equivalent but re-parameterized <see cref="AdvancedComplexMath.Faddeeva"/> function.</para>
        /// </remarks>
        /// <seealso cref="Erfc"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Error_function" />
        /// <seealso href="http://mathworld.wolfram.com/Erf.html" />
        public static double Erf (double x) {
            if (x < -1.5) {
                return (Erfc_ContinuedFraction(-x) - 1.0);
            } else if (x < 1.5) {
                // Use the series near the origin.
                return (Erf_Series(x));
            } else if (x <= Double.PositiveInfinity) {
                return (1.0 - Erfc_ContinuedFraction(x));
            } else {
                return (Double.NaN);
            }
        }

        /// <summary>
        /// Computes the complementary error function.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of erfc(<paramref name="x"/>) = 1 - erf(<paramref name="x"/>).</returns>
        /// <remarks>
        /// <para>The complementary error function can be defined by an integral, or by its relation to the error function.</para>
        /// <img src="../images/ErfcIntegral.png" />
        /// <para>The area in the tails of a bell curve ((<see cref="Meta.Numerics.Statistics.Distributions.NormalDistribution"/>)
        /// beyond &#x2213;z standard deviations from the mean is given by erfc(z/&#x221A;2).</para>
        /// <para>For small values of x, erfc(x) &#x2248; 1 to within floating-point accuracy. To obtain accurate values of erf(x) = 1 - erfc(x)
        /// in this region, use the <see cref="Erf" /> function.</para>
        /// <para>The inverse of the complementary error function is available as <see cref="InverseErfc(double)"/>.</para>
        /// <para>Values of the complementary error function for complex arguments can be obtained using the
        /// <see cref="AdvancedComplexMath.Faddeeva(Complex)"/> function.</para>
        /// </remarks>
        /// <seealso cref="Erf" />
        /// <seealso href="http://mathworld.wolfram.com/Erfc.html"/>
        public static double Erfc (double x) {
            if (x < -1.5) {
                return (2.0 - Erfc_ContinuedFraction(-x));
            } else if (x < 1.5) {
                return (1.0 - Erf_Series(x));
            } else if (x <= Double.PositiveInfinity) {
                return (Erfc_ContinuedFraction(x));
            } else {
                return (Double.NaN);
            }
        }

        // This is a straightforward implementation of the error series:
        //  erf(x) = \frac{2x}{\sqrt{\pi}} \sum_{k=0}{^\infty} \frac{(-1)^k x^{2k}}{k! (2k+1)}
        // Converges in about a dozen terms for x less than about 1, but still converges for higher x.

        // This is just the incomplete gamma series specialized for a -> 1/2, x -> x^2, but the specialization
        // speeds it up considerably, so we do it explicitly.

        private static double Erf_Series (double x) {
            double mx2 = -x * x;
            double t = 1.0;
            double s = t;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double s_old = s;
                t *= mx2 / k;
                s += t / (2 * k + 1);
                if (s == s_old) {
                    return (2.0 / Global.SqrtPI * x * s);
                }
            }
            throw new NonconvergenceException();
        }

        private static double Erfc_ContinuedFraction (double x) {

            Debug.Assert(x > 0.0);

            // Erfc underflows for x larger than about 27, so no need to compute continued fraction.
            // This check also avoids computing Infinity * 0 = NaN if x is Infinity.
            if (x > 28.0) return (0.0);

            // Compute the continued fraction.
            double x2 = x * x;
            double aa = 1.0;			// a_1
            double bb = x2 + 0.5;   	// b_1
            double D = 1.0 / bb;		// D_1 = b_0/b_1
            double Df = aa / bb;		// Df_1 = f_1 - f_0
            double f = 0.0 + Df;		// f_1 = f_0 + Df_1 = b_0 + Df_1
            for (int k = 1; k < Global.SeriesMax; k++) {
                double f_old = f;
                aa = -k * (k - 0.5);
                bb += 2.0;
                D = 1.0 / (bb + aa * D);
                Df = (bb * D - 1.0) * Df;
                f += Df;
                if (f == f_old) {
                    return (x / Global.SqrtPI * Math.Exp(-x2) * f);
                }
            }
            throw new NonconvergenceException();
        }

        /// <summary>
        /// Computes the inverse complementary error function.
        /// </summary>
        /// <param name="y">The value of erfc(x), which must lie between 0 and 1.</param>
        /// <returns>The corresponding argument x.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="y"/> lies outside [0,1].</exception>
        /// <seealso cref="Erf"/>
        /// <seealso cref="Erfc" />
        /// <seealso cref="InverseErf" />
        public static double InverseErfc (double y) {

            if ((y < 0.0) || (y > 1.0)) throw new ArgumentOutOfRangeException(nameof(y));

            if (y == 0.0) {
                return (Double.PositiveInfinity);
            } else if (y < 0.5) {
                return (InverseErfcByRefinement(y));
            } else {
                return (InverseErfSeries(1.0 - y));
            }
            
        }

        /// <summary>
        /// Computes the inverse error function.
        /// </summary>
        /// <param name="y">The error function value erf(x), which must lie between -1 and 1.</param>
        /// <returns>The corresponding argument x.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="y"/> is outside [-1,1].</exception>
        /// <seealso cref="Erf"/>
        /// <seealso cref="Erfc" />
        /// <seealso cref="InverseErfc" />
        public static double InverseErf (double y) {

            if (y < 0.0) {
                return (-InverseErf(-y));
            } else if (y < 0.5) {
                return (InverseErfSeries(y));
            } else if (y < 1.0) {
                return (InverseErfcByRefinement(1.0 - y));
            } else if (y == 1.0) {
                return (Double.PositiveInfinity);
            } else if (Double.IsNaN(y)) {
                return (Double.NaN);
            } else {
                throw new ArgumentOutOfRangeException(nameof(y));
            }

        }

        // There is no cancelation in this series, so accuracy is just a matter of how many terms we are willing to take.
        // For x ~ 0.1, 8 terms; for x ~ 0.25, 12 terms; for x ~ 0.5, 24 terms.

        // \erf^{-1}(x) = \sum_{k=0}^{\infty} \frac{c_k}{2k+1} \left( \frac{\sqrt{\pi}{2} x \right)^{2k+1}

        // c_k = \sum_{j=0}^{k-1} \frac{c_{j} c_{k-1-k}}{(j+1)(2j+1)}

        private static double InverseErfSeries (double x) {

            double z = Global.SqrtPI / 2.0 * x;
            double z2 = z * z;

            double s = 1.0;
            double z2k = 1.0;
            for (int k = 1; k < inverfSeriesCoefficients.Length; k++) {
                double s_old = s;
                z2k *= z2;
                s += inverfSeriesCoefficients[k] * z2k;
                if (s == s_old) { return (z * s); }
            }

            throw new NonconvergenceException();
        }

        private static readonly double[] inverfSeriesCoefficients = ComputeInverseErfSeriesCoefficients(24);

        private static double[] ComputeInverseErfSeriesCoefficients (int n) {
            double[] d = new double[n + 1];
            d[0] = 1.0;
            for (int k = 0; k < n; k++) {
                for (int j = 0; j <= k; j++) {
                    d[k + 1] += d[j] * d[k - j] / (j + 1) / (2 * j + 1);
                }
            }
            for (int k = 1; k <= n; k++) {
                d[k] /= (2 * k + 1);
            }
            return (d);
        }

        private static double InverseErfcByRefinement (double y) {

            // First we get an approximation to x using either the asymptotic expansion (good only for very small values)
            // or via a rational approximation. In either case, the cost is one log, plus one square root, plus a handfull of additional flops,
            // and the result is accurate to about 10^{-6}.

            double x;
            if (y < 1.0E-5) {
                x = InverseErfcAsymptoticExpansion(y);
            } else if (y < 0.55) {
                x = InverseErfcRationalApproximation(y);
            } else {
                throw new InvalidOperationException();
            }

            // Then we refine via Halley's method (http://en.wikipedia.org/wiki/Halley%27s_method), an extension of Newton's method.

            // Given f(x) = f0 + f1 x, Newton's method sets f(x) = 0 to obtain x = - f0 / f1.

            // Given f(x) = f0 + f1 x + 1/2 f2 x^2, Halley's method sets f(x) = 0 to obtain x = - f0 / (f1 + 1/2 f2 x) and uses the Netwon result
            // for x in the RHS to obtain x = - f0 / (f1 - 1/2 f2 f0 / f1).

            // Sufficiently close to the root, Newton's method is quadratic (significant digits double each iteration) and Halley's method is cubic (significant
            // digits triple each iteration). So an initial estimate with precision 10^{-4} should get us to full precision in two iterations.

            // For f_0 = erfc(x) - y, f_1 = - \frac{2}{\sqrt{\pi}} e^{-x^2}, and f_2 = \frac{4 x}{\sqrt{\pi}} e^{-x^2} = -2 x f_1.
            // Defining r = - f0 / f1 as the Newton estimate, the Halley estimate is simply r / (1 - x r). Since the second
            // derivative is so simple to compute from the first derivative, a Halley step is only trivially more costly than a Newton step.

            for (int i = 0; i < 2; i++) {
                double z = Erfc(x) - y;
                double r = z * Global.SqrtPI * Math.Exp(x * x) / 2.0;
                double dx = r / (1.0 - x * r);
                x += dx;
             }

            // We hard-code exactly two steps instead of testing for convergence because the iteration can get into loops of tiny changes instead
            // of producing a zero change, and we know that two iterations will be sufficient.

            return (x);

        }

        // Blair et al., Mathematics of Computation 30 (1976) 827 write that "by inverting the standard asymptotic series" for erf x in terms of 1/x
        // "we can derive an asymptotic expansion for inverf" and then give the expansion encoded here. Their expansion is the most accurate I've
        // found, and I would love to derive more terms, but I can't figure out how they derived it. The DLMF (http://dlmf.nist.gov/7.17) says
        // their expansion 7.17.3 "follows from Blair et al. (1976), after modifications" but they must have made an error in their modifications
        // because their accuracy sucks compared to Blair et al.

        private static double InverseErfcAsymptoticExpansion (double x) {

            double y = -Math.Log( Global.SqrtPI * x);
            double lny = Math.Log(y);

            double s = y - lny / 2.0;
            s += 1.0 / y * (lny / 4.0 - 1.0 / 2.0);
            double y2 = y * y;
            double lny2 = lny * lny;
            s += 1.0 / y2 * (lny2 / 16.0 - 3.0 / 8.0 * lny + 7.0 / 8.0);
            double y3 = y2 * y;
            double lny3 = lny2 * lny;
            s += 1.0 / y3 * (lny3 / 48.0 - 7.0 / 32.0 * lny2 + 17.0 / 16.0 * lny - 107.0 / 48.0);
            double y4 = y3 * y;
            double lny4 = lny3 * lny;
            s += 1.0 / y4 * (lny4 / 128.0 - 23.0 / 192.0 * lny3 + 29.0 / 32.0 * lny2 - 31.0 / 8.0 * lny + 1489.0 / 192.0);

            return (Math.Sqrt(s));

        }

        private static double InverseErfcRationalApproximation (double x) {

            // minimax (3,3) rational polynomial approximation
            // in region 3/4 < y < 5, corresponding to about 3.7E-6 < x < 0.75
            // with error less ~9E-7 throughout

            double y = Math.Sqrt(-2.0 * Math.Log(x));

            const double a0 = -0.0071034591712238851203;
            const double a1 = 0.047920059819193711592;
            const double a2 = 0.30517720195408962758;
            const double a3 = 0.48856023932132526440;

            const double b1 = 0.61364516662581878617;
            const double b2 = 0.67737432169146905075;
            const double b3 = 0.00056224293678634100951;

            return ((a0 + (a1 + (a2 + a3 * y) * y) * y) / (1.0 + (b1 + (b2 + b3 * y) * y) * y));

        }

        internal static double Probit (double P, double Q) {
            Debug.Assert(P + Q == 1.0);
            if (P < 0.25) {
                if (P == 0.0) {
                    return(Double.NegativeInfinity);
                } else {
                    return (-Global.SqrtTwo * InverseErfcByRefinement(2.0 * P));
                }
            } else if (P < 0.75) {
                return (Global.SqrtTwo * InverseErfSeries(2.0 * P - 1.0));
            } else {
                if (Q == 0.0) {
                    return (Double.PositiveInfinity);
                } else {
                    return (Global.SqrtTwo * InverseErfcByRefinement(2.0 * Q));
                }
            }

        }

        internal static double ApproximateInverseErf (double y) {
            double x = Global.SqrtPI * y / 2.0;
            double xx = x * x;
            double S = 1.0 + xx / 3.0 + 7.0 / 30.0 * xx * xx + 127.0 / 630.0 * xx * xx * xx;
            return (x * S);

        }

        internal static double ApproximateInverseErfc (double y) {
            double yy = y * y;
            double log = Math.Log(2.0 / Math.PI / yy);
            double S = log - Math.Log(log);
            return (Math.Sqrt(S / 2.0));
        }

        /// <summary>
        /// Computes the Dawson integral.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of F(x).</returns>
        /// <remarks>
        /// <para>The Dawson function is defined by the integral:</para>
        /// <img src="../images/DawsonIntegral.png" />
        /// <para>It is related to the error function for purely imaginary arguments.</para>
        /// <img src="../images/DawsonErfRelation.png" />
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Dawson_function"/>
        public static double Dawson (double x) {
            if (x < 0.0) {
                // Dawson is an odd function
                return (-Dawson(-x));
            } else if (x < 1.0) {
                // use the series expansion near the origin
                return (Dawson_Series(x));
            } else if (x < 10.0) {
                // use the Rybicki algorithm in between
                return (Dawson_Rybicki(x));
            } else if (x <= Double.PositiveInfinity) {
                // use the asymptotic expansion for large values
                return (Dawson_Asymptotic(x));
            } else {
                return (Double.NaN);
            }
        }

        // a series expansion for the Dawson integral near the origin
        // requires ~30 terms at x~2, ~20 terms at x~1, ~10 terms at x~0.5

        private static double Dawson_Series (double x) {
            double xx = -2.0 * x * x;
            double df = x;
            double f = df;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double f_old = f;
                df = df * xx / (2 * k + 1);
                f += df;
                if (f == f_old) {
                    return (f);
                }
            }
            throw new NonconvergenceException();
        }

        // as asymptotic expansion for the Dawson integral from large values
        // requires ~5 terms at x~50, ~15 terms at x~10, fails to converge below x~7

        private static double Dawson_Asymptotic (double x) {
            double xx = 2.0 * x * x;
            double df = 0.5 / x;
            double f = df;
            for (int k = 0; k < Global.SeriesMax; k++) {
                double f_old = f;
                df = df * (2 * k + 1) / xx;
                f += df;
                if (f == f_old) {
                    return (f);
                }
            }
            throw new NonconvergenceException();
        }

        // a rather strange algorithm by Rybicki, described in Numerical Recipes
        // it expresses Dawson's integral as a sum over exponentially suppressed terms with a suppression factor h
        // his result derives essentially from doing the integral numerically with a step size h
        // it converges in the limit h->0 and the error for non-zero h goes like e^{-(pi/2h)^2}; since the error
        // goes down exponentially with h, even quite moderate value of h give good results: we use h~0.25 to get full precision
        // the series is infinite, but since the terms die off exponentially you don't need many to get full precision: we need ~12 x 2

        private static double Dawson_Rybicki (double x) {

            int n0 = 2 * ((int)Math.Round(x / Dawson_Rybicki_h / 2.0));
            double x0 = n0 * Dawson_Rybicki_h;
            double y = x - x0;

            double f = 0.0;
            double b = Math.Exp(2.0 * Dawson_Rybicki_h * y);
            double bb = b * b;
            for (int k = 0; k < Dawson_Rybicki_coefficients.Length; k++) {
                double f_old = f;
                int m = 2 * k + 1;
                double df = Dawson_Rybicki_coefficients[k] * (b / (n0 + m) + 1.0 / b / (n0 - m));
                f += df;
                if (f == f_old) {
                    return (Math.Exp(-y * y) / Global.SqrtPI * f);
                }
                b = b * bb;
            }

            throw new NonconvergenceException();

        }

        private const double Dawson_Rybicki_h = 0.25;
        private static readonly double[] Dawson_Rybicki_coefficients = Compute_Dawson_Rybicki_Coefficients(0.25, 16);

        // pre-computes e^{-(h m)^2} for Rybicki algorithm

        private static double[] Compute_Dawson_Rybicki_Coefficients (double h, int n) {
            double[] coefficients = new double[n];
            for (int k = 0; k < n; k++) {
                int m = 2 * k + 1;
                double z = h * m;
                coefficients[k] = Math.Exp(-z * z);
            }
            return(coefficients);
        }


        // Fresnel integrals


        /// <summary>
        /// Computes the Fresnel cosine integral.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of C(x).</returns>
        /// <remarks>
        /// <para>The Fresnel cosine integral is defined as:</para>
        /// <img src="../images/FresnelCIntegral.png" />
        /// <para>It appears in wave optics in the calculation of diffraction patterns.</para>
        /// </remarks>
        /// <seealso cref="FresnelS"/>
        /// <seealso cref="Fresnel"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Fresnel_integral"/>
        public static double FresnelC (double x) {
            if (x < 0.0) {
                return (-FresnelC(-x));
            } else if (x < 2.0) {
                return (FresnelC_Series(x));
            } else if (x < fresnelLimit) {
                return (Fresnel_ContinuedFraction(x).Re);
            } else if (x <= Double.PositiveInfinity) {
                return (0.5);
            } else {
                return (Double.NaN);
            }
        }

        /// <summary>
        /// Computes the Fresnel sine integral.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of S(x).</returns>
        /// <remarks>
        /// <para>The Fresnel sine integral is defined as:</para>
        /// <img src="../images/FresnelSIntegral.png" />
        /// <para>It appears in wave optics in the calculation of diffraction patterns.</para>
        /// </remarks>
        /// <seealso cref="FresnelC"/>
        /// <seealso cref="Fresnel"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Fresnel_integral"/>
        public static double FresnelS (double x) {
            if (x < 0.0) {
                return (-FresnelS(-x));
            } else if (x < 2.0) {
                return (FresnelS_Series(x));
            } if (x < fresnelLimit) {
                return (Fresnel_ContinuedFraction(x).Im);
            } else if (x <= Double.PositiveInfinity) {
                return (0.5);
            } else {
                return (Double.NaN);
            }
        }

        /// <summary>
        /// Computes the Fresnel cosine and sine integrals.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value C(x) + i S(x).</returns>
        /// <remarks>
        /// <para>If you need both C(x) and S(x), it is more efficient to call this function than
        /// to call <see cref="FresnelC(double)"/> and <see cref="FresnelS(double)"/> seperately.</para>
        /// <para>Our definition of C(x) and S(x) agrees with that of Abromowitz and Stegun, the NIST Library of
        /// Mathematical Functions, Numerical Recipies, Mathematica, Matlab, SciPy, and any every other
        /// reference book and numerical software library I can find. Wikipedia, however uses an
        /// alternative definition.</para>
        /// <para>A plot of all values of this function in the complex plane as x ranges from
        /// negative infinity to positive infinity is called a Cornu spiral.</para>
        /// <para>The Fresnel function can be related to the complex error function along the line (1-&#x1D456;).</para>
        /// <img src="../images/FresnelErfRelation.png" />
        /// </remarks>
        /// <seealso cref="FresnelS"/>
        /// <seealso cref="FresnelC"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Fresnel_integral"/>
        /// <seealso href="http://mathworld.wolfram.com/CornuSpiral.html"/>
        public static Complex Fresnel (double x) {
            if (x < 0.0) {
                return (-Fresnel(-x));
            } else if (x < 2.0) {
                return (new Complex(FresnelC_Series(x), FresnelS_Series(x)));
            } else if (x < fresnelLimit) {
                return (Fresnel_ContinuedFraction(x));
            } else if (x <= Double.PositiveInfinity) {
                return (new Complex(0.5, 0.5));
            } else {
                return (Double.NaN);
            }
        }

        // For large x, C ~ \frac{1}{2} + \frac{1}{\pi x} \sin(\pi x^2 /2 ), S ~ frac{1}{2} - \frac{1}{\pi x} \cos(\pi x^2 / 2)
        // so if x is large enough that \frac{1}{\pi x} is indistinguishable from \frac{1}{2}, we can just return 1/2.
        
        // Even if we didn't have this limit, say because we were trying to compute C - 1/2 and S - 1/2, we would
        // have the problem that both the continued fraction and the direct asymptotic expansion require computing
        // trig functions of \pi x^2 /2, and this overflows for x larger than about 1.0E154.

        private const double fresnelLimit = 2.0 / Math.PI / 1.0E-16;

        // At one point, we used a Fresnel asymptotic series for x > 40, but it proved unnecessary.
        // It's neither faster nor more accurate than the continued fraction, and it can be used as
        // close to the origin.

        // C(z) = z \sum_{n=0}^{\infty} \frac{(-1)^n (\pi z^2 / 2)^{2n}}{(2n)! (4n + 1)}
        // Series requires ~10 terms at x~1, ~20 terms at x~2.

        private static double FresnelC_Series (double x) {
            double xx = Math.PI / 2.0 * x * x;
            double x4 = -xx * xx;
            double df = 1.0;
            double f = df;
            for (int two_k = 2; two_k < Global.SeriesMax; two_k += 2) {
                double f_old = f;
                df *= x4 / (two_k * (two_k - 1));
                f += df / (2 * two_k + 1);
                if (f == f_old) return (x * f);
            }
            throw new NonconvergenceException();
        }

        // S(z) = \frac{\pi z^3}{2} \sum_{n=0}^{\infty} \frac{(-1)^n (\pi z^2 / 2)^{2n}}{(2n + 1)! (4n + 3)}

        private static double FresnelS_Series (double x) {
            double xx = Math.PI / 2.0 * x * x;
            double x4 = -xx * xx;
            double df = 1.0;
            double f = df / 3.0;
            for (int two_k = 2; two_k < Global.SeriesMax; two_k += 2) {
                double f_old = f;
                df *= x4 / ((two_k + 1) * two_k);
                f += df / (2 * two_k + 3);
                if (f == f_old) return (x * xx * f);
            }
            throw new NonconvergenceException();
        }


        // Numerical Recipes quotes a continued fraction for C and S, derived frothe relation between C, S,
        // and the complex error function:
        //   C(x) + i S(x) = \frac{1 + i}{2} \erf(z) \quad z = \frac{\sqrt{\pi}}{2} (1 - i) x
        // and the continued fraction for erfc(z):
        //   \erfc(z) = \frac{2z}{\sqrt{\pi}} e^{-z^2} \left[
        //     \frac{1}{2z^2 + 1 -} \frac{1 \cdot 2}{\frac{2z^2 + 5 -} \frac{3 \cdot 4}{2z^2 + 9 -} \cdots
        //   \right]

        // Given the particular form of z, we can simplify this significantly. First note
        //   2z^2 = -i \pi x^2
        // Next note
        //   e^{-z^2} = e^{-i \pi x^2 / 2} = \cos(\pi x^2 / 2) - i \sin(\pi x^2 / 2)
        // So when we determine the continued fraction, we have
        //   \frac{1 + i}{2} \erf(z) = \frac{1 + i}{2} \left[1 - (1 - i) x (a + i b) \right]
        //     = \left( \frac{1}{2} - a x \right) + i \left( \frac{1}{2} - b x \right)

        // ~35 terms at x~2, ~12 terms at x~4, ~6 terms at x~8, ~3 terms for x~12 and above


        private static Complex Fresnel_ContinuedFraction (double x) {

            double xSquared = x * x;

            // It appears that (1) imaginary part of f doesn't change and
            // (2) real part of f merely increases at constant rate until "right" number is reached.
            // This is probably an indication that we could simplify the complex arithmetic further.

            // Use Lenz's method for the continued fraction.
            Complex a = 1.0;
            Complex b = new Complex(1.0, -Math.PI * xSquared);
            Complex D = 1.0 / b;
            Complex Df = a / b;
            Complex f = Df;
            for (int two_k = 2; two_k < Global.SeriesMax ; two_k += 2) {
                Complex f_old = f;
                a = -(two_k - 1) * two_k;
                b += 4.0;
                D = 1.0 / (b + a * D);
                Df = (b * D - 1.0) * Df;
                f += Df;
                if (f == f_old) {
                    double half_x_squared = 0.5 * xSquared;
                    Complex e = new Complex(MoreMath.CosPi(half_x_squared), MoreMath.SinPi(half_x_squared));
                    Complex g = e * f * x;
                    return new Complex(0.5 - g.Re, 0.5 - g.Im);
                }
            }
            throw new NonconvergenceException(); 
        }

    }


    public static partial class AdvancedComplexMath {

        /// <summary>
        /// Computes the complex Faddeeva function.
        /// </summary>
        /// <param name="z">The complex argument.</param>
        /// <returns>The complex value of w(z).</returns>
        /// <remarks>
        /// <para>The Faddeeva function w(z) is related to the error function with a complex argument.</para>
        /// <img src="../images/FaddeevaErfcRelation.png" />
        /// <para>It also has an integral representation.</para>
        /// <img src="../images/FaddeevaIntegral.png" />
        /// <para>For purely imaginary values, it reduces to the complementary error function (<see cref="AdvancedMath.Erfc"/>).
        /// For purely real values, it reduces to Dawson's integral (<see cref="AdvancedMath.Dawson"/>).</para>
        /// <para>It appears in the computation of the Voigt line profile function V(x;&#x3C3;,&#x3B3;).</para>
        /// <img src="../images/Voigt.png" />
        /// <para>Near the origin, w(z) &#x2248; 1. To accurately determine w(z) - 1 in this region, use the <see cref="Erf"/>
        /// function. Away from the origin near the large negative imaginary axis, the magnitude w(z) increases rapidly and
        /// may overflow.</para>
        /// <para>The image below shows the complex Faddeeva function near the origin, using domain coloring.</para>
        /// <img src="../images/ComplexFaddeevaPlot.png" />
        /// </remarks>
        /// <seealso cref="AdvancedComplexMath.Erf"/>
        /// <seealso cref="AdvancedMath.Erf" />
        /// <seealso cref="AdvancedMath.Erfc" />
        /// <seealso cref="AdvancedMath.Dawson"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Voigt_profile" />
        public static Complex Faddeeva (Complex z) {

            // use reflection formulae to ensure that we are in the first quadrant
            if (z.Im < 0.0) return (2.0 * ComplexMath.Exp(-ComplexMath.Sqr(z)) - Faddeeva(-z));
            if (z.Re < 0.0) return (Faddeeva(-z.Conjugate).Conjugate);

            double r = ComplexMath.Abs(z);
            if (r < 2.0) {
                // use series for small z
                return (ComplexMath.Exp(-ComplexMath.Sqr(z)) * (1.0 - Erf_Series(-ComplexMath.I * z)));
                //return (Faddeeva_Series(z));
            } else if ((z.Im < 0.1) && (z.Re < 30.0)) {
                // this is a special, awkward region
                // along the real axis, Re{w(x)} ~ e^{-x^2}; the Weideman algorithm doesn't compute this small number
                // well and the Laplace continued fraction misses it entirely; therefore very close to the real axis
                // we will use an analytic result on the real axis and Taylor expand to where we need to go.
                // unfortunately the Taylor expansion converges poorly for large x, so we drop this work-around near x~30,
                // when this real part becomes too small to represent as a double anyway
                double x = z.Re;
                double y = z.Im;
                return (Faddeeva_Taylor(new Complex(x, 0.0),
                                        Math.Exp(-x * x) + 2.0 * AdvancedMath.Dawson(x) / Global.SqrtPI * ComplexMath.I,
                                        new Complex(0.0, y)));
            } else if (r > 7.0) {
                // use Laplace continued fraction for large z
                return (Faddeeva_ContinuedFraction(z));
            } else {
                // use Weideman algorithm for intermediate region
                return (Faddeeva_Weideman(z));
            }
        }

        // there is a zero near (1.9915,-1.3548) -- check it out

        private static Complex Faddeeva_Weideman (Complex z) {
            Complex ZN = Faddeeva_Weideman_L + ComplexMath.I * z;
            Complex ZD = Faddeeva_Weideman_L - ComplexMath.I * z;
            Complex ZQ = ZN / ZD;
            Complex f = Faddeeva_Weideman_Coefficients[40];
            for (int k = 39; k > 0; k--) {
                f = f * ZQ + Faddeeva_Weideman_Coefficients[k];
            }
            Complex ZP = ZN * ZD;
            return (2.0 / ZP * f * ZQ + 1.0 / Global.SqrtPI / ZD);
        }

        private static readonly double Faddeeva_Weideman_L = Math.Sqrt(40.0 / Math.Sqrt(2.0));
        private static readonly double[] Faddeeva_Weideman_Coefficients = {
            3.0005271472811341147438, // 0
            2.899624509389705247492,
            2.616054152761860368947,
            2.20151379487831192991,
            1.725383084817977807050,
            1.256381567576513235243, // 5
            0.847217457659381821530,
            0.52665289882770863869581,
            0.2998943799615006297951,
            0.1550426380247949427170,
            0.0718236177907433682806, // 10
            0.0292029164712418670902,
            0.01004818624278342412539,
            0.002705405633073791311865,
            0.000439807015986966782752,
            -0.0000393936314548956872961, // 15
            -0.0000559130926424831822323,
            -0.00001800744714475095715480,
            -1.066013898494714388844e-6,
            1.483566113220077986810e-6,
            5.91213695189949384568e-7, // 20
            1.419864239993567456648e-8,
            -6.35177348504429108355e-8,
            -1.83156167830404631847e-8,
            3.24974651804369739084e-9,
            3.01778054000907084962e-9, // 25
            2.10860063470665179035e-10,
            -3.56323398659765326830e-10,
            -9.05512445092829268740e-11,
            3.47272670930455000726e-11,
            1.771449521401119186147e-11, // 30
            -2.72760231582004518397e-12,
            -2.90768834218286692054e-12,
            1.203145821938798755342e-13,
            4.53296667826067277389e-13,
            1.372562058671550042872e-14, // 35
            -7.07408626028685552231e-14,
            -5.40931028288214223366e-15,
            1.135768719899924165040e-14,
            1.128073562364402060469e-15,
            -1.89969494739492699566e-15 // 40
        };

        /*
        private static double Fadeeva_Weideman_L = Math.Sqrt(20.0 / Math.Sqrt(2.0));
        private static double[] Faddeeva_Weideman_Coefficients = {
            0.145204167876254758882, // 0
            0.136392115479649636711,
            0.112850364752776356308,
            0.081814663814142472504,
            0.051454522270816839073,
            0.027588397192576112169, // 5
            0.012224615682792673110,
            0.004202799708104555100,
            0.00094232196590738341762,
            0.000024250951611669019436,
           -0.000075616164647923403908, // 10
           -0.000025048115578198088771,
            1.35444636058911259445e-6,
            3.1103803514939499680e-6,
            4.3812443485181690079e-7,
           -3.1958626141517655223e-7, // 15
           -9.9512459304812101824e-8,
            3.38225585699344840075e-8,
            1.68263253379413070440e-8,
           -4.206345505308078461584e-9,
           -2.7027917878529318399463e-9 // 20
        };
        */

        // The power series for the error function (DLMF 7.6.1)
        //   erf(z) = \frac{2}{\sqrt{\pi}} \sum_{k=0}^{\infty} \frac{ (-1)^k z^{2k + 1} }{(2k+1) k!}
        //          = \frac{2}{\sqrt{\pi}} \left[ z - z^3 / 3 + z^5 / 10 - \cdots \right]
        // requires about 15 terms at |z| ~ 1, 30 terms at |z| ~ 2, 45 terms at |z| ~ 3

        private static Complex Erf_Series (Complex z) {
            Complex zp = 2.0 / Global.SqrtPI * z;
            Complex zz = - ComplexMath.Sqr(z);
            Complex f = zp;
            for (int k = 1; k < Global.SeriesMax; k++) {
                Complex f_old = f;
                zp *= zz / k;
                f += zp / (2 * k + 1);
                if (f == f_old) return (f);
            }
            throw new NonconvergenceException();
        }

        /// <summary>
        /// Computes the complex error function.
        /// </summary>
        /// <param name="z">The complex argument.</param>
        /// <returns>The value of erf(z).</returns>
        /// <remarks>
        /// <para>This function is the analytic continuation of the error function (<see cref="AdvancedMath.Erf"/>) to the complex plane.</para>
        /// <para>The image below shows the complex error function near the origin, using domain coloring.</para>
        /// <img src="../images/ComplexErfPlot.png" />
        /// <para>The complex error function is entire: it has no poles, cuts, or discontinuities anywhere in the complex plane.</para>
        /// <para>For pure imaginary arguments, erf(z) reduces to the Dawson integral (<see cref="AdvancedMath.Dawson"/>).</para>
        /// <para>Away from the origin near the real axis, the real part of erf(z) quickly approaches &#x0b1;1. To accurately determine
        /// the small difference erf(z) &#8723; 1 in this region, use the <see cref="Faddeeva"/> function. Away from the origin near
        /// the imaginary axis, the magnitude of erf(z) increases very quickly. Although erf(z) may overflow in this region, you
        /// can still accurately determine the value of the product erf(z) exp(z<sup>2</sup>) using the <see cref="Faddeeva"/>
        /// function.</para>
        /// </remarks>
        /// <seealso cref="AdvancedMath.Erf"/>
        /// <seealso cref="AdvancedMath.Dawson"/>
        /// <seealso cref="AdvancedComplexMath.Faddeeva"/>
        public static Complex Erf (Complex z) {

            double r = ComplexMath.Abs(z);

            if (r < 4.0) {
                // near the origin, use the series
                return (Erf_Series(z));
            } else {
                // otherwise, just compute from Faddeva
                if (z.Re < 0.0) {
                    // since Fadddeeva blows up for negative z.Re, use erf(z) = -erf(-z)
                    return (ComplexMath.Exp(-ComplexMath.Sqr(z)) * Faddeeva(-ComplexMath.I * z) - 1.0);
                } else {
                    return (1.0 - ComplexMath.Exp(-ComplexMath.Sqr(z)) * Faddeeva(ComplexMath.I * z));
                }
                // we don't do this near the origin because we would loose accuracy in the very small real parts there by subtracting from 1
            } 

        }

        // Continued fraction expansion (DLMF 7.9.3)
        //   w(z) = \frac{i}{\sqrt{\pi}} \left[ \frac{1}{z -} \frac{1/2}{z -} \frac{1}{z -} \frac{3/2}{z -} \frac{2}{z -} \cdots \right]
        // requires about 10 terms at |z|~10, 15 terms at |z|~7
        // for smaller z, still converges off the real axis, but fails to converge on the real axis for z~6 and below

        private static Complex Faddeeva_ContinuedFraction (Complex z) {
            Complex a = 1.0;			// a_1
            Complex b = z;	            // b_1
            Complex D = 1.0 / b;		// D_1 = b_0/b_1
            Complex Df = a / b;		    // Df_1 = f_1 - f_0
            Complex f = 0.0 + Df;		// f_1 = f_0 + Df_1 = b_0 + Df_1
            for (int k = 1; k < Global.SeriesMax; k++) {
                Complex f_old = f;
                a = -k / 2.0;
                D = 1.0 / (b + a * D);
                Df = (b * D - 1.0) * Df;
                f += Df;
                if (f == f_old) {
                    return (ComplexMath.I / Global.SqrtPI * f);
                }
            }
            throw new NonconvergenceException();
        }

        // given a Faddeeva value w0 at a point z0, this routine computes the Fadeeva value a distance dz away using a Taylor expansion
        // we use this near the real axis, where the Weideman algorithm doesn't give the e^(-x^2) real part of w(x) very accurately,
        // and the Laplace continued fraction misses it entirely

        private static Complex Faddeeva_Taylor (Complex z0, Complex w0, Complex dz) {
            // first order Taylor expansion
            Complex wp_old = w0;
            Complex wp = 2.0 * (ComplexMath.I / Global.SqrtPI - z0 * w0);
            Complex zz = dz;

            Complex w = w0 + wp * dz;
            // higher orders
            for (int k = 2; k < Global.SeriesMax; k++) {
                // remmeber the current value
                Complex w_old = w;
                // compute the next derivative
                Complex wp_new = -2.0 * (z0 * wp + (k - 1) * wp_old);
                wp_old = wp;
                wp = wp_new;
                // use it to generate the next term in the Taylor expansion
                zz = zz * dz / k;

                w = w_old + wp * zz;
                // test whether we have converged
                if (w == w_old) {
                    return (w);
                }
            }
            throw new NonconvergenceException();
        }

    }

}
