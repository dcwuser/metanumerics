using System;
using System.Diagnostics;


namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {


        /// <summary>
        /// Computes the error function.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of erf(<paramref name="x"/>).</returns>
        /// <remarks>
        /// <para>The error can be defined via a Gaussian integral.</para>
        /// <img src="../images/ErfIntegral.png" />
        /// <para>The area under a bell curve (<see cref="Meta.Numerics.Statistics.Distributions.NormalDistribution"/>) within &#x2213;z
        /// standard deviations of the mean is given by erf(z/&#x221A;2).</para>
        /// <para>For large values of x, erf(x) &#x2248; 1 to within floating-point accuracy. To obtain accurate values of erfc(x) = 1 - erf(x)
        /// in this range, use the <see cref="Erfc" /> function.</para>
        /// <para>The error function for complex arguments can be computed using <see cref="AdvancedComplexMath.Faddeeva"/>.</para>
        /// </remarks>
        /// <seealso cref="Erfc"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Error_function" />
        /// <seealso href="http://mathworld.wolfram.com/Erf.html" />
        public static double Erf (double x) {
            if (x < 0.0) {
                return (-LeftRegularizedGamma(0.5, x * x));
            } else {
                return (LeftRegularizedGamma(0.5, x * x));
            }
        }

        /// <summary>
        /// Computes the complementary error function.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of erfc(<paramref name="x"/>) = 1 - erf(<paramref name="x"/>).</returns>
        /// <remarks>
        /// <para>The complementary error function can be used to express the area in the tails of a Bell curve beyond a given distance from its center.</para>
        /// <para>It can be defined via an integral as erfc(x) = (2/&#x221A;&#x3C0;) <sub>x</sub>&#x222B;<sup>&#x221E;</sup>dt e<sup>-t<sup>2</sup></sup>.</para>
        /// <para>For small values of x, erfc(x) &#x2248; 1 to within floating-point accuracy. To obtain accurate values of erfc(x) = 1 - erf(x)
        /// in this region, use the <see cref="Erf" /> function.</para></remarks>
        /// <seealso cref="Erf" />
        public static double Erfc (double x) {
            if (x < 0.0) {
                return (1.0 + LeftRegularizedGamma(0.5, x * x));
            } else {
                return (RightRegularizedGamma(0.5, x * x));
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
        public static double InverseErf (double y) {
            if (Math.Abs(y) > 1.0) throw new ArgumentOutOfRangeException("y");
            if (y < 0.0) {
                return (-InverseErf(-y));
            } else {

                if (y == 1.0) return (Double.PositiveInfinity);

                // find y such that Erf(y) = x

                // make an initial guess; these get us there within about 10%
                double x;
                if (y < 0.8) {
                    // for small values, use series inversion
                    // the three-term series is good to within 5% for y~0.8, even better for smaller y
                    //double yy = y * y;
                    //x = (Global.SqrtPI / 2.0) * y * (1.0 + (Math.PI / 12.0) * yy * (1.0 + (7.0 * Math.PI / 40.0) * yy));
                    x = ApproximateInverseErf(y);
                } else {
                    // for large values, use a crude approximation
                    // this is pretty consistently 10% too high; how to correct?
                    //x = Math.Sqrt(-Math.Log((1.0 - y) * (1.0+ y)));

                    // this approximation is good to within 10% for y~0.8, 4% y~0.9.
                    // can we get the next term?
                    //double log = -Math.Log(Global.SqrtPI * (1.0 - y));
                    //x = Math.Sqrt(log) * (1.0 - Math.Log(log) / log / 4.0);

                    //x = Math.Sqrt(-Math.Log(y * Math.Sqrt(-Math.PI * Math.Log(y))));
                    //x = Math.Sqrt(Math.Log(Math.Sqrt(Math.PI)*(1.0-y)
                    x = ApproximateInverseErfc(1.0 - y);
                }

                // refine it via 2nd order Newton's method (aka Halley's method)
                // this typically takes 2-4 cycles to converge to full precision
                double dx = Double.MaxValue;
                for (int i = 1; i < 25; i++) {
                    double x_old = x;
                    double dx_old = dx;

                    double z = Erf(x) - y;
                    dx = -z / (2.0 / Global.SqrtPI * Math.Exp(-x * x) + x * z);

                    //Debug.WriteLine(String.Format("x={0}, z={1}, dx={2}", x, z, dx));
                    x += dx;

                    // return if we have converged
                    if (x == x_old) return (x);

                    // return if we are thrashing
                    if (Math.Abs(dx) >= Math.Abs(dx_old)) return (x);
                }

                throw new NonconvergenceException();

            }
        }


        private static double ApproximateInverseErf (double y) {
            double x = Global.SqrtPI * y / 2.0;
            double xx = x * x;
            double S = 1.0 + xx / 3.0 + 7.0 / 30.0 * xx * xx + 127.0 / 630.0 * xx * xx * xx;
            return (x * S);

        }

        private static double ApproximateInverseErfc (double y) {
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
            } else if (x > 10.0) {
                // use the asymptotic expansion for large values
                return (Dawson_Asymptotic(x));
            } else {
                // use the Rybicki algorithm in between
                return (Dawson_Rybicki(x));
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

        // a rather strange algorithm by Rybicki, described in Numerical Recipies
        // it expresses Dawson's integral as a sum over exponentially supressed terms with a suppression factor h
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
                    //Console.WriteLine(k);
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
            if (x < 0.0) return (-FresnelC(-x));
            if (x < 2.0) {
                return (FresnelC_Series(x));
            } if (x > 64.0) {
                Complex f = Fresnel_Asymptotic(x);
                return (f.Re);
            } else {
                Complex f = Fresnel_ContinuedFraction(x);
                return (f.Re);
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
            if (x < 0.0) return (-FresnelS(-x));
            if (x < 2.0) {
                return (FresnelS_Series(x));
            } if (x > 64.0) {
                Complex f = Fresnel_Asymptotic(x);
                return (f.Im);
            } else {
                Complex f = Fresnel_ContinuedFraction(x);
                return (f.Im);
            }
        }

        /// <summary>
        /// Computes the Fresnel cosine and sine integrals.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value C(x) + i S(x).</returns>
        /// <remarks>
        /// <para>A plot of all values of this function in the complex plane as x ranges from
        /// negative infinity to positive infinity is called a Cornu spiral.</para>
        /// <para>The Fresnel function can be related to the complex error function along the line (1-I).</para>
        /// <img src="../images/FresnelErfRelation.png" />
        /// </remarks>
        /// <seealso cref="FresnelS"/>
        /// <seealso cref="FresnelC"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Fresnel_integral"/>
        /// <seealso href="http://mathworld.wolfram.com/CornuSpiral.html"/>
        public static Complex Fresnel (double x) {
            if (x < 0.0) return (-Fresnel(-x));
            if (x < 2.0) {
                return (new Complex(FresnelC_Series(x), FresnelS_Series(x)));
            } else if (x > 64.0) {
                return (Fresnel_Asymptotic(x));
            } else {
                return (Fresnel_ContinuedFraction(x));
            }
        }
        
        // the use of the asymptotic series doesn't appear to be strictly necessary; the
        // continued fraction converges quickly to the right result in just a few terms
        // for arbitrarily large x; ahha! this is because the continued fraction is just the
        // asymptotic series in that limit, down to the evaluation of sin(pi x^2/2)

        // series requires ~10 terms at x~1, ~20 terms at x~2

        private static double FresnelC_Series (double x) {

            // compute the zero-order value
            double df = x;

            // pre-compute the factor the series is in
            double x4 = x * x * Math.PI / 2.0;
            x4 = x4 * x4;

            // add corrections as needed
            double f = df;
            for (int n = 1; n < Global.SeriesMax; n++) {
                double f_old = f;
                df = - df * x4 / (2 * n) / (2 * n - 1);
                f = f + df / (4 * n + 1);
                if (f == f_old) return (f);
            }

            throw new NonconvergenceException();
        }

        private static double FresnelS_Series (double x) {

            // compute the zero-order value
            double df = x * x * x * Math.PI / 2.0;

            // pre-compute the factor the series is in
            double x4 = x * x * Math.PI / 2.0;
            x4 = x4 * x4;

            // add corrections as needed
            double f = df / 3.0;
            for (int n = 1; n < Global.SeriesMax; n++) {
                double f_old = f;
                df = -df * x4 / (2 * n + 1) / (2 * n);
                f = f + df / (4 * n + 3);
                if (f == f_old) return (f);
            }

            throw new NonconvergenceException();
        }

        // ~3 terms at x~50, ~6 terms at x~10, ~10 terms at x~7,
        // fails to converge to double precision much below that

        private static Complex Fresnel_Asymptotic (double x) {

            double x2 = x * x * Math.PI;
            double x4 = x2 * x2;

            double f = 1.0;
            double g = 1.0;

            double d = 1.0;
            int n = 1;
            while (true) {

                double f_old = f;
                d = -d * (4 * n - 1) / x4;
                f += d;
                double g_old = g;
                d = d * (4 * n + 1);
                g += d;

                if ((f == f_old) && (g == g_old)) break;

                n++;
                if (n > Global.SeriesMax) throw new NonconvergenceException();
            }

            double px = Math.PI * x;
            f = f / px;
            g = g / x2 / px;

            double xx = x * x / 4.0;
            double sin = Sin(0.0, xx);
            double cos = Cos(0.0, xx);

            double C = 0.5 + f * sin - g * cos;
            double S = 0.5 - f * cos - g * sin;

            return (new Complex(C, S));

        }

        // ~35 terms at x~2, ~12 terms at x~4, ~6 terms at x~8, ~3 terms for x~12 and above

        private static Complex Fresnel_ContinuedFraction (double x) {

            double px2 = Math.PI * x * x;

            Complex z = Global.SqrtPI * x / 2.0 * (1.0 - ComplexMath.I);
            //Console.WriteLine("z = {0}", z);
            //Console.WriteLine("z^2 = {0}", z * z);

            // investigate this carefully: it appears that (1) imaginary part of f doesn't change and
            // (2) real part of f merely increases at constant rate until "right" number is reached

            Complex a = 1.0;			            // a_1
            Complex b = new Complex(1.0, -px2);     // b_1
            Complex D = 1.0 / b;                    // D_1 = b_0 / b_1
            Complex Df = a / b;		                // Df_1 = f_1 - f_0
            Complex f = 0.0 + Df;		            // f_1 = f_0 + Df_1 = b_0 + Df_1
            for (int k = 1; true ; k++) {
                //Console.WriteLine("f={0}", b);
                Complex f_old = f;
                a = -(2 * k - 1) * (2 * k);
                b = b + 4.0;
                D = 1.0 / (b + a * D);
                Df = (b * D - 1.0) * Df;
                f += Df;
                if (f == f_old) break;
                if (k > Global.SeriesMax) throw new NonconvergenceException();
            }
            //Console.WriteLine(f);
            double xx = x * x / 4.0;
            Complex j = new Complex(1.0, -1.0);
            Complex e = new Complex(Cos(0.0, xx), Sin(0.0, xx));
            Complex erfc = j * e * f * x;
            //Console.WriteLine("erfc={0}", erfc);
            Complex erf = 1.0 - erfc;
            //Console.WriteLine("erf={0}", erf);
            return (erf * new Complex(1.0, 1.0) / 2);   
            
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
        /// <para>For purely imaginary values, it can be reduced to the error function. For purely real values, it can be reduced to Dawson's integral.</para>
        /// <para>It appears in the computation of the Voigt line profile function V(x;&#x3C3;,&#x3B3;).</para>
        /// </remarks>
        /// <seealso cref="AdvancedMath.Erf" />
        /// <seealso cref="AdvancedMath.Erfc" />
        /// <seealso cref="AdvancedMath.Dawson"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Voigt_profile" />
        public static Complex Faddeeva (Complex z) {

            // use reflection formulae to ensure that we are in the first quadrant
            if (z.Im < 0.0) return (2.0 * ComplexMath.Exp(-z * z) - Faddeeva(-z));
            if (z.Re < 0.0) return (Faddeeva(-z.Conjugate).Conjugate);

            double r = ComplexMath.Abs(z);
            if (r < 2.0) {
                // use series for small z
                return (Faddeeva_Series(z));
            } else if ((z.Im < 0.1) && (z.Re < 30.0)) {
                // this is a special, awkward region
                // along the real axis, Re{w(x)} ~ e^{-x^2}; the Weideman algorthm doesen't compute this small number
                // well and the Laplace continued fraction misses it entirely; therefore very close to the real axis
                // we will use an analytic result on the real axis and Taylor expand to where we need to go.
                // unfortunately the Taylor expansion converges poorly for large x, so we drop this work-arround near x~30,
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

        // requires about 15 terms at |z|~1, 30 terms at |z|~2, 45 terms at |z|~3
        // converges at about the same rate along both real and imaginary axis

        private static Complex Faddeeva_Series (Complex z) {
            Complex df = 2.0 / Global.SqrtPI * ComplexMath.I * z;
            Complex f = 1.0 + df;
            Complex zz = z * z;
            for (int i = 1; i < 250; i++) {
                Complex f_old = f;
                df = df * zz / i;
                f += df / (2 * i + 1);
                if (f == f_old) {
                    return (ComplexMath.Exp(-zz) * f);
                }
            }
            throw new NonconvergenceException();
        }

        // requires about 10 terms at |z|~10, 15 terms at |z|~7
        // for smaller z, still converges off the real axis, but fails to converge on the real axis for z~6 and below

        private static Complex Faddeeva_ContinuedFraction (Complex z) {
            Complex a = 1.0;			// a_1
            Complex b = z;	// b_1
            Complex D = 1.0 / b;		// D_1 = b_0/b_1
            Complex Df = a / b;		// Df_1 = f_1 - f_0
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
