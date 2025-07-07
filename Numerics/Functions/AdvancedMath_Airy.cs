using Meta.Numerics.Extended;
using System;
using System.Diagnostics;
using System.Runtime.Serialization;
using System.Text;

namespace Meta.Numerics.Functions
{
    public static partial class AdvancedMath
    {

        internal static readonly double Ai0 = 1.0 / (Math.Pow(3.0, 2.0 / 3.0) * AdvancedMath.Gamma(2.0 / 3.0));

        internal static readonly double AiPrime0 = -1.0 / (Math.Pow(3.0, 1.0 / 3.0) * AdvancedMath.Gamma(1.0 / 3.0));

        internal static readonly double Bi0 = 1.0 / (Math.Pow(3.0, 1.0 / 6.0) * AdvancedMath.Gamma(2.0 / 3.0));

        internal static readonly double BiPrime0 = Math.Pow(3.0, 1.0 / 6.0) / AdvancedMath.Gamma(1.0 / 3.0);

        /// <summary>
        /// Computes the Airy function of the first kind.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value Ai(x).</returns>
        /// <remarks>
        /// <para>For information on the Airy functions, see <see cref="AdvancedMath.Airy"/>.</para>
        /// </remarks>
        /// <seealso cref="AiryBi"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Airy_functions" />
        public static double AiryAi (double x) {
            if (x <= -14.0) {
                return Airy_Asymptotic_Negative(-x).FirstSolutionValue;
            } else if (x < -3.0) {
                double y = 2.0 / 3.0 * Math.Pow(-x, 3.0 / 2.0);
                SolutionPair s = Bessel(1.0 / 3.0, y);
                return 0.5 * Math.Sqrt(-x) * (s.FirstSolutionValue - s.SecondSolutionValue / Global.SqrtThree);
            } else if (x < 2.0) {
                // Close to the origin, use a power series.
                // Convergence is slightly better for negative x than for positive, so we use it further out to the left of the origin.
                return AiryAi_Series(x);
                // The definitions in terms of bessel functions become 0 X infinity at x=0, so we can't use them directly here anyway.
            } else if (x < 14.0) {
                double y = 2.0 / 3.0 * Math.Pow(x, 3.0 / 2.0);
                return Math.Sqrt(x / 3.0) / Math.PI * ModifiedBesselK(1.0 / 3.0, y);
            } else if (x < 108.0) {
                SolutionPair s = Airy_Asymptotic_Positive_Scaled(x, out double xi);
                return Math.Exp(-xi) * s.FirstSolutionValue;
            } else if (x <= Double.PositiveInfinity) {
                // Beyond 108 the function underflows.
                return 0.0;
            } else {
                Debug.Assert(Double.IsNaN(x));
                return x;
            }
        }

        // The Maclaurin series for Ai (DLMF 9.4.1) is:
        //   Ai(x) = Ai(0) [ 1 + 1/3! x^3 + 1 * 4 / 6! x^6 + 1 * 4 * 7 / 9! x^9 + \cdots] +
        //           x Ai'(0) [ 1 + 2 / 4! x^3 + 2 * 5 / 7! x^6 + 2 * 5 * 8 / 10! x^9 + \cdots ]

        private static double AiryAi_Series (double x) {

            // Problematic for positive x once exponential decay sets in; don't use for x > 2

            double p = Ai0;
            double q = AiPrime0 * x;
            double f = p + q;

            double x3 = x * x * x;
            for (int k = 0; k < Global.SeriesMax; k += 3) {
                double f_old = f;
                p *= x3 / ((k + 2) * (k + 3));
                q *= x3 / ((k + 3) * (k + 4));
                f += p + q;
                if (f == f_old) {
                    return f;
                }
            }
            throw new NonconvergenceException();

        }

        /// <summary>
        /// Computes the Airy function of the second kind.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value Bi(x).</returns>
        /// <remarks>
        /// <para>For information on the Airy functions, see <see cref="AdvancedMath.Airy"/>.</para>
        /// <para>While the notation Bi(x) was chosen simply as a natural complement to Ai(x), it has influenced the common
        /// nomenclature for this function, which is now often called the "Bairy function".</para>
        /// </remarks>        
        /// <seealso cref="AiryAi"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Airy_functions" />
        public static double AiryBi (double x) {
            if (x <= -14.0) {
                return Airy_Asymptotic_Negative(-x).SecondSolutionValue;
            } else if (x < -3.0) {
                // Get from Bessel functions.
                double y = 2.0 / 3.0 * Math.Pow(-x, 3.0 / 2.0);
                SolutionPair s = Bessel(1.0 / 3.0, y);
                return -0.5 * Math.Sqrt(-x) * (s.FirstSolutionValue / Global.SqrtThree + s.SecondSolutionValue);
            } else if (x < 5.0) {
                // The Bi series is better than the Ai series for positive values, because it blows up rather than cancelling down.
                // It's also slightly better for negative values, because a given number of oscillations occur further out.  
                return AiryBi_Series(x);
            } else if (x < 14.0) {
                // Get from modified Bessel functions.
                double y = 2.0 / 3.0 * Math.Pow(x, 3.0 / 2.0);
                SolutionPair s = ModifiedBessel(1.0 / 3.0, y);
                return Math.Sqrt(x) * (2.0 / Global.SqrtThree * s.FirstSolutionValue + s.SecondSolutionValue / Math.PI);
            } else if (x < 108.0) {
                SolutionPair s = Airy_Asymptotic_Positive_Scaled(x, out double xi);
                return Math.Exp(xi) * s.SecondSolutionValue;
            } else if (x <= Double.PositiveInfinity) {
                return Double.PositiveInfinity;
            } else {
                Debug.Assert(Double.IsNaN(x));
                return x;
            }
        }

        private static double AiryBi_Series (double x) {

            double p = Bi0;
            double q = x * BiPrime0;
            double f = p + q;

            double x3 = x * x * x;
            for (int k = 0; k < Global.SeriesMax; k += 3) {
                double f_old = f;
                p *= x3 / ((k + 2) * (k + 3));
                q *= x3 / ((k + 3) * (k + 4));
                f += p + q;
                if (f == f_old) return f;
            }
            throw new NonconvergenceException();

        }


        /// <summary>
        /// Computes both Airy functions and their derivatives.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The values of Ai(x), Ai'(x), Bi(x), and Bi'(x).</returns>
        /// <remarks>
        /// <para>Airy functions are solutions to the Airy differential equation:</para>
        /// <img src="../images/AiryODE.png" />
        /// <para>The Airy functions appear in quantum mechanics in the semi-classical WKB solution to the wave functions in a potential.</para>
        /// <para>For negative arguments, Ai(x) and Bi(x) are oscillatory. For positive arguments, Ai(x) decreases exponentially and Bi(x) increases
        /// exponentially with increasing x.</para>
        /// <para>This method simultaneously computes both Airy functions and their derivatives. If you need both Ai and Bi, it is faster to call
        /// this method once than to call <see cref="AiryAi(double)"/> and <see cref="AiryBi(double)"/> separately. If on, the other hand, you need
        /// only Ai or only Bi, and no derivative values, it is faster to call the appropriate method to compute the one you need.</para>
        /// </remarks>
        /// <seealso cref="AiryAi"/>
        /// <seealso cref="AiryBi"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Airy_functions" />
        public static SolutionPair Airy (double x) {
            if (x <= -14.0) {
                return Airy_Asymptotic_Negative(-x);
            } else if (x <= -2.0) {
                // Map to Bessel functions for negative values
                double z = 2.0 / 3.0 * Math.Pow(-x, 3.0 / 2.0);
                SolutionPair p = Bessel(1.0 / 3.0, z);
                double a = (p.FirstSolutionValue - p.SecondSolutionValue / Global.SqrtThree) / 2.0;
                double b = (p.FirstSolutionValue / Global.SqrtThree + p.SecondSolutionValue) / 2.0;
                double sx = Math.Sqrt(-x);
                return new SolutionPair(
                    sx * a, (x * (p.FirstSolutionDerivative - p.SecondSolutionDerivative / Global.SqrtThree) - a / sx) / 2.0,
                    -sx * b, (b / sx - x * (p.FirstSolutionDerivative / Global.SqrtThree + p.SecondSolutionDerivative)) / 2.0
                );
            } else if (x <= 2.0) {
                // Use series near origin
                return Airy_Series(x);
            } else if (x < 14.0) {
                // Map to modified Bessel functions for positive values
                double z = 2.0 / 3.0 * Math.Pow(x, 3.0 / 2.0);
                SolutionPair p = ModifiedBessel(1.0 / 3.0, z);
                double a = 1.0 / (Global.SqrtThree * Math.PI);
                double b = 2.0 / Global.SqrtThree * p.FirstSolutionValue + p.SecondSolutionValue / Math.PI;
                double sx = Math.Sqrt(x);
                return new SolutionPair(
                    a * sx * p.SecondSolutionValue, a * (x * p.SecondSolutionDerivative + p.SecondSolutionValue / sx / 2.0),
                    sx * b, x * (2.0 / Global.SqrtThree * p.FirstSolutionDerivative + p.SecondSolutionDerivative / Math.PI) + b / sx / 2.0
                );
            } else if (x < 108.0) {
                SolutionPair s = Airy_Asymptotic_Positive_Scaled(x, out double xi);
                double e = Math.Exp(xi);
                return new SolutionPair(s.FirstSolutionValue / e, s.FirstSolutionDerivative / e, s.SecondSolutionValue * e, s.SecondSolutionDerivative * e);
            } else if (x <= Double.PositiveInfinity) {
                return new SolutionPair(0.0, -0.0, Double.PositiveInfinity, Double.PositiveInfinity);
            } else {
                Debug.Assert(Double.IsNaN(x));
                return new SolutionPair(x, x, x, x);
            }

            // NR recommends against using the Ai' and Bi' expressions obtained from simply differentiating the Ai and Bi expressions
            // They give instead expressions involving Bessel functions of different orders. But their reason is to avoid the cancellations
            // among terms ~1/x that get large near x ~ 0, and their method would require two Bessel evaluations.
            // Since we use the power series near x ~ 0, we avoid their cancelation problem and can get away with a single Bessel evaluation.

            // It might appear that we can optimize by not computing both \sqrt{x} and x^{3/2} explicitly, but instead computing \sqrt{x} and
            // then (\sqrt{x})^3. But the latter method looses accuracy, presumably because \sqrt{x} compresses range and cubing then expands
            // it, skipping over the double that is actually closest to x^{3/2}. So we compute both explicitly.

        }


        private static SolutionPair Airy_Series (double x) {

            // compute k = 0 terms in f' and g' series, and in f and g series
            // compute terms to get k = 0 terms in a' and b' and a and b series

            double g = -AiPrime0;
            double ap = -g;
            double bp = g;

            double f = Ai0;
            g *= x;
            double a = f - g;
            double b = f + g;

            // we will need to multiply by x^2 to produce higher terms, so remember it
            double x2 = x * x;

            for (int k = 1; k < Global.SeriesMax; k++) {

                // remember old values
                double a_old = a;
                double b_old = b;
                double ap_old = ap;
                double bp_old = bp;

                // compute 3k
                double tk = 3 * k;

                // kth term in f' and g' series, and corresponding a' and b' series
                f *= x2 / (tk - 1);
                g *= x2 / tk;
                ap += (f - g);
                bp += (f + g);

                // kth term in f and g series, and corresponding a and b series
                f *= x / tk;
                g *= x / (tk + 1);
                a += (f - g);
                b += (f + g);

                // check for convergence
                if ((a == a_old) && (b == b_old) && (ap == ap_old) && (bp == bp_old)) {
                    return (new SolutionPair(
                        a, ap,
                        Global.SqrtThree * b, Global.SqrtThree * bp
                    ));
                }
            }

            throw new NonconvergenceException();

        }

        private static SolutionPair Airy_Asymptotic_Positive_Scaled (double x, out double xi) {

            Debug.Assert(x > 0.0);

            xi = 2.0 / 3.0 * Math.Pow(x, 3.0 / 2.0);
        
            // We need the sequences:
            //   s_a = a_0 + a_1 + a_2 + a_3 + \cdots
            //   t_a = a_0 - a_1 + a_2 - a_3 + \cdots
            // and the corresponding series of b-terms.

            double sa = 1.0;
            double ta = 1.0;
            double sb = 1.0;
            double tb = 1.0;

            double a = 1.0;
            for (int k = 1; k <Global.SeriesMax; k++) {

                double sa_old = sa;
                double ta_old = ta;
                double sb_old = sb;
                double tb_old = tb;

                a *= (0.5 * (k - 1) + 5.0 / 72.0 / k) / xi;
                double b = -(1.0 + 2.0 / (6 * k - 1)) * a;

                sa += a;
                sb += b;
                if (k % 2 == 0) {
                    ta += a;
                    tb += b;
                } else {
                    ta -= a;
                    tb -= b;
                }

                if (sa == sa_old && ta == ta_old && sb == sb_old && tb == tb_old) {
                    double q = Math.Pow(x, 1.0 / 4.0);
                    return (new SolutionPair(
                        0.5 / Global.SqrtPI / q * ta,
                        -0.5 / Global.SqrtPI * q * tb,
                        1.0 / Global.SqrtPI / q * sa,
                        1.0 / Global.SqrtPI * q * sb
                    ));
                } 
            }

            throw new NonconvergenceException();

        }

        private static SolutionPair Airy_Asymptotic_Negative (double x) {

            Debug.Assert(x > 0.0);

            double xi = 2.0 / 3.0 * Math.Pow(Math.Abs(x), 3.0 / 2.0);

            // We need the sequences:
            //   u_a = a_0 + a_1 - a_2 - a_3 + \cdots
            //   v_a = a_0 - a_1 - a_2 + a_3 + \cdots
            // and the corresponding series of b-terms.

            double ua = 1.0;
            double va = 1.0;
            double ub = 1.0;
            double vb = 1.0;

            double a = 1.0;
            for (int k = 1; k < Global.SeriesMax; k++) {

                double ua_old = ua;
                double va_old = va;
                double ub_old = ub;
                double vb_old = vb;

                a *= (0.5 * (k - 1) + 5.0 / 72.0 / k) / xi;
                double b = -(1.0 + 2.0 / (6 * k - 1)) * a;

                switch (k % 4) {
                    case 0:
                        ua += a; va += a; ub += b; vb += b;
                        break;
                    case 1:
                        ua += a; va -= a; ub += b; vb -= b;
                        break;
                    case 2:
                        ua -= a; va -= a; ub -= b; vb -= b;
                        break;
                    case 3:
                        ua -= a; va += a; ub -= b; vb += b;
                        break;
                    default:
                        throw new NotImplementedException();
                }

                if (ua == ua_old && va == va_old && ub == ub_old && vb == vb_old) {
                    double c = MoreMath.Cos(xi);
                    double s = MoreMath.Sin(xi);
                    double q = Math.Pow(x, 1.0 / 4.0);
                    double n = 1.0 / Global.SqrtTwoPI / q;
                    double m = 1.0 / Global.SqrtTwoPI * q;
                    return new SolutionPair(
                        n * (s * ua + c * va), m * (s * vb - c * ub),
                        n * (c * ua - s * va), m * (s * ub + c * vb)
                    );
                }

            }

            throw new NonconvergenceException();

        }

        /// <summary>
        /// Computes the requested zero of the Airy Ai function.
        /// </summary>
        /// <param name="k">The index of the zero.</param>
        /// <returns>The <paramref name="k"/>th value of x for which Ai(x) = 0.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> is less than 1.</exception>
        /// <seealso cref="AiryAi(double)"/>
        public static double AiryAiZero (int k) {
            return (BesselMath.AiryAiZero(k));
        }

        /// <summary>
        /// Computes the requested zero of the Airy Bi (Bairy) function.
        /// </summary>
        /// <param name="k">The index of the zero.</param>
        /// <returns>The <paramref name="k"/>th value of x for which Bi(x) = 0.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> is less than 1.</exception>
        /// <seealso cref="AiryBi(double)"/>
        public static double AiryBiZero (int k) {
            return (BesselMath.AiryBiZero(k));
        }

    }


    public static partial class AdvancedComplexMath {

        /// <summary>
        /// Computes the complex Airy function of the first kind.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The value of Ai(z).</returns>
        /// <remarks>
        /// <para>The conplex Airy function is entire, meaning it has no poles or cuts anywhere in the complex plane.</para>
        /// </remarks>
        public static Complex AiryAi (Complex z) {

            if (-3.0 < z.Re && z.Re < 2.0 && Math.Abs(z.Im) < 4.0 ) {
                // Sufficiently near the origin, use the power series.
                return AiryAi_Series(z);
            } else if (ComplexMath.Abs(z) < 12.0) {
                // In the intermediate region, use integration.
                return AiryAi_Integration(z);
            } else {
                // Sufficiently far from the origin, use the asymototic series.
                return AiryAi_Asymptotic(z);
            }

        }

        /// <summary>
        /// Computes the complex Airy function of the second kind.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The value of Bi(z).</returns>
        /// <remarks>
        /// <para>The conplex Bairy function is entire, meaning it has no poles or cuts anywhere in the complex plane.</para>
        /// </remarks>
        public static Complex AiryBi (Complex z) {
            if (-3.0 < z.Re && z.Re < 5.0 && Math.Abs(z.Im) < 4.0) {
                // Sufficiently near the origin, use the power series.
                return AiryBi_Series(z);
            } else if (ComplexMath.Abs(z) < 12.0) {
                // In the intermediate region, fall back to the relation to Ai.
                return r12 * AiryAi(r3 * z) + r12.Conjugate * AiryAi(r3.Conjugate * z);
            } else {
                // Sufficiently far from the origin, use the asymptotic series.
                return AiryBi_Asymptotic(z);
            }
        }


        // Gil, Segura, and Temme, "Computing complex Airy functions by numerical quadrature"
        // Numerical Algorithms 30 (2002) 11-23

        private static Complex AiryAi_Integration(Complex z) {

            // Our integral only applies for |arg(z)| < 2/3 \pi, so if we are in the wrong 1/3rd of the complex
            // plane, use this relation to get us out of there.
            if (Math.Abs(z.Im) < Global.SqrtThree * -z.Re) {
                return -(r3 * AiryAi_Integration(r3 * z) + r3.Conjugate * AiryAi_Integration(r3.Conjugate * z));
            }
            Debug.Assert(Math.Abs(ComplexMath.Arg(z)) <= 2.0 / 3.0 * Math.PI);

            Complex xi = 2.0 / 3.0 * ComplexMath.Pow(z, 3.0 / 2.0);

            Complex s;
            if (z.Re >= 0.0) {
                // DLMF 9.5.8
                // Ai(z) = \int_0^{\infty} \! \left( 2 + \frac{t}{\xi} \right)^{-1/6} t^{-1/6} e^{-t} \, dt 
                Complex one_over_xi = 1.0 / xi;
                s = GaussLaguerre(t => ComplexMath.Pow(2.0 + one_over_xi * t, -1.0 / 6.0));
            } else {
                // The expression (2 + \frac{t}{\xi}) has a zero for real t, and therefore the usual integrand has a pole,
                // for \theta = arg(z) = \pm 2/3 \pi. For smaller |z|, this pole occurs at smaller t, less supressed by e^{-t}.
                // Even for nearby \theta, for which the pole becomes a peak, the integrand is still highly non-polynomial
                // so Gauss-Laguerre has large errors. For z = -1 + \sqrt{3} i, for example, we get only 4 accurate
                // digits of Ai(z).
                // Gil et al suggest the change of integration variable t -> t (1 + i \tan \tau), which essentially
                // means integrating along a different ray in the complex plane from zero to infinity, which avoids
                // the singularity. The price, as you discover after the change of variable, it an osciallting
                // factor e^{-i t \tan \tau}, so you don't want \tan \tau to get too big. Gil et all suggest
                // \tau = 3/2 (\theta - \pi / 2), which means \theta = \pi / 2 -> \tau = 0 and \theta = 3 \pi / 2 ->
                // \tau = \pi / 4. This is nice because (i) |\tan \tau| < 1 and (ii) at \theta = \pi / 2 it agrees
                // with the old expression, so the imaginary axis becomes a natural transition point.
                // There may be better changes of variable, but I haven't tried to optimize; their suggestion works.
                double theta = ComplexMath.Arg(z);
                double tau = Math.Sign(theta) * 3.0 * (Math.Abs(theta) - Math.PI / 2.0) / 2.0;
                double tanTau = Math.Tan(tau);
                Complex r = new Complex(1.0, tanTau);
                Complex one_over_xiPrime = r / xi;
                s = GaussLaguerre(t => ComplexMath.Pow(2.0 + one_over_xiPrime * t, -1.0 / 6.0) * ComplexMath.Exp(-tanTau * t * Complex.I));
                s *= ComplexMath.Pow(r, 5.0 / 6.0);
            }
            s /= Global.SqrtPI * AdvancedMath.Gamma(5.0 / 6.0) * ComplexMath.Exp(xi) * ComplexMath.Pow(48.0 * xi, 1.0 / 6.0);

            return s;

        }

        private static readonly Complex r3 = new Complex(-0.5, 0.5 * Global.SqrtThree); // e^{2 \pi i / 3} = e^{i \pi 2/3}

        private static readonly Complex r6 = new Complex(0.5, 0.5 * Global.SqrtThree); // e^{2 \pi i / 6} = e^{i \pi / 3}

        private static readonly Complex r12 = new Complex(0.5 * Global.SqrtThree, 0.5); // e^{2 \pi i / 12} = e^{i \pi / 6}

        private static Complex GaussLaguerre (Func<double, Complex> f) {
            Complex s = Complex.Zero;
            Debug.Assert(gaussLaguerreAbcissas.Length == gaussLaguerreWeights.Length);
            for (int i = 0; i < gaussLaguerreAbcissas.Length; i++) {
                double t = gaussLaguerreAbcissas[i];
                Complex u = f(t);
                s += gaussLaguerreWeights[i] * u;
            }
            return s;
        }

        // I got these Gauss-Laguerre absissias and weights for \alpha = -1/6 from https://keisan.casio.com/exec/system/1281279441

        private static readonly double[] gaussLaguerreAbcissas = new double[] {
            0.035395864628857473341,
            0.213221067973112801407,
            0.5437031179102081500529,
            1.027709813466397791865,
            1.66640897470266755049,
            2.46135161401298910305,
            3.41449892069320827029,
            4.528250905488224631257,
            5.805481666804404269353,
            7.24958316010105625397,
            8.86451930505549785447,
            10.65489271569047635826,
            12.62602704095225320732,
            14.7840688844501833364,
            17.13611463893863901777,
            19.69036950114667397932,
            22.45634870646121706382,
            25.44513509005754032126,
            28.66971317766305134724,
            32.14540937639244886562,
            35.89048263433744709797,
            39.92693406675568300285,
            44.28164484760133505655,
            48.98802364101213113421,
            54.08847832852900608478,
            59.63828965538021818053,
            65.7120216743068648735,
            72.41490018813192702663,
            79.90499555543955427487,
            88.44262309973055869568,
            98.5256453490892114923,
            111.4344240837636221002
        };

        private static readonly double[] gaussLaguerreWeights = new double[] {
            0.171523294924177243413,
            0.265541055251247591941,
            0.261612743720697231084,
            0.199873080038620509854,
            0.124324866435969689392,
            0.064151475171879700905,
            0.027689764505614639628,
            0.01003323933072850002218,
            0.00305465020190647387487,
            7.8075236437140582058E-4,
            1.67171402008372942307E-4,
            2.9887718568438169971E-5,
            4.4426681724328769955E-6,
            5.46141624850188518657E-7,
            5.5168591512386702208E-8,
            4.544584666340893105619E-9,
            3.02555613613217453064E-10,
            1.61071992178314859744E-11,
            6.7717426155550971066E-13,
            2.21509883100294082875E-14,
            5.538383156152087434814E-16,
            1.03598569005705957861E-17,
            1.41224276995732571164E-19,
            1.3579426050820722803E-21,
            8.8367153703821723821E-24,
            3.6874306944681702237E-26,
            9.1759547194454262391E-29,
            1.22938512559029336904E-31,
            7.6101892480080210749E-35,
            1.6940913721216430991E-38,
            8.4317439269185431821E-43,
            2.8203249181741957545E-48
        };

        private static Complex AiryAi_Series (Complex z) {
            // DLMF 9.4.1
            // |z|~4 needs k~50, looses ~2 digits of accuracy
            // |z|~8 needs k~90, looses ~3 digits of accuracy
            // Does better for z < 0 than z > 0, because the latter is trying to create an exponential fall-off via cancellation among terms.
            return Airy_Series(AdvancedMath.Ai0, AdvancedMath.AiPrime0, z);
        }

        private static Complex AiryBi_Series(Complex z) {
            // DLMF 9.4.3
            return Airy_Series(AdvancedMath.Bi0, AdvancedMath.BiPrime0, z);
        }

        // DLMF 9.4.1 and 9.4.3
        // Ai and Bi share the same series structure; they only have diffierent coefficients of the two sub-series.
        // So re-use the same code for both.

        private static Complex Airy_Series (double a, double b, Complex z) {
            Complex p = a;
            Complex q = b * z;
            Complex f = p + q;
            Complex z3 = ComplexMath.Sqr(z) * z;
            for (int k = 0; k < Global.SeriesMax; k += 3) {
                Complex f_old = f;
                p *= z3 / ((k + 2) * (k + 3));
                q *= z3 / ((k + 3) * (k + 4));
                f += p + q;
                if (f == f_old) return f;
            }
            throw new NonconvergenceException();
        }

        private static Complex AiryAi_Asymptotic (Complex z) {

            // The "left" series is mathematically valid for |\arg(z)| < \pi, i.e. right up until the negative real axis.
            // The "right" series is mathematically valid for |\arg(-z)| < 2/3 \pi, i.e. 1/6 of the way into the positive real quadrants.
            // Together they more than cover the complex plane. Pick an easy to evaluate cross-over equally far from either limit.

            if (Math.Abs(ComplexMath.Arg(z)) < 2.0 / 3.0 * Math.PI) {
                // DLMF 9.7.5
                Complex xi = 2.0 / 3.0 * ComplexMath.Pow(z, 3.0 / 2.0);
                Airy_Asymptotic_Subseries(xi, out _, out Complex u1, out _, out _);
                Complex e = ComplexMath.Exp(xi);
                Complex q = ComplexMath.Pow(z, 1.0 / 4.0);
                return 0.5 / Global.SqrtPI / q / e * u1;
            } else {
                // DLMF 9.7.9
                z = -z;
                Complex xi = 2.0 / 3.0 * ComplexMath.Pow(z, 3.0 / 2.0);
                Airy_Asymptotic_Subseries(xi, out _, out _, out Complex u12, out Complex u23);
                Complex c = ComplexMath.Cos(xi);
                Complex s = ComplexMath.Sin(xi);
                Complex q = ComplexMath.Pow(z, 1.0 / 4.0);
                return 1.0 / Global.SqrtPI / Global.SqrtTwo / q * (c * u12 + s * u23);
            }
     
        }


        private static Complex AiryBi_Asymptotic_Positive (Complex z) {
            Complex xi = 2.0 / 3.0 * ComplexMath.Pow(z, 3.0 / 2.0);
            Airy_Asymptotic_Subseries(xi, out Complex u0, out _, out _, out _);
            Complex e = ComplexMath.Exp(xi);
            Complex q = ComplexMath.Pow(z, 1.0 / 4.0);
            return 1.0 / Global.SqrtPI / q * e * u0;
        }

        private static Complex AiryBi_Asymptotic_Negative (Complex z) {
            z = -z;
            Complex xi = 2.0 / 3.0 * ComplexMath.Pow(z, 3.0 / 2.0);
            Airy_Asymptotic_Subseries(xi, out _, out _, out Complex u12, out Complex u23);
            Complex c = ComplexMath.Cos(xi);
            Complex s = ComplexMath.Sin(xi);
            Complex q = ComplexMath.Pow(z, 1.0 / 4.0);
            return 1.0 / Global.SqrtPI / Global.SqrtTwo / q * (c * u23 - s * u12);
        }

        private static Complex AiryBi_Asymptotic_Third (Complex z) {
            // Rather than rotating z by \pi / 3 and then computing xi, etc., we preserve
            // symmetries and accuracy better by expliciting applying the known effect of the rotation.
            Complex xi = 2.0 / 3.0 * ComplexMath.Pow(z, 3.0 / 2.0);
            Complex d = new Complex(0.0, 0.5 * Global.LogTwo);
            Complex r = r12;
            if (z.Im < 0.0) {
                xi *= Complex.I;
                z *= r6;
                d = -d;
                r = r.Conjugate;
            } else {
                xi *= -Complex.I;
                z *= r6.Conjugate;
            }
            Airy_Asymptotic_Subseries(xi, out _, out _, out Complex u12, out Complex u23);
            Complex c = ComplexMath.Cos(xi - d);
            Complex s = ComplexMath.Sin(xi - d);
            Complex q = ComplexMath.Pow(z, 1.0 / 4.0);
            return 1.0 / Global.SqrtPI * (r / q) * (c * u12 + s * u23);
        }

        private static Complex AiryBi_Asymptotic (Complex z) {
            // The Bi asymptotic series have smaller ranges of validity than the Ai asymptotic series.
            // \pm 1/3 \pi around the positive real axis, \pm 2/3 \pi around the negative real axis, and
            // \pm 2/3 \pi around the rays \arg(z) = \pm 1/3 \pi that divide the first two regions.
            // We need all three in order not to require evaluating a series at the limit of its validity.
           if (z.Re < 0.0) {
                return AiryBi_Asymptotic_Negative(z);
            } else if (Math.Abs(z.Im) < z.Re / Global.SqrtThree) {
                return AiryBi_Asymptotic_Positive(z);
            } else {
                return AiryBi_Asymptotic_Third(z);
            }
        }

        private static void Airy_Asymptotic_Subseries (Complex xi, out Complex u0, out Complex u1, out Complex u12, out Complex u23) {

            // Convergence to double accuracy fails for |z| < 9
            // For |z| ~ 9 needs ~28 terms, looses 2d of accuracy
            // For |z| ~ 10 needs ~24 terms, looses 1d of accuracy

            Complex x = (1.0 / 216.0) / xi;
            Complex du = Complex.One;

            u0 = Complex.One;
            u1 = Complex.One;
            u12 = Complex.One;
            u23 = Complex.One;

            for (int k = 1; k < Global.SeriesMax; k++) {

                Complex u0_old = u0;
                Complex u1_old = u1;
                Complex u12_old = u12;

                int six_k = 6 * k;
                du *= x * ((six_k - 5) * (six_k - 3) * (six_k - 1)) / (k * (2 * k - 1));

                u0 += du;
                switch (k % 4) {
                    case 0:
                        u1 += du;
                        u12 += du;
                        u23 += du;
                        break;
                    case 1:
                        u1 -= du;
                        u12 -= du;
                        u23 += du;
                        break;
                    case 2:
                        u1 += du;
                        u12 -= du;
                        u23 -= du;
                        break;
                    case 3:
                        u1 -= du;
                        u12 += du;
                        u23 -= du;
                        break;
                }

                if (u0_old == u0 && u1_old == u1 && u12_old == u12) {
                    return;
                }

            }

            throw new NonconvergenceException();

        }

    }


}
