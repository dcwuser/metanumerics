using System;
using System.Diagnostics;
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
                return (Airy_Asymptotic_Negative(-x).FirstSolutionValue);
            } else if (x < -3.0) {
                double y = 2.0 / 3.0 * Math.Pow(-x, 3.0 / 2.0);
                SolutionPair s = Bessel(1.0 / 3.0, y);
                return (Math.Sqrt(-x) / 2.0 * (s.FirstSolutionValue - s.SecondSolutionValue / Global.SqrtThree));
            } else if (x < 2.0) {
                // Close to the origin, use a power series.
                // Convergence is slightly better for negative x than for positive, so we use it further out to the left of the origin.
                return (AiryAi_Series(x));
                // The definitions in terms of bessel functions become 0 X infinity at x=0, so we can't use them directly here anyway.
            } else if (x < 14.0) {
                double y = 2.0 / 3.0 * Math.Pow(x, 3.0 / 2.0);
                return (Math.Sqrt(x / 3.0) / Math.PI * ModifiedBesselK(1.0 / 3.0, y));
            } else if (x < 108.0) {
                SolutionPair s = Airy_Asymptotic_Positive_Scaled(x, out double xi);
                return (Math.Exp(-xi) * s.FirstSolutionValue);
            } else if (x <= Double.PositiveInfinity) {
                return (0.0);
            } else {
                Debug.Assert(Double.IsNaN(x));
                return (x);
            }
        }

        // The Maclaurin series for Ai (DLMF 9.4.1) is:
        //   Ai(x) = Ai(0) [ 1 + 1/3! x^3 + 1 * 4 / 6! x^6 + 1 * 4 * 7 / 9! x^9 + \cdots] +
        //           x Ai'(0) [ 1 + 2 / 4! x^3 + 2 * 5 / 7! x^6 + 2 * 5 * 8 / 10! x^9 + \cdots ]

        private static double AiryAi_Series (double x) {

            // problematic for positive x once exponential decay sets in; don't use for x > 2
            Debug.Assert(x <= 2.0);

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
                    return (f);
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
                return (Airy_Asymptotic_Negative(-x).SecondSolutionValue);
            } else if (x < -3.0) {
                // Get from Bessel functions.
                double y = 2.0 / 3.0 * Math.Pow(-x, 3.0 / 2.0);
                SolutionPair s = Bessel(1.0 / 3.0, y);
                return (-Math.Sqrt(-x) / 2.0 * (s.FirstSolutionValue / Global.SqrtThree + s.SecondSolutionValue));
            } else if (x < 5.0) {
                // The Bi series is better than the Ai series for positive values, because it blows up rather than cancelling down.
                // It's also slightly better for negative values, because a given number of oscillations occur further out.  
                return (AiryBi_Series(x));
            } else if (x < 14.0) {
                // Get from modified Bessel functions.
                double y = 2.0 / 3.0 * Math.Pow(x, 3.0 / 2.0);
                SolutionPair s = ModifiedBessel(1.0 / 3.0, y);
                return (Math.Sqrt(x) * (2.0 / Global.SqrtThree * s.FirstSolutionValue + s.SecondSolutionValue / Math.PI));
            } else if (x < 108.0) {
                SolutionPair s = Airy_Asymptotic_Positive_Scaled(x, out double xi);
                return (Math.Exp(xi) * s.SecondSolutionValue);
            } else if (x <= Double.PositiveInfinity) {
                return (Double.PositiveInfinity);
            } else {
                Debug.Assert(Double.IsNaN(x));
                return (x);
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
                if (f == f_old) return (f);
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
                return (Airy_Asymptotic_Negative(-x));
            } else if (x <= -2.0) {
                // Map to Bessel functions for negative values
                double z = 2.0 / 3.0 * Math.Pow(-x, 3.0 / 2.0);
                SolutionPair p = Bessel(1.0 / 3.0, z);
                double a = (p.FirstSolutionValue - p.SecondSolutionValue / Global.SqrtThree) / 2.0;
                double b = (p.FirstSolutionValue / Global.SqrtThree + p.SecondSolutionValue) / 2.0;
                double sx = Math.Sqrt(-x);
                return (new SolutionPair(
                    sx * a, (x * (p.FirstSolutionDerivative - p.SecondSolutionDerivative / Global.SqrtThree) - a / sx) / 2.0,
                    -sx * b, (b / sx - x * (p.FirstSolutionDerivative / Global.SqrtThree + p.SecondSolutionDerivative)) / 2.0
                ));
            } else if (x <= 2.0) {
                // Use series near origin
                return (Airy_Series(x));
            } else if (x < 14.0) {
                // Map to modified Bessel functions for positive values
                double z = 2.0 / 3.0 * Math.Pow(x, 3.0 / 2.0);
                SolutionPair p = ModifiedBessel(1.0 / 3.0, z);
                double a = 1.0 / (Global.SqrtThree * Math.PI);
                double b = 2.0 / Global.SqrtThree * p.FirstSolutionValue + p.SecondSolutionValue / Math.PI;
                double sx = Math.Sqrt(x);
                return (new SolutionPair(
                    a * sx * p.SecondSolutionValue, a * (x * p.SecondSolutionDerivative + p.SecondSolutionValue / sx / 2.0),
                    sx * b, x * (2.0 / Global.SqrtThree * p.FirstSolutionDerivative + p.SecondSolutionDerivative / Math.PI) + b / sx / 2.0
                ));
            } else if (x < 108.0) {
                SolutionPair s = Airy_Asymptotic_Positive_Scaled(x, out double xi);
                double e = Math.Exp(xi);
                return (new SolutionPair(s.FirstSolutionValue / e, s.FirstSolutionDerivative / e, s.SecondSolutionValue * e, s.SecondSolutionDerivative * e));
            } else if (x <= Double.PositiveInfinity) {
                return (new SolutionPair(0.0, -0.0, Double.PositiveInfinity, Double.PositiveInfinity));
            } else {
                Debug.Assert(Double.IsNaN(x));
                return (new SolutionPair(x, x, x, x));
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
        public static double AiryAiZero (int k) {
            return (BesselMath.AiryAiZero(k));
        }

        /// <summary>
        /// Computes the requested zero of the Airy Bi (Bairy) function.
        /// </summary>
        /// <param name="k">The index of the zero.</param>
        /// <returns>The <paramref name="k"/>th value of x for which Bi(x) = 0.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> is less than 1.</exception>
        public static double AiryBiZero (int k) {
            return (BesselMath.AiryBiZero(k));
        }

    }

    public static partial class AdvancedComplexMath {

        public static Complex AiryAi_Series (Complex z) {

            Complex p = AdvancedMath.Ai0;
            Complex q = AdvancedMath.AiPrime0 * z;
            Complex f = p + q;

            Complex z3 = ComplexMath.Sqr(z) * z;
            for (int k = 0; k < Global.SeriesMax; k += 3) {
                Complex f_old = f;
                p *= z3 / ((k + 2) * (k + 3));
                q *= z3 / ((k + 3) * (k + 4));
                f += p + q;
                if (f == f_old) {
                    return (f);
                }
            }
            throw new NonconvergenceException();

        }

        public static Complex AiryAi_Asymptotic (Complex z) {

            Debug.Assert(ComplexMath.Abs(z) >= 9.0);

            if (z.Re >= 0.0) {
                Complex xi = 2.0 / 3.0 * ComplexMath.Pow(z, 3.0 / 2.0);
                Airy_Asymptotic_Subseries(xi, out Complex u0, out Complex v0, out Complex u1, out Complex v1);
                Complex e = ComplexMath.Exp(xi);
                Complex q = ComplexMath.Pow(z, 1.0 / 4.0);
                return (0.5 / Global.SqrtPI / q / e * u1);
            } else {
                z = -z;
                Complex xi = 2.0 / 3.0 * ComplexMath.Pow(z, 3.0 / 2.0);
                Airy_Asymptotic_Subseries(xi, out Complex u0, out Complex v0, out Complex u1, out Complex v1);
                Complex c = ComplexMath.Cos(xi);
                Complex s = ComplexMath.Sin(xi);

                throw new NotImplementedException();
            }
     
        }

        private static void Airy_Asymptotic_Subseries (Complex xi, out Complex u0, out Complex v0, out Complex u1, out Complex v1) {

            Complex du = 1.0;

            u0 = 1.0;
            v0 = 1.0;
            u1 = 1.0;
            v1 = 1.0;

            for (int k = 1; k < Global.SeriesMax; k++) {

                Complex u0_old = u0;
                Complex u1_old = u1;

                du *= (6 * k - 5) * (6 * k - 3) * (6 * k - 1) / (k * (2 * k - 1) * 216.0) / xi;
                Complex dv = (6 * k + 1) * du / (1 - 6 * k);

                u0 += du;
                v0 += dv;
                if (k % 2 == 0) {
                    u1 += du;
                    v1 += dv;
                } else {
                    u1 -= du;
                    v1 -= dv;
                }

                if (u0_old == u0 && u1_old == u1) {
                    return;
                }

            }

            throw new NonconvergenceException();

        }

    }
}
