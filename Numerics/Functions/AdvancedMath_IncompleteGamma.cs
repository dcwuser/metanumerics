using System;
using System.Diagnostics;
using System.IO;
using System.Text;

namespace Meta.Numerics.Functions {
    public static partial class AdvancedMath {

        // Incomplete Gamma is covered in DLMF Chapter 8 (https://dlmf.nist.gov/8) and NR Chapter 6.2

        // Main functions a small (left, lower) gamma, big (right, upper) gamma, and their regularized versions P and Q.

        // The error function is a special case.

        // We should also look into negative values and inversion.

        // "Efficient and accurate algorithms for the computation and inversion of the incomplete gamma function ratios"
        // https://arxiv.org/abs/1306.1754
        // Division of 1st quadrant, \alpha(x), Q series, another aysmptotic fomrulation

        // "Computation of the incomplete gamma function for negative values of the argument" 2016
        // https://arxiv.org/abs/1608.04152

        /// <summary>
        /// Computes the normalized lower (left) incomplete gamma function.
        /// </summary>
        /// <param name="a">The shape parameter, which must be non-negative.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of P(a,x) = &#x3B3;(a,x) / &#x393;(a).</returns>
        /// <remarks><para>The incomplete Gamma function is obtained by carrying out the Gamma function integration from zero to some
        /// finite value x, instead of to infinity. The function is normalized by dividing by the complete integral, so the
        /// function ranges from 0 to 1 as x ranges from 0 to infinity.</para>
        /// <para>For sufficiently high values of x, this function becomes 1 within floating point precision. To determine its deviation from 1
        /// in this region, use the complementary function <see cref="RightRegularizedGamma"/>.</para>
        /// <para>For a=&#x3BD;/2 and x=&#x3C7;<sup>2</sup>/2, this function is the CDF of the &#x3C7;<sup>2</sup> distribution with &#x3BD; degrees of freedom.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="a"/> is negative, or <paramref name="x"/> is negative.</exception>
        /// <seealso cref="RightRegularizedGamma" />
        public static double LeftRegularizedGamma(double a, double x) {
            RegularizedIncompleteGamma(a, x, out double P, out _);
            return P;
        }

        /// <summary>
        /// Computes the normalized upper (right) incomplete gamma function.
        /// </summary>
        /// <param name="a">The shape parameter, which must be positive.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of Q(a, x) = &#x393;(a,x) / &#x393;(a).</returns>
        /// <remarks>
        /// <para>This function is the complement of the left regularized incomplete gamma function <see cref="LeftRegularizedGamma"/>.
        /// Their values sum to one.</para>
        /// <para>For sufficiently low values of x, this function becomes 1 within floating point precision. To determine its deviation from 1
        /// in this region, use the complementary function <see cref="LeftRegularizedGamma"/>.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="a"/> is negative, or <paramref name="x"/> is negative.</exception>
        /// <seealso cref="LeftRegularizedGamma"/>
		public static double RightRegularizedGamma(double a, double x) {
            RegularizedIncompleteGamma(a, x, out _, out double Q);
            return Q;
        }

        // NR says to transition between P series and Q continue fraction at x ~ a + 1.
        // For a > 1 this is fine (and may indeed be computationally more efficient computationally).
        // But the actual median m where P = Q = 1/2 is lower, a - 1/3 < m < a, and for a < 1 this
        // can have a impact on accuracy because we end up computing P ~ 1 when we should be computing
        // Q << 1. Transition at x ~ a better, but still not good enough for a << 1,
        // for which m ~ e^{-\gamma} 2^{-1/a} << a. Temme suggests using a ~ ln(1/2) / ln(1/2 x) for x < 1,
        // which accounts for fact that median occurs at x << a for a << 1. But then there is a new
        // problem: CF for Q converges poorly to not at all for x < 1. (NR's criteria ensured that
        // we never computed Q for x < 1 and maybe that's why they chose it.) So we need a new
        // method for Q when x < 1 and a < m. Temme suggested using what we call the Q series.

        // There is a series for P good for sufficiently small x. There is a continued fraction
        // for Q good for sufficiently large x. For a > 1 but not too big, these
        // together cover all x. For big a, through, there is a area in between where neither
        // converges practically. In that regime, gamma distribution is approximately Gaussian
        // with \mu ~ a, \sigma ~ \sqrt{a}. Temme gives a method that starts from Gaussian and corrects it.
        // For a < 1, the Q continued fraction does not converge so we use a series for Q derived
        // via analytic manipulations from an alternative series for P.

        private static readonly double expmeg = Math.Exp(-AdvancedMath.EulerGamma);

        private static double EstimatedGammaMedian (double a) {
            // For a > 1, median \nu \approx a - 1/3. For a << 1, \nu \approx e^{-\gamma} 2^{-1/a}.
            // These meet at (a,\nu) = (0.46, 0.13) (with Q ~ 0.58), so transition between them there.
            return a < 0.46 ? expmeg * Math.Pow(2.0, -1.0 / a) : a - 1.0 / 3.0;
        }

        private static void RegularizedIncompleteGamma(double a, double x, out double P, out double Q) {
            if (a < 0.0) throw new ArgumentOutOfRangeException(nameof(a));
            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));
            if ((a > 128.0) && (Math.Abs(x - a) < 0.25 * a)) {
                Gamma_Temme(a, x, out P, out Q);
            } else {
                // We should compute the smaller of P or Q, and the median seperates these two regimes.
                if (x < EstimatedGammaMedian(a)) {
                    // left of median, so compute P
                    P = AdvancedMath.PoissonProbability(x, a) * IncompleteGamma_PSeries(a, x);
                    Q = 1.0 - P;
                } else {
                    // right of median, so compute Q
                    if (x > 1.5) {
                        Q = a * AdvancedMath.PoissonProbability(x, a) * IncompleteGamma_QContinuedFraction(a, x);
                    } else {
                        Q = RegularizedIncompleteGammaQ_Series(a, x);
                    }
                    P = 1.0 - Q;
                }
            }
            Debug.Assert(0.0 <= P && P <= 1.0);
            Debug.Assert(0.0 <= Q && Q <= 1.0);
        }

        /// <summary>
        /// Computes the lower (left) incomplete (small) gamma function.
        /// </summary>
        /// <param name="a">The shape parameter, which must be non-negative.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of &#x3B3;(a,x).</returns>
        /// <remarks><para>The incomplete Gamma function is obtained by carrying out the Gamma function integration from zero to some
        /// finite value x, instead of to infinity. It therefore ranges from zero to &#x3B3;(a) as x varies from zero to infinity.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="a"/> is negative, or <paramref name="x"/> is negative.</exception>
        /// <seealso cref="LeftRegularizedGamma" />
        /// <seealso cref="UpperIncompleteGamma"/>
        public static double LowerIncompleteGamma(double a, double x) {
            if (a < 0.0) throw new ArgumentOutOfRangeException(nameof(a));
            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));
            if ((a > 128.0) && (Math.Abs(x - a) < 0.25 * a)) {
                Gamma_Temme(a, x, out double P, out _);
                return AdvancedMath.Gamma(a) * P;
            } else {
                if (x < EstimatedGammaMedian(a)) {
                    return LowerIncompleteGamma_Series(a, x);
                } else {
                    double completeGamma = AdvancedMath.Gamma(a);
                    double G;
                    if (x > 1.5) {
                        G = UpperIncompleteGamma_ContinuedFraction(a, x);
                    } else {
                        G = UpperIncompleteGamma_Series(a, x);
                    }
                    Debug.Assert(0.0 <= G && G <= completeGamma);
                    return completeGamma - G;
                }
            }
        }

        /// <summary>
        /// Computes the upper (right) incomplete (big) gamma function.
        /// </summary>
        /// <param name="a">The shape parameter, which must be non-negative.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of &#x393;(a,x).</returns>
        /// <remarks><para>The upper incomplete Gamma function is obtained by carrying out the Gamma function integration from x to infinity.
        /// It therefore ranges from &#x393;(a) to zero as x varies from zero to infinity.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="a"/> is negative, or <paramref name="x"/> is negative.</exception>
        /// <seealso cref="RightRegularizedGamma" />
        /// <seealso cref="LowerIncompleteGamma"/>
        public static double UpperIncompleteGamma(double a, double x) {
            if (a < 0.0) throw new ArgumentOutOfRangeException(nameof(a));
            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));
            if ((a > 128.0) && (Math.Abs(x - a) < 0.25 * a)) {
                Gamma_Temme(a, x, out _, out double Q);
                return AdvancedMath.Gamma(a) * Q;
            } else {
                if (x < EstimatedGammaMedian(a)) {
                    double bigGamma = AdvancedMath.Gamma(a);
                    double littleGamma = LowerIncompleteGamma_Series(a, x);
                    Debug.Assert(0.0 <= littleGamma && littleGamma <= bigGamma);
                    return bigGamma - littleGamma;
                } else {
                    // G is smaller, so compute G
                    if (x > 1.5) {
                        return UpperIncompleteGamma_ContinuedFraction(a, x);
                    } else {
                        return UpperIncompleteGamma_Series(a, x);
                    }
                }
            }
        }

        // NR 6.2.5, DLMF 8.7.1	
        // P(a, x) = \gamma(a, x) \Gamma(a) given by a series in x
        //   P(a, x) = \frac{\gamma(a, x)}{\Gamma(a)} = e^{-x} x^a \sum_{k=0} \frac{x^k}{\Gamma(a + 1 + k)}
        //           = e^{-x} x^a \left[ \frac{1}{\Gamma(a + 1)} + \frac{x}{\Gamma(a+2)} + \frac{x^2}{\Gamma(a+3)} + \cdots \right]
        //           = \frac{e^{-x} x^a}{\Gamma(a + 1)} \left[ 1 + \frac{x}{(a+1)} + \frac{x^2}{(a+1)(a+2)} + \cdots \right]
        // For P pre-factor is poisson probability. For \gamma(a, x) = P(a, x) \Gamma(a), pre-factor is e^{-x} x^a / a.
        // We have specialized code for either pre-factor.

        private static double LowerIncompleteGamma_Series(double a, double x) {
            return AdvancedMath.ExpTimesPower(x, a) / a * IncompleteGamma_PSeries(a, x);
        }

        // S_P(a, x) = \sum_{k=0}^{\infty} x^k / (a + 1)_k = 1 + \frac{x}{a + 1} + \frac{x^2}{(a+1)(a+2)} + \cdots
        // Note S_P(a, 0) = 1, and the value grows as x increases.

        // \gamma(a, x) = \frac{e^{-x} a^x}{a} S_P(a, x), and since P(a, x) = \frac{\gamma(a,x)}{\Gamma(a)},
        // P(a, x) = e^{-x} = \frac{e^{-x} a^x}{\Gamma(a+1)} S_P(a,x)

        private static double IncompleteGamma_PSeries(double a, double x) {
            double ap = a + 1.0;
            double ds = x / ap;
            double s = 1.0 + ds;
            for (int k = 2; k < Global.SeriesMax; k++) {
                double s_old = s;
                ap += 1.0;
                ds *= x / ap;
                s += ds;
                if (s == s_old) {
                    return s;
                }
            }
            throw new NonconvergenceException();
        }

        // Continued fraction for Q from NR 6.2.7
        //    Q(a,x) = \frac{e^{-x} x^a}{\Gamma(a)} \left( \frac{1}{x - a + 1 -} \frac{1 (1 - a)}{x - a + 3 -} \frac{2 (2 - a)}{x - a + 5 -} \right)
        // Covergence depends mostly on x (~25 terms for x~4, ~45 terms for x~2, ~90 terms for x~1, ~150 terms for x~1/2).
        // For Q, use the pre-factor \frac{e^{-x} x^a}{\Gamma(a)}. For \Gamma, use the pre-factor e^{-x} x^a.

        private static double UpperIncompleteGamma_ContinuedFraction(double a, double x) {
            if (x == Double.PositiveInfinity) return 0.0;
            return AdvancedMath.ExpTimesPower(x, a) * IncompleteGamma_QContinuedFraction(a, x);
        }

        private static double IncompleteGamma_QContinuedFraction(double a, double x) {

            // Entering this loop with bb infinite (as caused e.g. by infinite x) will cause a
            // NonconvergenceException instead of the expected convergence to zero, so detect
            // and return early.
            if (Double.IsPositiveInfinity(x)) return 0.0;

            double aa = 1.0;            // a_1
            double bb = x - a + 1.0;	// b_1
            Debug.Assert(bb > 0.0);
            double D = 1.0 / bb;        // D_1 = b_0 / b_1
            double Df = aa / bb;        // Df_1 = f_1 - f_0
            double f = Df;              // f_1 = f_0 + Df_1 = b_0 + Df_1
            for (int k = 1; k < Global.SeriesMax; k++) {
                double f_old = f;
                aa = -k * (k - a);
                bb += 2.0;
                D = 1.0 / (bb + aa * D);
                Df = (bb * D - 1.0) * Df;
                f += Df;
                if (f == f_old) {
                    return f;
                }
            }
            throw new NonconvergenceException();
        }


        // DLMF 8.7.1
        // The continued fraction does not converge for x < ~1, so we need an alternative.
        // Start from an alternative series for P.
        //   P(a, x) = \frac{x^a}{\Gamma(a)} \sum_{k=0}^{\infty} \frac{(-x)^k}{k! (a + k)}
        // For a -> 0,  the first term x^{a} / \Gamma(a + 1) -> 1, so remaining terms tend
        // toward a series for Q.

        // More carefully,
        //   Q(a, x) = 1 - P(a, x) = Q1(a, x) + Q2(a, x)
        // with
        //   Q1(a, x) = 1 - \frac{x^a}{\Gamma(a + 1)}
        //   Q2(a, x) = - \frac{x^a}{\Gamma(a)} \sum{k=1}^[\infty} \frac{(-x)^k}{k! (a + k)}
        
        // To compute Q1, use
        //   \frac{x^a}{\Gamma(a + 1)} = \exp{a \ln x - \ln\Gamma(1 + a)}
        // All terms proprotional to a for a small, so we do a->0 well, and
        //   Q1(a, x) = -expm1( a ln(x) - LogGammaOnePlus(a) )
        // To compute Q2, just use the series.


        // We need to accurately compute (i.e. without loss of accuracy due to subtraction)
        // how much the first term differs from 1.
        //   1 - \frac{x^a}{\Gamma(a + 1)} = \frac{\Gamma(a + 1) - x^a}{\Gamma(a + 1)}
        // Concentrate on numerator.
        //   \Gamma(a + 1) - x^a = [\Gamma(a + 1) - 1] + [1 - x^a]
        //   = [e^{\ln\Gamma(a + 1)} - 1] + [1 - e^{a \ln x}]
        //   = expm1(\ln\Gamma(a + 1)) - expm1(a \ln x)
        // For small a, both arguments are proportional to and can be accurately computed.

        // For \gamma = Q \Gamma(a), we need above result divided by a. For this to work
        // for a = 0, we need to divide analytically. Introduce "reduced" functions,
        //   expm1(z) = z rxpm1(z)
        //   \ln\Gamma(z + 1) = z rlng1p(z)
        // so previous numerator becomes
        //   a rlng1p(a) rexpm1(a rlng1p(a)) - a ln(x) rexpm1(a ln(x))
        // which can be divided analytically by a.

        private static void Stuff (double a, double x, out double F, out double gamma_one_plus_a, out double x_to_a) {
            double t1 = Math.Log(x);
            double t2 = GammaSeries.ReducedLogGammaOnePlus(a);
            double t = t1 - t2;
            F = -MoreMath.ReducedExpMinusOne(a * t) * t;
            x_to_a = Math.Exp(a * t1);
            gamma_one_plus_a = Math.Exp(a * t2);
        }

        private static double RegularizedIncompleteGammaQ_Series(double a, double x) {
            Debug.Assert(0.0 <= a && a <= 2.5);
            Debug.Assert(0.0 <= x && x <= a + 1.5);
            double t1 = a * Math.Log(x);
            double t2;
            if (a <= 0.5) {
                t2 = GammaSeries.LogGammaOnePlus(a);
            } else if (a <= 1.5) {
                t2 = GammaSeries.LogGammaTwoPlus(a - 1.0);
            } else {
                t2 = LogGamma(a + 1.0); // Do better here to avoid multiple cancelling adds/subtracts
            }
            double t = t1 - t2; // There is some cancellation here if x > 1.
            double q1 = -MoreMath.ExpMinusOne(t);
            double q2 = a * Math.Exp(t) * IncompleteGamma_QSeries(a, x);
            return q1 + q2;
            /*
                double lng1pa = (a <= 0.5) ? GammaSeries.LogGammaOnePlus(a) : GammaSeries.LogGammaTwoPlus(a - 1.0);
            double u1 = MoreMath.ExpMinusOne(lng1pa);
            double u2 = -MoreMath.ExpMinusOne(alnx);
            double u = (u1 + u2) / Math.Exp(lng1pa);
            return u + a * AdvancedMath.PowerOverFactorial(x, a) * IncompleteGamma_QSeries(a, x);
            */
        }

        private static double UpperIncompleteGamma_Series(double a, double x) {
            Debug.Assert(0.0 <= a && a <= 1.5);
            double lnx = Math.Log(x);
            double rlng1pa = (a <= 0.5) ? GammaSeries.ReducedLogGammaOnePlus(a) : GammaSeries.LogGammaTwoPlus(a - 1.0) / a;
            double u1 = rlng1pa * MoreMath.ReducedExpMinusOne(a * rlng1pa);
            double u2 = -lnx * MoreMath.ReducedExpMinusOne(a * lnx);
            double u = u1 + u2;
            return u + Math.Pow(x, a) * IncompleteGamma_QSeries(a, x);
        }

        private static double IncompleteGamma_QSeries(double a, double x) {
            double t = x;
            double s = t / (a + 1.0);
            for (int k = 2; k < Global.SeriesMax; k++) {
                double s_old = s;
                t *= -x / k;
                s += t / (a + k);
                if (s == s_old) {
                    return s;
                }
            }
            throw new NonconvergenceException();
        }

        // For large a, for x ~ a, the convergence of both the series and the continued fraction for the incomplete gamma
        // function is very slow; the problem is that the kth term goes like x/(a+k), and for large a, adding k
        // makes little difference until k ~ a.

        // In this region, NR uses ~15-point Gaussian quadrature of the peaked integrand, which should be about as good
        // as one iteration of the 15-point Gauss-Kronrod integrator used by our adaptive integrator. But when I tried
        // our adaptive integrator, it wanted to subdivide and repeat, requiring hundreds of evaluations to achieve full
        // accuracy. This leads me to believe that the NR algorithm is unlikely to achieve full accuracy.

        // So in this region we use instead a rather strange expansion due to Temme. It looks simple enough at first:
        //   P = \frac{1}{2} \erfc(-\zeta) -  R     Q = \frac{1}{2} \erfc(\zeta) / 2 + R
        // where \zeta \approx \frac{x-a}{\sqrt{2a}}, corresponding to the Normal(a,\sqrt{a}) approximation,
        // and R is a correction series.

        // The first oddity is that \zeta is not quite \frac{x-a}{\sqrt{2a}}. Instead \zeta^2 = \eta^2 a / 2, where
        //   \frac{1}{2} \eta^2 = D(x/a) = x/a - 1 - \ln(x/a)
        // or, in terms of z = \frac{x-a}{a} = x/a - 1,
        //   \frac{1}{2} \eta^2 = D(1 + z) = z - \ln(1 + z)
        // This last form is computationally useful because we can compute D accurately as z \rightarrow 0 using ln1p.

        // R can be expressed as a series in 1/a
        //   R = \frac{\exp(-\zeta^2)}{\sqrt{2 \pi a}} \sum_{k=0}^{\infty} \frac{C_k(\eta)}{a^k}
        // Here the C_k are functions of \eta which Temme derives. Unfortunately, their analytic expressions
        // suffer from multiple canceling singularities as \eta \rightarrow 0. Therefore they must be expanded analytically
        // in their argument before being used numerically.

        // In our experiments, the expansion variable the exhibited fastest convergence was u = \ln(x / a) = \ln(1 + z),
        // in terms of which
        //   z = \exp(u) - 1
        //   \frac{1}{2} \eta^2 = \exp(u) - 1 - u  
        //   R = \frac{\exp(-z^2)}{\sqrt{2 \pi a}} \sum_{n=0}^{\infty} \sum_{m=0}^{\infty} D_{n,m} \frac{u^m}{a^n}
        // So effectively we compute a double power series in u and 1/a.

        // The first coefficient function is
        //   C_0 = \frac{1}{z} - \frac{1}{\eta}
        //       = \frac{1}{e^u - 1} - \frac{1}{\sqrt{2(e^u - 1 - u)}}
        //       = -1/3 + 1/12 u - 1/1080 u^2 - 19/12960 u^3 + \cdots
        // Higher C can be derived using
        //   \eta C_n = \frac{d}{d\eta} C_{n-1} + \frac{\eta}{z} \gamma_n
        // where \gamma_n are coefficients that appear in the asymptotic expansion of \Gamma^*. For example
        //   C_1 = \frac{1}{\eta^3} - \frac{1}{z^3} - \frac{1}{z^2} - \frac{1}{12z}
        //       = -1/540 - 1/288 u + 25/12096 u^2 + \cdots

        // There turns out to be a bit of a shortcut for computing these series. At every order
        //   C_n = (-1)^{n+1} \frac{(2n-1)!!}{\eta^{2n+1}} + \frac{1}{z} terms
        // and the 1/z terms are only there to cancel the singularities as z->0 of the 1/\eta term;
        // all the real information is encoded in the non-singular part of the 1/\eta term.
        // So we can just expand the 1/\eta term in powers of z, drop all the singular terms, then
        // re-expand the non-singular terms in terms of u using z = e^u - 1. This is how we have derived
        // the D_{n,m}. (Note that this is not the same as directly expanding in u and dropping singular terms;
        // if you do that, you get wrong results, because the 1/z terms do contribute terms non-singular in u.)

        // We record enough terms to obtain full accuracy when a > 100 and z < 0.25.

        // https://www.ams.org/journals/mcom/1975-29-132/S0025-5718-1975-0387674-2/S0025-5718-1975-0387674-2.pdf (?)
        // https://www.researchgate.net/publication/243594366_The_Asymptotic_Expansion_of_the_Incomplete_Gamma_Functions
        // https://www.ams.org/journals/mcom/1992-58-198/S0025-5718-1992-1122079-8/S0025-5718-1992-1122079-8.pdf

        // Note https://arxiv.org/abs/1803.07841 shows an asymptotic series in 1/\sqrt{a} with more easily evaluated
        // coefficients. I tested it and the convergence, as one might have expected, really is much poorer than Temme's.

        private static readonly double[][] TemmeD = new double[][] {
            new double[] { - 1.0 / 3.0, 1.0 / 12.0, - 1.0 / 1080.0, - 19.0 / 12960.0, 1.0 / 181440.0, 47.0 / 1360800.0,
                1.0 / 32659200.0, - 221.0 / 261273600.0, - 281.0 / 155196518400.0,  857.0 / 40739086080.0, 1553.0 / 40351094784000.0 },
            new double[] { -1.0 / 540.0, - 1.0 / 288.0, 25.0 / 12096.0, - 223.0 / 1088640.0, - 89.0 / 1088640.0,
                757.0 / 52254720.0,  445331.0 / 155196518400.0, - 1482119.0 / 2172751257600.0, - 7921307.0 / 84737299046400.0 },
            new double[] { 25.0 / 6048.0, - 139.0 / 51840.0, 101.0 / 311040.0, 1379.0 / 7464960.0, - 384239.0 / 7390310400.0,
                - 1007803.0 / 155196518400.0, 88738171.0 / 24210656870400.0, 48997651.0 / 484213137408000.0 },
            new double[] { 101.0 / 155520.0, 571.0 / 2488320.0, - 3184811.0 / 7390310400.0, 36532751.0 / 310393036800.0,
                10084279.0 / 504388684800.0, - 82273493.0 / 5977939968000.0 },
            new double[] { - 3184811.0 / 3695155200.0, 163879.0 / 209018880.0, - 2745493.0 / 16303472640.0,
                - 232938227.0 / 2934625075200.0, 256276123.0 / 5869250150400.0 },
            new double[] { - 2745493.0 / 8151736320.0, - 5246819.0 / 75246796800.0, 119937661.0 / 451480780800.0,
                - 294828209.0 / 2708884684800.0 },
            new double[] { 119937661.0 / 225740390400.0, - 534703531.0 / 902961561600.0 },
            new double[] { 8325705316049.0 / 24176795811840000.0 , 4483131259.0 / 86684309913600.0 }
        };

        private static void Gamma_Temme(double a, double x, out double P, out double Q) {

            double u = Math.Log(x / a);

            // compute argument of error function, which is almost (x-a)/sqrt(a)

            double zee = (x - a) / a;
            Debug.Assert(Math.Abs(zee) <= 0.25);
            double zetaSquared = a * (zee - MoreMath.LogOnePlus(zee));
            double zeta;
            if (zee >= 0.0) {
                zeta = Math.Sqrt(zetaSquared);
                Q = 0.5 * AdvancedMath.Erfc(zeta);
                P = 1.0 - Q;
            } else {
                zeta = -Math.Sqrt(zetaSquared);
                P = 0.5 * AdvancedMath.Erfc(-zeta);
                Q = 1.0 - P;
            }

            double dz = 1.0;
            double z = dz;
            for (int i = 3; true; i++) {
                if (i > Global.SeriesMax) throw new NonconvergenceException();
                double z_old = z;
                dz *= u / i;
                z += dz;
                if (z == z_old) break;
            }
            z = u * Math.Sqrt(a * z / 2.0);

            // the first approximation is just the almost-Gaussian one

            if (z > 0) {
                Q = AdvancedMath.Erfc(z) / 2.0;
                P = 1.0 - Q;
            } else {
                P = AdvancedMath.Erfc(-z) / 2.0;
                Q = 1.0 - P;
            }

            // compute Temme's correction to the Gaussian approximation

            double R0 = Math.Exp(-z * z) / Math.Sqrt(Global.TwoPI * a);

            double S0 = 0.0;
            double ai = 1.0;
            for (int i = 0; i < TemmeD.Length; i++) {
                double dS = 0.0;
                double uj = 1.0;
                for (int j = 0; j < TemmeD[i].Length; j++) {
                    dS += TemmeD[i][j] * uj;
                    uj *= u;
                }
                S0 += dS / ai;
                ai *= a;
            }

            double R = R0 * S0;
            Q += R;
            P -= R;

        }

    }
}
