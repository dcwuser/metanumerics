using System;
using System.Diagnostics;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Functions {
    public static partial class AdvancedMath {

        /// <summary>
        /// Computes the Gauss hypergeometric function. 
        /// </summary>
        /// <param name="a">The first upper parameter.</param>
        /// <param name="b">The second upper parameter.</param>
        /// <param name="c">The lower parameter.</param>
        /// <param name="x">The argument, which must be less than or equal to one.</param>
        /// <returns>The value of <sub>2</sub>F<sub>1</sub>(a, b; c; x).</returns>
        /// <remarks>
        /// <para>The Gauss Hypergeometric function is defined by the hypergeometric series:</para>
        /// <img src="../images/Hypergeometric2F1.png" />
        /// <para>For generic values of a, b, and c, the Gauss hypergeometric function becomes complex for x > 1.
        /// However, there are specific cases, most commonly for negative integer values of a and b, for which the
        /// function remains real in this range.</para>
        /// <para>For some arguments, the value of the function depends on the order in which the arguments are
        /// regarded as approaching their values. For example, in general, if x = 0, F = 1. On the other hand,
        /// in general, if c = -1, F is infinite. If x = 0 and c = -1, the value depends on which limit is taken first.
        /// </para>
        /// <para>Our implementation does not achieve full precision in all parameter regions.
        /// For |x| &lt; 1/2 and |a|, |b| &lt; 10, we achieve approximately full precision.
        /// For |x| &gt; 1/2, we loose about one decimal digit of precision, and for
        /// every order of magnitude that |a| or |b| exceeds 10, we loose about one decimal digit of precision.
        /// </para>
        /// </remarks>
        /// <seealso href="https://en.wikipedia.org/wiki/Hypergeometric_function"/>
        public static double Hypergeometric2F1 (double a, double b, double c, double x) {

            // The correct order of all these limits is debatable.

            if (x == 0.0) return (1.0);
            if ((a == 0.0) || (b == 0.0)) return (1.0);

            // Check for the polynomial case.
            double polynomial_n = -1.0;
            double polynomial_b = 0.0;
            if (IsNonPositiveInteger(a)) {
                polynomial_n = -a;
                polynomial_b = b;
            }
            if (IsNonPositiveInteger(b) && ((polynomial_n < 0.0) || (b > -polynomial_n))) {
                polynomial_n = -b;
                polynomial_b = a;
            }
            if (polynomial_n >= 0) {
                return (Hypergeometric2F1_Polynomial((int) polynomial_n, polynomial_b, c, x));
                // Can we apply a transformation to eliminate the cancelation problem?
            }

            // Check for c canceling with a or b.
            // (This check must occur after the polynomial check to agree with Mathematica's result for -3,-1/2,-3,
            // which appears to be required to fulfill recurrence. Basically, it looks like there is some tension
            // between fulfilling recurrence and fulfilling quadratic transforms for polynomial cases; look into this.)
            if (c == a) return (Math.Pow(1.0 - x, -b));
            if (c == b) return (Math.Pow(1.0 - x, -a));

            // We have dealt with the polynomial case and the case where c cancels with a or b, so
            // if we are here and c is a non-positive-integer, the function blows up.
            if (IsNonPositiveInteger(c)) return (Double.NaN);

            if (x < -1.0) {

                // x -> 1/(1-x) maps (-\inf, -1) -> (0, 1/2)
                double xPrime = 1.0 / (1.0 - x);

                double bma = b - a;
                double m = Math.Round(bma);
                double e = bma - m;

                if (m < 0) {
                    return (Hypergeometric2F1_Series_OneOverOneMinusX(b, -(int) m, -e, c, xPrime));
                } else {
                    return (Hypergeometric2F1_Series_OneOverOneMinusX(a, (int) m, e, c, xPrime));
                }

            } else if (x < -0.25) {

                // x -> x/(x-1) maps (-1, 0) -> (1/2, 0)
                double xPrime = x / (x - 1.0);

                // Use transformed series with smallest a, b values
                if (Math.Max(Math.Abs(a), Math.Abs(c - b)) <= Math.Max(Math.Abs(b), Math.Abs(c - a))) {
                    return (Math.Pow(1.0 - x, -a) * Hypergeometric2F1_Series(a, c - b, c, xPrime));
                } else {
                    return (Math.Pow(1.0 - x, -b) * Hypergeometric2F1_Series(c - a, b, c, xPrime));
                }

            } else if (x <= 0.5) {

                // Use series with smallest a, b values
                if (Math.Max(Math.Abs(a), Math.Abs(b)) <= Math.Max(Math.Abs(c - a), Math.Abs(c - b))) {
                    return (Hypergeometric2F1_Series(a, b, c, x));
                } else {
                    // This transform doesn't work in the polynomial case! (e.g. -1,-3/2,-2,1/2)
                    return (Math.Pow(1.0 - x, c - a - b) * Hypergeometric2F1_Series(c - a, c - b, c, x));
                }

            } else if (x <= 1.0) {

                // x -> 1-x maps (1/2, 1) -> (0, 1/2)
                double xPrime = 1.0 - x;

                // When xPrime is exactly zero, we can have some problems inside the transformed series routine.
                // Really we should deal with this by handling the x = 1 case analytically, but it's actually
                // not trivial to do so. So for now, just use a very very tiny value of xPrime.
                if (x == 1.0) xPrime = 1.0 / Double.MaxValue;

                double cab = c - a - b;
                double m = Math.Round(cab);
                double e = cab - m;

                if (m < 0) {
                    return (Hypergeometric2F1_Series_OneMinusX(c - a, c - b, -(int) m, -e, xPrime) * Math.Pow(xPrime, cab));
                } else {
                    return (Hypergeometric2F1_Series_OneMinusX(a, b, (int) m, e, xPrime));
                }

            } else {

                throw new ArgumentOutOfRangeException(nameof(x));

            }

        }

        private static bool IsNonPositiveInteger (double x) {
            return ((x <= 0.0) && (Math.Round(x) == x));
        }

        private static double Hypergeometric2F1_Series (double a, double b, double c, double x) {

            if ((Math.Abs(a) + 1.0) * (Math.Abs(b) + 1.0) * Math.Abs(x) > 16.0 * (Math.Abs(c) + 1.0)) {
                // We are in danger of not converging. Can we do anything about it? Use a recurrence, for example?
            }

            double df = a * b / c * x;
            double f = 1.0 + df;

            for (int km1 = 1; km1 < Global.SeriesMax; km1++) {
                double f_old = f;
                df *= (a + km1) * (b + km1) / (c + km1) * x / (km1 + 1);
                f += df;
                if (f == f_old) return (f);
            }

            throw new NonconvergenceException();

        }

        private static double Hypergeometric2F1_Polynomial (int m, double b, double c, double x) {

            // This can go wrong because sign of terms can alternate.
            // Is there anything we can do to mitigate this problem?

            Debug.Assert(m >= 0);
            Debug.Assert(!IsNonPositiveInteger(b) || (b <= -m));

            if (IsNonPositiveInteger(c) && (c > -m)) return (Double.NaN);

            double t = (-m) * b / c * x;
            double f = 1.0 + t;

            for (int k = 2; k <= m; k++) {
                int km1 = k - 1;
                t *= (km1 - m) * (b + km1) / (c + km1) * x / k;
                f += t;
            }

            return (f);

        }

        // The point of this method is to compute
        //   e G_{e}(x) = \frac{1}{\Gamma(x)} - \frac{1}{\Gamma(x + e)}
        // To do this, we will re-use machinery that we developed to accurately compute the Pochhammer symbol
        //   (x)_e = \frac{\Gamma(x + e)}{\Gamma(x)}
        // To do this, we use the reduced log Pochhammer function L_{e}(x).
        //   \ln((x)_e) = e L_{e}(x)
        // To see why we developed this function, see the Pochhamer code. The Lanczos apparatus allows us to compute
        // L_{e}(x) accurately, even in the small-e limit. To see how look at the Pochhammer code.

        // To connect G_{e}(x) to L_{e}(x), write
        //   e G_{e}(x) = \frac{(x)_e - 1}{\Gamma(x + e)}
        //              = \frac{\exp(\ln((x)_e)) - 1}{\Gamma(x + e)}
        //              = \frac{\exp(e L_{e}(x)) - 1}{\Gamma(x + e)}
        //     G_{e}(x) = \frac{E_{e}(L_{e}(x))}{\Gamma(x + e)}
        // where e E_{e}(x) = \exp(e x) - 1, which we also know how to compute accurately even in the small-e limit.

        // This deals with G_{e}(x) for positive x. But L_{e}(x) and \Gamma(x + e) still blow up for x or x + e
        // near a non-positive integer, and our Lanczos machinery for L_{e}(x) assumes positive x. To deal with
        // the left half-plane, use the reflection formula
        //   \Gamma(z) \Gamma(1 - z) = \frac{\pi}{\sin(\pi z)}
        // on both Gamma functions in the definition of G_{e}(x) to get
        //   e G_{e}(x) = \frac{\sin(\pi x)}{\pi} \Gamma(1 - x) - \frac{\sin(\pi x + \pi e)}{\pi} \Gamma(1 - x - e)
        // Use the angle addition formula on the second \sin and the definition of the Pochhammer symbol
        // to get all terms proportional to one Gamma function with a guaranteed positive argument.
        //   \frac{e G_{e}(x)}{\Gamma(1 - x - e)} =
        //     \frac{\sin(\pi x)}{\pi} \left[ (1 - x - e)_{e} - \cos(\pi e) \right] -  \frac{\cos(\pi x) \sin(\pi e)}{\pi}
        // We need the RHS ~e to for small e. That's manifestly true for the second term because of the factor \sin(\pi e).
        // It's true for the second term because (1 - x - e)_{e} and \cos(\pi e) are both 1 + O(e), but to avoid cancelation
        // we need to make it manifest. Write
        //   (y)_{e} = \exp(e L_{e}(y)) - 1 + 1 = e E_{e}(L_{e}(y)) + 1
        // and
        //   1 - \cos(\pi e) = 2 \sin^2(\half \pi /e)
        // Now we can divide through by e.
        //   \frac{G_{e}(x)}{\Gamma(y)} =
        //     \sin(\pi x) \left[ \frac{E_{e}(L_{e}(y))}{\pi} + \frac{\sin^2(\half \pi e)}{\half \pi e} \right] -
        //     \cos(\pi x) \frac{\sin(\pi e)}{\pi e}
        // and everything can be safely computed.

        // This is a different approach than the one in the Michel & Stoitsov paper. Their approach also used
        // the Lanczos evaluation of the Pochhammer symbol, but had some deficiencies.
        // For example, for e <~ 1.0E-15 and x near a negative integer, it gives totally wrong
        // answers, and the answers loose accuracy for even larger e. This is because the
        // computation relies on a ratio of h to Gamma, both of which blow up in this region.

        private static double NewG(double x, double e) {

            Debug.Assert(Math.Abs(e) <= 0.5);

            // It would be better to compute G outright from Lanczos, rather than via h. Can we do this?

            // Also, we should probably pick larger of 1 - x and 1 - x - e to use as argument of
            // Gamma function factor.

            double y = x + e;
            if ((x < 0.5) || (y < 0.5)) {

                double h = MoreMath.ReducedExpMinusOne(Lanczos.ReducedLogPochhammer(1.0 - y, e), e);

                if (e == 0.0) {
                    double t = MoreMath.SinPi(x) * h / Math.PI - MoreMath.CosPi(x);
                    return (AdvancedMath.Gamma(1.0 - y) * t);
                } else {
                    double s = MoreMath.SinPi(e) / (Math.PI * e);
                    double s2 = MoreMath.Sqr(MoreMath.SinPi(e / 2.0)) / (Math.PI * e / 2.0);
                    double t = MoreMath.SinPi(x) * (h / Math.PI + s2) - MoreMath.CosPi(x) * s;
                    return (AdvancedMath.Gamma(1.0 - y) * t);
                }

            }

            return (MoreMath.ReducedExpMinusOne(Lanczos.ReducedLogPochhammer(x, e), e) / AdvancedMath.Gamma(x + e));

        }

        // Our approach to evaluating the transformed series is taken from Michel & Stoitsov, "Fast computation of the
        // Gauss hypergeometric function with all its parameters complex with application to the Poschl-Teller-Ginocchio
        // potential wave functions" (https://arxiv.org/abs/0708.0116). Michel & Stoitsov had a great idea, but their
        // exposition leaves much to be desired, so I'll put in a lot of detail here.

        // The basic idea is an old one: use the linear transformation formulas (A&S 15.3.3-15.3.9) to map all x into
        // the region [0, 1/2]. The x -> (1-x) transformation, for example, looks like
        //   F(a, b, c, x) =
        //     \frac{\Gamma(c) \Gamma(c-a-b)}{\Gamma(c-a) \Gamma(c-b)} F(a, b, a+b-c+1, 1-x) +
        //     \frac{\Gamma(c) \Gamma(a+b-c)}{\Gamma(a) \Gamma(b)} F(c-a, c-b, c-a-b, 1-x) (1-x)^{c-a-b}

        // When c-a-b is close to an integer, though, there is a problem. Write c = a + b + m + e, where m is a positive integer
        // and |e| <= 1/2. The transformed expression becomes:
        //   \frac{F(a, b, c, x)}{\Gamma(c)} =
        //     \frac{\Gamma(m + e)}{\Gamma(b + m + e) \Gamma(a + m + e)} F(a, b, 1 - m - e, 1 - x) +
        //     \frac{\Gamma(-m - e)}{\Gamma(a) \Gamma(b)} F(b + m + e, a + m + e, 1 + m + e, 1 - x) (1-x)^{m + e}
        // In the first term, the F-function blows up as e->0 (or, if m=0, \Gamma(m+e) blows up), and in the second term 
        // \Gamma(-m-e) blows up in that limit. By finding the divergent O(1/e) and the sub-leading O(1) terms, it's not too
        // hard to show that the divgences cancel, leaving a finite result, and to derive that finite result for e=0.
        // (A&S gives the result, and similiar ones for the divergent limits of other linear transformations.)
        // But we still have a problem for e small-but-not-zero. The pre-limit expressions will have large cancelations.
        // We can't ignore O(e) and higher terms, but developing a series in e in unworkable -- the higher derivatives
        // rapidly become complicated and unwieldy. No expressions in A&S get around this problem,  but we will now
        // develop an approach that does.

        // Notice the divergence of F(a, b, 1 - m - e, 1 - x) is at the mth term, where (1 - m - e)_{m} ~ e. Pull out
        // the finite sum up to the (m-1)th term
        //   \frac{F_0}{\Gamma(c)} = \frac{\Gamma(m+e)}{\Gamma(b + m + e) \Gamma(a + m + e)}
        //     \sum_{k=0}^{m-1} \frac{(a)_k (b)_k}{(1 - m - e)_{k}} \frac{(1-x)^k}{k!}
        // The remainder, which contains the divergences, is:
        //   \frac{F_1}{\Gamma(c)} =
        //     \frac{\Gamma(m + e)}{\Gamma(b + m + e) \Gamma(a + m + e)} \sum_{k=0}^{\infty} \frac{(1-x)^{m + k}}{\Gamma(1 + m + k)}
        //     \frac{\Gamma(a + m + k) \Gamma(b + m +  k) \Gamma(1 - m - e}{\Gamma(a) \Gamma(b) \Gamma(1 - e + k)} +
        //     \frac{\Gamma(-m - e)}{\Gamma(a) \Gamma(b)} \sum_{k=0}^{\infty} \frac{(1-x)^{m + e + k}}{\Gamma(1 + k)}
        //     \frac{\Gamma(b + m + e + k) \Gamma(a + m + e + k) \Gamma(1 + m + e)}{\Gamma(b + m + e) \Gamma(a + m + e) \Gamma(1 + m + e + k)}
        // where we have shifted k by m in the first sum. Use the \Gamma reflection formulae
        //   \Gamma(m + e) \Gamma(1 - m - e) = \frac{\pi}{\sin(\pi(m + e))} = \frac{(-1)^m \pi}{\sin(\pi e)}
        //   \Gamma(-m - e) \Gamma(1 + m + e) = \frac{\pi}(\sin(-\pi(m + e))} = -\frac{(-1)^m \pi}{\sin(\pi e)}
        // to make this
        //   \frac{F_1}{\Gamma(c)} =
        //     \frac{(-1)^m \pi}{\sin(\pi e)} \sum_{k=0}^{\infty} \frac{(1-x)^{m + k}}{\Gamma(a) \Gamma(b) \Gamma(a + m + e) \Gamma(b + m + e)}
        //     \left[ \frac{\Gamma(a + m + k) \Gamma(b + m + k)}{\Gamma(1 + k - e) \Gamma(1 + m + k)} -
        //            \frac{\Gamma(a + m + k + e) \Gamma(b + m + k + e)}{\Gamma(1 + k) \Gamma(1 + m + k + e)} (1-x)^e \right]
        // Notice that \frac{\pi}{\sin(\pi e)} diverges like ~1/e. And that the two terms in parenthesis contain exactly the
        // same products of \Gamma functions, execpt for having their arguments shifted by e. Therefore in the e->0
        // limit their leading terms must cancel, leaving terms ~e, which will cancel the ~1/e divergence, leaving a finite result.

        // We would like to acomplish this cancelation analytically. This isn't too hard to do for e=0. Just write out a Taylor
        // series for \Gamma(z + e), keeping only terms up to O(e). The O(1) terms cancel, the e in front of the O(e) terms gets
        // absorbed into a finite prefactor \frac{\pi e}{\sin(\pi e)}, and we have a finite result. A&S gives the resulting expression.
        // The trouble is for e small-but-not-zero. If we try to evaluate the terms directly, we get cancelations between large terms,
        // leading the catastrophic loss of precision. If we try to use Taylor expansion, we need all the higher derivatives,
        // not just the first one, and the expressions rapidly become so complex and unwieldy as to be unworkable.

        // A good solution, introduced by Forrey, and refined by Michel & Stoistov, is to use finite differences instead
        // of derivatives. If we can express the difference betwen \Gamma(z) and \Gamma(z + e) as a function of z and e that
        // we can compute, then we can analytically cancel the divergent parts and be left with a finite expression involving
        // our finite difference function instead of an infinite series of Taylor series terms. For e=0, the finite difference
        // is just the first derivative, but for non-zero e, it implicitly sums the contributions of all Taylor series terms.

        // The finite difference function to use is:
        //   e G_{e}(z) = \frac{1}{\Gamma(z)} - \frac{1}{\Gamma(z+e)}
        // I played around with a few others, e.g. the perhaps more obvious choice \frac{\Gamma(z+e)}{\Gamma(z)} = 1 + e P_{e}(z),
        // but the key advantage of G_{e}(z) is that it is perfectly finite even for non-positive-integer values of z and z+e,
        // because it uses the recriprocol \Gamma function. (I actually had a mostly-working algorithm using P_{e}(z), but it
        // broke down at non-positive-integer z, because P_{e}(z) itself still diverged for those values.)

        // For a discussion of how to actually compute G_{e}(z), refer to the method notes.

        // The next trick Michel & Stoistov use is to first concentrate just on the k=0 term. The relevent factor is
        //   t = \frac{1}{\Gamma(a) \Gamma(b) \Gamma(a + m + e) \Gamma(b + m + e)}
        //       \left[ \frac{\Gamma(a + m) \Gamma(b + m)}{\Gamma(1 - e) \Gamma(1 + m)} -
        //              \frac{\Gamma(a + m + e) \Gamma(b + m + e)}{\Gamma(1) \Gamma(1 + m + e)} (1-x)^e \right]
        //     = \frac{\Gamma(a + m) \Gamma(b + m)}{\Gamma(a) \Gamma(b)}
        //       \left[ \frac{1}{\Gamma(a + m + e) \Gamma(b + m + e) \Gamma(1 + m) \Gamma(1 - e)} -
        //              \frac{1}{\Gamma(a + m) \Gamma(b + m) \Gamma(1 + m + e) \Gamma(1)} (1-x)^e \right]
        // In the second step, we have put all the \Gamma functions we will need to compute for e-shifted arguments
        // in the denominator, which makes it easier to apply our definition of G_e(z), since there they are also
        // in the deonominator.

        // Before we begin using G_e(z), let's isolate the e-dependence of the (1-x)^e factor. Write
        //   (1-x)^e = \exp(e \ln(1-x) ) = 1 + [ \exp(e \ln(1-x)) - 1 ] = 1 + e E_e(\ln(1-x))
        // where e E_e(z) = \exp(e \ln(1-x)) - 1. We could continue like this, using G_(e) to
        // eliminate every \Gamma(z + e) in favor of G_e(z) and \Gamma(z), but by doing so we
        // would end up with terms containing two and more explicit powers of e, and products
        // of different G_e(z). That would be perfectly correct, but we end up with a nicer
        // expression if we instead "peel off" only one e-shifted function at a time, like this...
        //   \frac{1}{\Gamma(a + m + e) \Gamma(b + m + e) \Gamma(1 + m) \Gamma(1 - e)} =
        //     \frac{1}{\Gamma(a + m + e) \Gamma(b + m + e) \Gamma(1 + m) \Gamma(1)} +
        //     \frac{e G_{-e}(1)}{\Gamma(a + m + e) \Gamma(b + m + e) \Gamma(1 + m)}
        //   \frac{1}{\Gamma(a + m) \Gamma(b + m) \Gamma(1 + m + e)} =
        //     \frac{1}{\Gamma(a + m + e) \Gamma(b + m) \Gamma(1 + m + e)} +
        //     \frac{e G_e(a + m)}{\Gamma(b + m) \Gamma(1 + m + e)}
        //   \frac{1}{\Gamma(a + m + e) \Gamma(b + m) \Gamma(1 + m + e)} =
        //     \frac{1}{\Gamma(a + m + e) \Gamma(b + m + e) \Gamma(1 + m + e)} +
        //     \frac{e G_e(b + m)}{Gamma(a + m + e) \Gamma(1 + m + e)}
        //   \frac{1}{\Gamma(a + m + e) \Gamma(b + m + e) \Gamma(1 + m + e)} =
        //     \frac{1}{\Gamma(a + m + e) \Gamma(b + m + e) \Gamma(1 + m)} -
        //     \frac{e G_e(1 + m)}{\Gamma(a + m + e) \Gamma(b + m + e)}
        // Putting this all together, we have
        //   t = \frac{1}{\Gamma(a + m + e) \Gamma(b + m + e)} \left[ \frac{e G_{-e}(1)}{\Gamma(1 + m) + e G_e(1 + m) \right]
        //     - \frac{1}{\Gamma(1 + m + e)} \left[ \frac{e G_e(a + m)}{\Gamma(b + m)} + \frac{e G_e(b+m)}{\Gamma(a + m + e)} \right]
        //     - \frac{1}{\Gamma(a + m + e) \Gamma(b + m + e) \Gamma(1 + m + e)} e E_e(\ln(1-x))
        // which is, as promised, proportional to e. For some a, b, m, and e, some of these \Gamma functions will blow up, but they are
        // all in the denoninator, so that will just zero some terms. The G_e(z) that appear in the numerator are finite for all z.

        // So now we have the k=0 term. What about higher k terms? We could repeat this analysis, carrying along the k's,
        // and get an expression involving G_e(z) and E_e(z) for each k. Michel & Stoistov's last trick is to realize
        // we don't have to do this, but can instead use our original expressions for each term as a ratio of \Gamma
        // functions to derive a recurrence. Let u_k be the first term, v_k be the second term, so t_k = u_k + v_k.
        // Let r_k = u_{k+1} / u_{k} and s_k = v_{k+1} / v_{k}. It's easy to write down r_k and s_k because they
        // follow immediately from the \Gamma function recurrence.
        //    r_{k} = \frac{(a + m + k)(b + m + k)}{(1 + k - e)(1 + m + k)}
        //    s_{k} = \frac{(a + m + k + e)(b + m + k + e)}{(1 + k)(1 + m + k + e)}
        // Notice that r_k and s_k are almost equal, but not quite: they differ by O(e). To advance t_k, use
        //    t_{k+1} = u_{k+1} + v_{k+1} = r_k u_k + s_k v_k = s_k (u_k + v_k) + (r_k - s_k) u_k
        //            = s_k t_k + d_k * u_k
        // where d_k = r_k - s_k, which will be O(e), since r_k and s_k only differ by e-shifted arguments.

        // In the x -> (1-x), x -> 1/x, x -> x / (1-x), and x -> 1 - 1/x linear transformations, canceling divergences
        // appear when some arguments of the transformed functions are non-positive-integers.

        private static double Hypergeometric2F1_Series_OneOverOneMinusX (double a, int m, double e, double c, double x1) {

            Debug.Assert(m >= 0);
            Debug.Assert(Math.Abs(e) <= 0.5);
            Debug.Assert(Math.Abs(x1) <= 0.75);

            double b = a + m + e;

            double g_c = AdvancedMath.Gamma(c);
            //double rg_a = 1.0 / AdvancedMath.Gamma(a);
            double rg_b = 1.0 / AdvancedMath.Gamma(b);
            double rg_cma = 1.0 / AdvancedMath.Gamma(c - a);
            //double rg_cmb = 1.0 / AdvancedMath.Gamma(c - b);

            // Pochhammer product, keeps track of (a)_k (c-b)_k (x')^{a + k}
            double p = Math.Pow(x1, a);

            double f0 = 0.0;
            if (m > 0) {

                f0 = p;

                double q = 1.0;
                for (int k = 1; k < m; k++) {
                    int km1 = k - 1;
                    p *= (a + km1) * (c - b + km1) * x1;
                    q *= (k - m - e) * k;
                    f0 += p / q; 
                }

                f0 *= g_c * rg_b * rg_cma * AdvancedMath.Gamma(m + e);
                p *= (a + (m - 1)) * (c - b + (m - 1)) * x1;
            }

            // Now compute the remaining terms with analytically canceled divergent parts.

            double t = rg_b * rg_cma * (NewG(1.0, -e) / AdvancedIntegerMath.Factorial(m) + NewG(m + 1, e)) -
                1.0 / AdvancedMath.Gamma(1 + m + e) * (NewG(a + m, e) / AdvancedMath.Gamma(c - a - e) + NewG(c - a, -e) / AdvancedMath.Gamma(b)) -
                MoreMath.ReducedExpMinusOne(Math.Log(x1), e) / AdvancedMath.Gamma(a + m) / AdvancedMath.Gamma(c - a - e) / AdvancedMath.Gamma(m + 1 + e);
            t *= p;

            double f1 = t;

            double u = p * rg_b * rg_cma / AdvancedMath.Gamma(1.0 - e) / AdvancedIntegerMath.Factorial(m);

            for (int k = 0; k < Global.SeriesMax; k++) {

                double f1_old = f1;

                int k1 = k + 1;
                int mk1 = m + k1;
                double amk = a + m + k;
                double amke = amk + e;
                double cak = c - a + k;
                double cake = cak - e;
                double k1e = k1 - e;
                double mk1e = mk1 + e;

                double r = amk * cake / k1e / mk1;
                double s = amke * cak / mk1e / k1;

                // Compute (r - s) / e analytically because leading terms cancel
                double d = (amk * cake / mk1 - amk - cake - e + amke * cak / k1) / mk1e / k1e;

                t = (s * t + d * u) * x1;

                f1 += t;

                if (f1 == f1_old) {
                    f1 *= ReciprocalSincPi(e) * g_c;
                    if (m % 2 != 0) f1 = -f1;
                    return (f0 + f1);
                }

                u *= r * x1;

            }

            throw new NonconvergenceException();

        }

        private static double Hypergeometric2F1_Series_OneMinusX (double a, double b, int m, double e, double x1) {

            Debug.Assert(m >= 0);
            Debug.Assert(Math.Abs(e) <= 0.5);
            Debug.Assert(Math.Abs(x1) <= 0.75);

            double c = a + b + m + e;

            // Compute all the gammas we will use.
            double g_c = AdvancedMath.Gamma(c);
            double rg_am = 1.0 / AdvancedMath.Gamma(a + m);
            double rg_bm = 1.0 / AdvancedMath.Gamma(b + m);
            double rg_ame = 1.0 / AdvancedMath.Gamma(a + m + e);
            double rg_bme = 1.0 / AdvancedMath.Gamma(b + m + e);
            double rg_m1e = 1.0 / AdvancedMath.Gamma(m + 1 + e);

            // Pochhammer product, keeps track of (a)_m (b)_m (x')^m
            double p = 1.0;

            // First compute the finite sum, which contains no divergent terms even for e = 0.
            double f0 = 0.0;
            if (m > 0) {

                double t0 = 1.0;
                f0 = t0;
                for (int k = 1; k < m; k++) {
                    int km1 = k - 1;
                    p *= (a + km1) * (b + km1) * x1;
                    t0 *= 1.0 / (1.0 - m - e + km1) / k;
                    f0 += t0 * p;
                }

                f0 *= g_c * rg_bme * rg_ame * AdvancedMath.Gamma(m + e);
                p *= (a + (m - 1)) * (b + (m - 1)) * x1;

            }

            // Now compute the remaining terms with analytically canceled divergent parts.

            double t = rg_bme * rg_ame * (NewG(1.0, -e) / AdvancedIntegerMath.Factorial(m) + NewG(m + 1, e)) -
                rg_m1e * (NewG(a + m, e) * rg_bme + NewG(b + m, e) * rg_am) -
                MoreMath.ReducedExpMinusOne(Math.Log(x1), e) * rg_am * rg_bm * rg_m1e;


            t *= p;
            double f1 = t;
            double u = p * rg_bme * rg_ame / AdvancedMath.Gamma(1.0 - e) / AdvancedIntegerMath.Factorial(m);

            for (int k = 0; k < Global.SeriesMax; k++) {

                double f1_old = f1;

                // Compute a bunch of sums we will use.
                int k1 = k + 1;
                int mk = m + k;
                int mk1 = mk + 1;
                double k1e = k1 - e;
                double amk = a + mk;
                double bmk = b + mk;
                double amke = amk + e;
                double bmke = bmk + e;
                double mk1e = mk1 + e;
 
                // Compute the ratios of each term. These are close, but not equal for e != 0.
                double r = amk * bmk / mk1 / k1e;
                double s = amke * bmke / mk1e / k1;

                // Compute (r - s) / e, with O(1) terms of (r - s) analytically canceled.
                double d = (amk * bmk / mk1 - (amk + bmk + e) + amke * bmke / k1) / mk1e / k1e;

                // Advance to the next term, including the correction for s != t.
                t = (s * t + d * u) * x1;

                f1 += t;

                if (f1 == f1_old) {
                    f1 *= ReciprocalSincPi(e) * g_c;
                    if (m % 2 != 0) f1 = -f1;
                    return (f0 + f1);
                }

                // Advance the u term, which we will need for the next iteration.
                u *= r * x1;

            }

            throw new NonconvergenceException();

        }

        private static double ReciprocalSincPi (double e) {
            Debug.Assert(Math.Abs(e) <= 1.0);
            if (e == 0.0) {
                return (1.0);
            } else {
                double x = Math.PI * e;
                return (x / Math.Sin(x));
            }
        }

    }
}
