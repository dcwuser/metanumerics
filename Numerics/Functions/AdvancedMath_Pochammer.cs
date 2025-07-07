using System;
using System.Diagnostics;

namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {

        /// <summary>
        /// Computes the logarithm of the the specified Pochammer symbol.
        /// </summary>
        /// <param name="x">The value of the base argument, which must be non-negative.</param>
        /// <param name="y">The value of the exponent argument, which must be greater than or equal to -<paramref name="x"/>.</param>
        /// <returns>The value of ln (x)<sub>y</sub>.</returns>
        /// <remarks>
        /// <para>While <see cref="Pochhammer(double, double)"/> does not overflow for much of the parameter space for which
        /// <see cref="Gamma(double)"/> of <paramref name="x"/> or (<paramref name="x"/> + <paramref name="y"/>) would, there
        /// are still areas of parameter space for which it does. Additionally, for small <paramref name="y"/>, <see cref="Pochhammer(double, double)"/>
        /// can be 1 to within floating-point accuracy, but you may still want to be able to quantify how much it differs from 1. This
        /// method overcomes both of those limitations by providing a non-overflowing value of the logarithm of the Pochammer symbol
        /// over its entire parameter space.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative or (<paramref name="x"/> + <paramref name="y"/>) is negative.</exception>
        public static double LogPochhammer (double x, double y) {

            double z = x + y;
            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));
            if (z < 0.0) throw new ArgumentOutOfRangeException(nameof(y));

            if (y == 0.0) {
                return 0.0;
            } else if (x == 0.0) {
                return Double.NegativeInfinity;
            } else if (y == 1.0) {
                return Math.Log(x);
            }

            // Deal with infinite arguments (and NaNs?)
            // if x is infinite is it also infinite?
            if (Double.IsPositiveInfinity(y)) return Double.PositiveInfinity;

            // The trick in all these cases is to find a form with good behavior (i.e. no numerical
            // cancellation) as y -> 0. I have done this for Stirling, Lanczos, and Series implementations
            // of Gamma.

            if (x < GammaSeries.Limit) {
                if ( z < GammaSeries.Limit) {
                    return GammaSeries.LogPochhammerTwoPlus(x, y) - MoreMath.LogOnePlus(y / x) - MoreMath.LogOnePlus(y / (1.0 + x));
                } else {
                    return LogGamma(z) - LogGamma(x);
                }
            } else if (x >= 16.0 && y >= 16.0) {
                return Stirling.LogPochhammer(x, y);
            } else if (Math.Abs(x - 1.0) < GammaSeries.Limit && Math.Abs(z - 1.0) < GammaSeries.Limit & z < 3.0 - x) {
                return GammaSeries.LogPochhammerTwoPlus(x - 1.0, y) - MoreMath.LogOnePlus(y / x);
            } else if (Math.Abs(x - 2.0) < GammaSeries.Limit && Math.Abs(z - 2.0) < GammaSeries.Limit) {
                return GammaSeries.LogPochhammerTwoPlus(x - 2.0, y);
            } else {
                return Lanczos.LogPochhammer(x, y);
            }
            /*
            if (x > 16.0 && z > 16.0) {
                return Stirling.LogPochhammer(x, y);
            } else if (x < 2.5) {

                // If x is small (near 1 or 2) we will need to use the series development around an integer value
                // to preserve accuracy.
                int n = (int)Math.Round(x);
                double dx = x - n;

                // If y is small enough, we can analytically take the difference of the series for
                // \ln\Gamma(x+y) and \ln\Gamma(x). This is necessary to perserve accuracy as y -> 0.
                // The precise allowed regions are overlapping parallelograms in dx, y space, centered
                // on each integer to develop the series around.
                // For a series radius of d_max, this is always possible for y > 2 d_max -1 (0.5 for d_max = 0.75).
                // Depending on dx, it is sometimes possible for y up do 2 d_max( 1.5 for d_max = 0.75).
                if (Math.Abs(y) < 1.5) {
                    if (Math.Abs(dx + y) < 0.75) return GammaSeries.LogPochhammer(n, dx, y);
                    double dxp = dx - 1.0;
                    if (Math.Abs(dxp) < 0.75 && Math.Abs(dxp + y) < 0.75) return GammaSeries.LogPochhammer(n + 1, dxp, y);
                    double dxm = dx + 1.0;
                    if (Math.Abs(dxm) < 0.75 && Math.Abs(dxm + y) < 0.75) return GammaSeries.LogPochhammer(n - 1, dxm, y);
                }
                // This code's choice of which n to develop around isn't optimal. In some cases,
                // we could get a smaller maximum deviation by picking a different integer. For example,
                // if dx = 0.4 and y = 0.3, we develop around n with d = 0.4 + 0.3 = 0.7, but
                // it would have been better to develop around n+1 with d = 1.0 - 0.4 = 0.6.
                // The logic for an optimal pick is complicated; consider doing it later.

                // This should only happen if |y| > 2d-1. Still x is small so we need to use the
                // series development of the denominator.
                Debug.Assert(Math.Abs(y) >= 0.5);
                return LogGamma(z) - GammaSeries.LogGamma(n, dx);
                //if (z < 2.5) {
                //    return GammaSeries.LogPochhammer(x, y);
                //} else {
                //    return LogGamma(z) - LogGamma(x);
                //}
            } else {
                return Lanczos.LogPochhammer(x, y);
            }
            */
        }

        // Note that for x < ~1.5 and y > 0 there is an additional line, at increasing y for decreasing x, at which (x)_y = 1 and there fore \ln (x)_y = 0
        // E.g. x = 0.5, y = ~2.365. This shows up via non-analytic cancellation between terms, and therefore we do loose accuracy along this line.

        // There are several transformation available for (x)_y = \frac{\Gamma(x + y)}{\Gamma(x)}.
        //   \Gamma(x) (x)_y = \Gamma(y) (y)_x
        // follows straightforwardly from the definition.
        //   (x)_y = 1 / (x + y)_{-y}
        // follows from writing x = x + y - y in denominator.
        // From application of -z \Gamma(-z) \Gamma(z) = \frac{\pi}{\sin(\pi z)} 

        /// <summary>
        /// Computes the Pochammer symbol (x)<sub>y</sub>.
        /// </summary>
        /// <param name="x">The first argument.</param>
        /// <param name="y">The second argument.</param>
        /// <returns>The value of (x)<sub>y</sub>.</returns>
        /// <remarks>
        /// <para>The Pochhammer symbol is defined as a ratio of Gamma functions.</para>
        /// <para>For positive integer y, this is equal to a rising factorial product.</para>
        /// <para>Note that while the Pochhammer symbol notation is very common, combinatorialists sometimes
        /// use the same notation for a falling factorial.</para>
        /// <para>If you need to compute a ratio of Gamma functions, in most cases you should compute it by calling this
        /// function instead of taking the quotient of two calls to the <see cref="Gamma(double)"/> function,
        /// because there is a large fraction of the parameter space where Gamma(x+y) and Gamma(x) overflow,
        /// but their quotient does not, and is correctly computed by this method.</para>
        /// <para>This function accepts both positive and negative values for both x and y.</para>
        /// </remarks>
        /// <seealso href="http://mathworld.wolfram.com/PochhammerSymbol.html"/>
        public static double Pochhammer(double x, double y) {

            // Deal with zero arguments first, both for efficiency and to avoid corner-cases later.
            if (y == 0.0) {
                return 1.0;
            } else if (y == 1.0) {
                return x;
            }
            Debug.Assert(y != 0.0);

            // To deal with x and z = x + y negative, we will need to transform to positive arguments.
            // The relevent regions are:
            //    x  y  z  transform
            //    +  +  +  none
            //    +  -  +  none
            //    +  -  -  transform \Gamma(z)
            //    -  +  +  transform \Gamma(x)
            //    -  +  -  transform \Gamma(z) and \Gamma(x) to (-x)_{-y}
            //    -  -  -  transform \Gamma(z) and \Gamma(x) to (-x)_{-y}
            double z = x + y;
            if (Double.IsNaN(z)) return z;

            if (x <= 0.0) {

                // If x is a negative integer then (x)_y is zero _unless_ y is also a negative integer.
                // This occurs beause x being a negative integer makes the Gamma(x) in the denominiator
                // blow up. But if x + y is a negative integer then Gamma(z) in the numerator also
                // blows up. What remains after cancelling infinities is (-1)^{-y} / (1-x)_{-y}
                double x0 = Math.Round(x);
                if (x0 == x) {
                    // Negative (or zero) integer x.
                    if (z <= 0.0) {
                        RangeReduction.ReduceByOnes(y, out long y0, out double y1);
                        if (y1 == 0.0) {
                            // (-m)_{-n} = (1)^n / (m+1)_n
                            return (y0 % 2 == 0 ? 1 : -1) / Pochhammer(1 - x0, -y0);
                            // Would be nice to replace this with integer Pochhammer when we have one
                        }
                    }

                    return 0.0;
                }

                if (z <= 0.0) {
                    // - - - -> + + +
                    // - + - -> + - +
                    return x / z * MoreMath.SinPi(x) / MoreMath.SinPi(z) / Pochhammer(-x, -y);
                } else {
                    // - + +
                    return -x * MoreMath.SinPi(x) / Math.PI * AdvancedMath.Gamma(-x) * AdvancedMath.Gamma(z);
                }
            } else if (z < 0.0) {
                // + - -
                return Math.PI / MoreMath.SinPi(z) / (-z * AdvancedMath.Gamma(x) * AdvancedMath.Gamma(-z));
            }
            // Note that it is also possible to transform - + - -> + + + by recognizing 1/(-x)_{-y} = (-z)_{y},
            // and similiarly + - + -> + + + by recognizing (x)_{y} = 1 / (z)_{-y}.
            // You might think this is nice because it allows us to always end up with + + +.
            // Note, however that it transforms a calculation with inputs x & y into a calculation with
            // inputs z & y, and the computation z = x + y introduces cancellation error in these cases
            // because x & y have opposite signs. I discovered this the hard way. We can work with
            // negative y as long as z is positive.

            // We have dealt with negative arguments and are now assured x and z are both positive.
            Debug.Assert(x >= 0.0);
            Debug.Assert(z >= 0.0);

            // If y is a small natural number, it will be faster to just do that many multiplications.
            if (0.0 < y && y < 32.0) {
                double yRound = Math.Round(y);
                if (y == yRound) {
                    double p = x;
                    for (int i = ((int)yRound) - 1; i > 0; i--) {
                        x += 1.0;
                        p *= x;
                    }
                    return p;
                }
            }


            if (x < GammaSeries.Limit) {
                if (z < GammaSeries.Limit) {
                    return MoreMath.SinPi(x) / MoreMath.SinPi(z) * (1.0 - z) / (1.0 - x) / GammaSeries.PochhammerTwoPlus(-x, -y);
                } else {
                    return MoreMath.SinPi(x) / Math.PI * GammaSeries.GammaOnePlus(-x) * AdvancedMath.Gamma(z);
                    // Amazingly, this isn't good enough, I guess because it uses computed z instead of y.
                }
            } else if (y > 172.0) {
                // As long as x isn't tiny, (x)_y overflows for all y greater than 172 (i.e. values for which Gamma overflows).
                // In these cases return early not just to be faster, but also to avoid NaN from intermediate 0 * infinity,
                // observed for example when Pochhammer(100,Infinity) went into Stirling code. 
                Debug.Assert(x > 0.5);
                return Double.PositiveInfinity;
            } else if (x >= 16.0 && z >= 16.0) {
                return Stirling.Pochhammer(x, y);
            } else if (Math.Abs(x - 1.0) < GammaSeries.Limit && Math.Abs(z - 1.0) < GammaSeries.Limit && z < 3.0 - x) {
                // Pochhammer (1 + x)_y
                // The additional test z < 3 - x is there because above that line the series around 2 will be closer.
                return x / z * GammaSeries.PochhammerTwoPlus(x - 1.0, y);
            } else if (Math.Abs(x - 2.0) < GammaSeries.Limit && Math.Abs(z - 2.0) < GammaSeries.Limit) {
                // Pochhammer (2 + x)_y
                return GammaSeries.PochhammerTwoPlus(x - 2.0, y);
            } else {
                // We are not in the series or Stirling regiemes, so fall back to Lanczos.
                // It it possible that one of the component gamma functions is in a regieme where we could compute
                // it more accurately, but the other one is not, so it wouldn't buy us much to do so.
                return Lanczos.Pochammer(x, y);
            }
            /*
            if (x > 16.0 && z > 16.0) {

                // If x and z are both large enough, Stirling is straightforward, fast, and accurate.
                // However it produces 0 * infinity = NaN for infinite y, so avoid calculation and
                // just return infinity if y is big enough that overflow is guaranteed.
                if (y < 1000.0) {
                    return Stirling.Pochhammer(x, y);
                } else {
                    return Double.PositiveInfinity;
                }
                // Keeping all leading order terms in Stirling approxmation, the criteria for overflow
                // would be (z - 1/2) ln z - (x - 1/2) ln x - y > 710. This is highly accurate, but at
                // the cost of two logs + 6 flops. 

            } else if (x < 2.5) {

                // If x is small, we want to use the series for the \Gamma(x) in the denominator to maintain accuracy

                // If y is also small enough, we can use series for \Gamma(x+y) too
                int n = (int)Math.Round(x);
                double dx = x - n;
                if (Math.Abs(y) < 1.5) {
                    if (Math.Abs(dx + y) < 0.75) return GammaSeries.Pochhammer(n, dx, y);
                    double dxp = dx - 1.0;
                    if (Math.Abs(dxp) < 0.75 && Math.Abs(dxp + y) < 0.75) return GammaSeries.Pochhammer(n + 1, dxp, y);
                    double dxm = dx + 1.0;
                    if (Math.Abs(dxm) < 0.75 && Math.Abs(dxm + y) < 0.75) return GammaSeries.Pochhammer(n - 1, dxm, y);
                }

                // Now we have x small, but y can still be almost anywhere
                Debug.Assert(Math.Abs(y) >= 0.5);

                if (z < 16.0) {
                    // If y is small enough, it's okay to just use \Gamma(x+y) / \Gamma(x) directly.
                    return Gamma(z) / GammaSeries.Gamma(n, dx);
                    // Neither numerator nor denominator will overflow for 0<x<2.5 and 0.5<y<16
                    // As y get bigger, though, we loose accuacry when computing z = x + y, so we don't want
                    // to use this for arbitrarily large y.
                } else if (z < 172.0) {
                    // For large y, we can use (x)_y \Gamma(x) = (y)_x \Gamma(y) to transform to
                    // an expression that uses y directly and is guaranteed to be in the Stirling regime.
                    return Stirling.Gamma(y) / GammaSeries.Gamma(n, dx) * Stirling.Pochhammer(y, x);
                    // Neither Gamma(x) nor (y)_x overflow, so no cancelling overflows here
                } else {
                    // For y this big, (x)_y is guaranteed to overflow
                    return (Double.PositiveInfinity);
                }

                // Fix this for y that would overflow Gamma(y)
                //if (y > 18.0) {
                //    return Stirling.Gamma(y) / GammaSeries.Gamma(n, dx) * Stirling.Pochhammer(y, x);
                //}

                //Debug.Assert(Math.Abs(y) >= 0.5);
                //if (y > 40.0) {
                //    return Lanczos.Pochammer(x, y);
                //} else {
                //    return Gamma(z) / GammaSeries.Gamma(n, dx);
                //}
                // This ratio can be inf/inf = NaN if n=0 and dx is sub-normal
                // This looses accuracy for y >> x. E.g. x = 0.44095141326363707 and y = 54 has error ~10^(-14)
                // even though Gamma(z) error is only ~10^(-15) and Gamma(x) error ~10^(-17). Problem is
                // presumably error introduced by forming z = x + y. I don't see any analytical cancellation
                // by using Stirling for Gamma(x+y) and series for Gamma(x).

            } else {
                // We have x > 2.5 and are not in the Stirling regime, so use Lanczos
                return Lanczos.Pochammer(x, y);
            }
            */
        }

    }
}
