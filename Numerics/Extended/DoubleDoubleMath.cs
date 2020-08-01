using System;
using System.Collections.Generic;
using System.Text;

namespace Meta.Numerics.Extended {

    // By analogy with Math, functions should go here. By analogy with BCL's Complex and BigInteger,
    // they should go on structure.

    internal static partial class DoubleDoubleMath {

        private static DoubleDouble Sin_Series (DoubleDouble x) {
            DoubleDouble mx2 = -x * x;
            DoubleDouble t = x;
            DoubleDouble s = t;
            for (int k = 3; k < Global.SeriesMax; k += 2) {
                DoubleDouble s_old = s;
                t *= mx2 / (k * (k - 1));
                s += t;
                if (s == s_old) {
                    return s;
                }
            }
            throw new NonconvergenceException();
        }

        private static DoubleDouble Cos_Series (DoubleDouble x) {
            DoubleDouble mx2 = -x * x;
            DoubleDouble t = 0.5 * mx2;
            DoubleDouble s = DoubleDouble.One + t;
            for (int k = 4; k < Global.SeriesMax; k += 2) {
                DoubleDouble s_old = s;
                t *= mx2 / (k * (k - 1));
                s += t;
                if (s == s_old) {
                    return s;
                }
            }
            throw new NonconvergenceException();
        }

        // Decompose x = \frac{\pi}{2} ( n + z1)
        // For now, this doesn't use extended precision value of pi, so
        // looses digits in z1 for higher x.

        private static void RangeReduction (DoubleDouble x, out long n, out DoubleDouble z1) {

            DoubleDouble z = x * TwoOverPi;
            double zHi = (double) z;
            double zInteger = Math.Round(zHi);

            n = (long) zInteger;
            z1 = z - zInteger;

        }

        private static readonly DoubleDouble PiOverTwo = 0.5 * DoubleDouble.Pi;

        private static readonly DoubleDouble TwoOverPi = 1.0 / PiOverTwo;

        /// <summary>
        /// Computes the sine of a double double real number.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of sin(x).</returns>
        /// <remarks>
        /// <para>The result of this function has essentially full precision for 
        /// arguments within the first few multiples of 2 pi, but slowly looses
        /// precision in higher digits as the argument becomes large.</para>
        /// </remarks>
        public static DoubleDouble Sin (DoubleDouble x) {

            RangeReduction(x, out long z0, out DoubleDouble z1);
            DoubleDouble x1 = z1 * PiOverTwo;

            switch (MoreMath.Mod(z0, 4L)) {
                case 0:
                    return Sin_Series(x1);
                case 1:
                    return Cos_Series(x1);
                case 2:
                    return -Sin_Series(x1);
                case 3:
                    return -Cos_Series(x1);
                default:
                    throw new InvalidOperationException();
            }

        }

        public static DoubleDouble Cos (DoubleDouble x) {

            RangeReduction(x, out long z0, out DoubleDouble z1);
            DoubleDouble x1 = z1 * PiOverTwo;

            switch (MoreMath.Mod(z0, 4L)) {
                case 0L:
                    return Cos_Series(x1);
                case 1L:
                    return -Sin_Series(x1);
                case 2L:
                    return -Cos_Series(x1);
                case 3L:
                    return Sin_Series(x1);
                default:
                    throw new InvalidOperationException();
            }

        }

    }
}
