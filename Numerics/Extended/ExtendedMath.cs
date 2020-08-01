using System;
using System.Diagnostics;

namespace Meta.Numerics.Extended {

    // http://web.mit.edu/tabbott/Public/quaddouble-debian/qd-2.3.4-old/docs/qd.pdf

    internal static class ExtendedMath {

        public static bool IsNotFinite (double x) {
            ulong storage = (ulong) BitConverter.DoubleToInt64Bits(x);
            return ((storage & 0x7ff0000000000000) == 0x7ff0000000000000);
        }

        public static void FastTwoSum (double a, double b, out double sum, out double err) {
            // The following diagram illustrates how this algorithm recovers the digits of b that are lost
            // when computing (a+b):
            //         a = aaaa
            //         b =   bbbb
            // s = a + b = ssss
            // u = s - a =   bb00
            // e = b - u =     bb

            Debug.Assert(Math.Abs(a) >= Math.Abs(b) || Double.IsNaN(a) || Double.IsNaN(b));

            sum = a + b;
            double u = sum - a;
            err = b - u;
        }

        public static void TwoSum (double a, double b, out double sum, out double err) {
            if (Math.Abs(a) >= Math.Abs(b)) {
                FastTwoSum(a, b, out sum, out err);
            } else {
                FastTwoSum(b, a, out sum, out err);
            }
        }

        public static void ThreeSum (double a, double b, double c, out double sum, out double err) {
            double t1, t2, t3;
            ExtendedMath.TwoSum(a, b, out t1, out t2);
            ExtendedMath.TwoSum(c, t1, out sum, out t3);
            ExtendedMath.TwoSum(t2, t3, out err, out t1);
        }

        // 2^27 + 1 is the constant to split a 52-bit mantissa into two 26-bit parts
        private const double c = (double) ((1 << 27) + 1);

        public static void Split (double x, out double hi, out double lo) {
            double t = c * x;
            double u = x - t;
            hi = t + u;
            lo = x - hi;
        }

        public static void TwoProduct (double a, double b, out double product, out double err) {
            // Multiplying two n-digit numbers can produce a numer with up to 2n digits.
            // If the numbers are 52-bit doubles, the the product contains up to 104 bits,
            // but we loose the bottom 52. But if the numbers are just 26 bits but stored
            // in a 52-bit double container, then we obtain the exact product.

            // If we use Veltkamp splitting to split inputs a and b into hi and low
            // parts that each use only half the available bits of a double, then
            // products of those parts are exact and can be used to obtain an exact
            // product.

            // a * b = (aHi + aLo) * (bHi + bLo) = aHi * bHi + aHi * bLo + aLo * bHi + aLo * bLo
            //       = (aHi * bHi) + (aHi *bLo + aLo * bHi + aLo * bLo)

            // aaaa * bbbb = (pp00 + 00qq) * (rr00 + 00ss)
            //             = pp00 * rr00 + (pp00 * 00ss + 00qq * rr00) + (00qq * 00ss)
            //             = uuuu0000 + 00vvvv00 + 0000wwww
            //              = (uuuu + 00vv) + (vv00 + wwww)
            // The first part can be obtained 

            product = a * b;

            double aHi, aLo;
            Split(a, out aHi, out aLo);

            double bHi, bLo;
            Split(b, out bHi, out bLo);

            err = ((aHi * bHi - product) + aHi * bLo + aLo * bHi) + aLo * bLo;
        }

    }

}
