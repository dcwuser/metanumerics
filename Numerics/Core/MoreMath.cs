using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

namespace Meta.Numerics {

    /// <summary>
    /// Contains addtional basic math operations.
    /// </summary>
    /// <remarks>
    /// <para>The <see cref="System.Math"/> class defines many basic math operations, but a few that are important for optimal numerical
    /// practice are missing. They are defined by this class.</para>
    /// </remarks>
    public static class MoreMath {

        /// <summary>
        /// Rasises an argument to an integer power.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <param name="n">The power.</param>
        /// <returns>The value of x<sup>n</sup>.</returns>
        /// <remarks>
        /// <para>Low integer powers can be computed by optimized algorithms much faster than the general
        /// alrogithm for an arbitrary real power employed by <see cref="System.Math.Pow"/>.</para>
        /// </remarks>
        public static double Pow (double x, int n) {

            if (n < 0) return (1.0 / Pow(x, -n));

            switch(n) {
                case 0:
                    // we follow convention that 0^0 = 1
                    return (1.0);
                case 1:
                    return (x);
                case 2:
                    // 1 multiply
                    return (x * x);
                case 3:
                    // 2 multiplies
                    return (x * x * x);
                case 4: {
                    // 2 multiplies
                    double x2 = x * x;
                    return (x2 * x2);
                }
                case 5: {
                    // 3 multiplies
                    double x2 = x * x;
                    return (x2 * x2 * x);
                }
                case 6: {
                    // 3 multiplies
                    double x2 = x * x;
                    return (x2 * x2 * x2);
                }
                case 7: {
                    // 4 multiplies
                    double x3 = x * x * x;
                    return (x3 * x3 * x);
                }
                case 8: {
                    // 3 multiplies
                    double x2 = x * x;
                    double x4 = x2 * x2;
                    return (x4 * x4);
                }
                case 9: {
                    // 4 multiplies
                    double x3 = x * x * x;
                    return (x3 * x3 * x3);
                }
                case 10: {
                    // 4 multiplies
                    double x2 = x * x;
                    double x4 = x2 * x2;
                    return (x4 * x4 * x2);
                }
                case 12: {
                    // 4 multiplies
                    double x3 = x * x * x;
                    double x6 = x3 * x3;
                    return (x6 * x6);
                }
                case 16: {
                    // 4 multiplies
                    double x2 = x * x;
                    double x4 = x2 * x2;
                    double x8 = x4 * x4;
                    return (x8 * x8);
                }
                // that's all the cases do-able in 4 or fewer multiplies
                default:
                    return (Math.Pow(x, n));
            }


        }

        // an internal method for squaring

        internal static double Pow2 (double x) {
            return (x * x);
        }

        /// <summary>
        /// Computes the length of a right triangle's hypotenuse.
        /// </summary>
        /// <param name="x">The length of one side.</param>
        /// <param name="y">The length of another side.</param>
        /// <returns>The length of the hypotenuse, sqrt(x<sup>2</sup> + y<sup>2</sup>).</returns>
        /// <remarks>
        /// <para>The length is computed accurately, even in cases where
        /// x<sup>2</sup> or y<sup>2</sup> would overflow.</para>
        /// </remarks>
        public static double Hypot (double x, double y) {
            if ((x == 0.0) && (y == 0.0)) {
                return (0.0);
            } else {
                double ax = Math.Abs(x);
                double ay = Math.Abs(y);
                if (ax > ay) {
                    double r = y / x;
                    return (ax * Math.Sqrt(1.0 + r * r));
                } else {
                    double r = x / y;
                    return (ay * Math.Sqrt(1.0 + r * r));
                }
            }
        }

        /// <summary>
        /// Computes e<sup>x</sup>-1.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of e<sup>x</sup>-1.</returns>
        /// <remarks>
        /// <para>If x is close to 0, then e<sup>x</sup> is close to 1, and computing e<sup>x</sup>-1 by by subtracting one from
        /// e<sup>x</sup> as computed by the <see cref="Math.Exp"/> function will be subject to severe loss of significance due to
        /// cancelation. This method maintains full precision for all values of x by switching to a series expansion for values of
        /// x near zero.</para>
        /// </remarks>
        public static double ExpMinusOne (double x) {
            if (Math.Abs(x) < 0.125) {
                // For small x, use the series e^{x} = \sum_{k=0}^{\infty} \frac{x^k}{k!} = 1 + x + x^2 / 2! + x^3 / 3! + \cdots
                double dr = x;
                double r = dr;
                for (int k = 2; k < Global.SeriesMax; k++) {
                    double r_old = r;
                    dr *= x / k;
                    r += dr;
                    if (r == r_old) return (r);
                }
                throw new NonconvergenceException();
            } else {
                return (Math.Exp(x) - 1.0);
            }
        }

        /// <summary>
        /// Computes log(1+x).
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of log(1+x).</returns>
        /// <remarks>
        /// <para>If x is close to 0, computing log(1+x) by first adding one and then taking the log can result in a loss of accuracy.
        /// This function maintains full precision of all values of x by switching to a series expansion for values of x near zero.</para>
        /// </remarks>
        public static double LogOnePlus (double x) {
            if (Math.Abs(x) < 0.125) {
                // For small x, use the series \log(1-x) = - \sum_{k=1}^{\infty} \frac{x^k}{k}. 
                double xk = x;
                double f = xk;
                for (int k = 2; k < Global.SeriesMax; k++) {
                    double f_old = f;
                    xk *= -x;
                    f += xk / k;
                    if (f == f_old) return (f);
                }
                throw new NonconvergenceException();
            } else {
                return (Math.Log(1.0 + x));
            }
        }

        /// <summary>
        /// Computes x<sup>2</sup>.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The square of the argument.</returns>
        /// <remarks>
        /// <para>There is nothing numericaly sophisticated inside this function; it is simply for programmers' convenience. Given a complicated expression
        /// that needs to be squared, it is nice to be able to wrap it in a simple function call instead of explicitly assigning it to a new variable and then
        /// multiplying that variable by itself.</para>
        /// </remarks>
        public static double Sqr (double x) {
            return (x * x);
        }

        /// <summary>
        /// Computes the sine of the given argument to full significance over the full range of 
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of sin(x).</returns>
        /// <remarks>
        /// <para>This method addresses several subtle shortcommings of the trigonometric functions of the <see cref="System.Math"/> class.
        /// One shortcoming, quite striking but rarely encountered, is that <see cref="Math.Sin"/> returns entirely wrong results very large
        /// arguments; for x larger than about 10<sup>20</sup>, it simply returns the argument as the function value. Another shortcomming,
        /// more commonly encountered but often unnoticed, is that values are not precise near zeros of the function.</para>
        /// <para>The underlying issue is that 
        /// (I have no idea
        /// why the base class library designers did not at least choose to return <see cref="Double.NaN"/> so as to signal to the user
        /// that the result should not be trusted.) 
        /// </para>
        /// <para>
        /// This method 
        /// </para>
        /// <para>For typical arguments, say between 10<sup>-4</sup>-4 and 10<sup>4</sup>, the extra cost of this function over
        /// <see cref="Math.Sin"/> is just a couple of comparisons and a single floating point operation; less than 0.1% of
        /// arguments in this range are then routed to our much slower, higher-accuracy algorithm. We therefore suggest that,
        /// for general use, you prefer this method over <see cref="Math.Sin"/>; only in very unusual situations where (i) you are guaranteed
        /// never to encounter very large arguments, (ii) full precision values are not required, and (iii) the run-time of your
        /// application is critical and dominated by trigonometric calculations, should you prefer the base class library method.</para>
        /// </remarks>
        public static double Sin (double x) {
            if (x < 0.0) {
                return (-Sin(-x));
            } else if (x < 1.0) {
                // If x is small enough not to cross a zero, use the built-in function
                return (Math.Sin(x));
            } else if (x < RangeReduction.xLimit) {
                // If x is in the intermediate region, try the built-in function but switch to our range-reduction
                // algorithm if the result is too small.
                double y = Math.Sin(x);
                double ym = x / RangeReduction.xLimit;
                if (Math.Abs(y) < ym) {
                    return (RangeReduction.Sin(x));
                } else {
                    return (y);
                }
            } else {
                // If x is beyond the range of the built-in function eintirely, use our range-reduction algorithm
                return (RangeReduction.Sin(x));
            }
        }


        public static double Cos (double x) {

            // Ensure x is non-negative; this is easy because cosine is an even function.
            if (x < 0.0) x = - x;

            if (x < 1.0) {
                // If x is small enough not to cross a zero, use the built-in function
                return (Math.Cos(x));
            } else if (x < RangeReduction.xLimit) {
                // If x is in the intermediate region, try the built-in function but switch to our range-reduction
                // algorithm if the result is too small.
                double y = Math.Cos(x);
                double ym = x / RangeReduction.xLimit;
                if (Math.Abs(y) < ym) {
                    return (RangeReduction.Cos(x));
                } else {
                    return (y);
                }
            } else {
                // If x is beyond the range of the built-in function eintirely, use our range-reduction algorithm
                return (RangeReduction.Cos(x));
            }
        }

    }


    internal static class RangeReduction {

        // The number of bits to keep in each part when doing Veltkamp splitting

        private const int bitsPerPart = 26;

        // The corresponding constant used in the Veltkamp splitting algorithm

        private const double K = (1 << bitsPerPart) + 1;

        // The binary represenation of 2/\pi. We will use this to find z = (2\pi) x.

        // The largest double value is 2^1024 = 1.7E308. Since we need to get around
        // 54 bits after the decimal point in order to get full precision, we need
        // 2/\pi to 1024+54 = 1078 bits, at least.

        private const string twoOverPi =
            "101000101111100110000011011011100100111001000100000101010010100111111100001001110" +
            "1010111110100011111010100110100110111011100000011011011011000101001010110011" +
            "0010011110001000011100100000100000111111110010100010110001110101011110111101" +
            "0111011110001010110000110110111001001000110111000111010010000100100110111010" +
            "0101110000000000110010010010010111011101010000010011101000110010010000111001" +
            "1111110000111011110101100011100101100010010100110100111001111101110100010000" +
            "0100011010111110101001011101011101101000100100001001110100110011100011100000" +
            "0100110101101000101111101111110010000010011100110010001110101100011100110000" +
            "0110101001100111001111101001001110010000100010111111000101110111101111110010" +
            "0101000001110110001111111111000100101111111111111011110000001011001100000001" +
            "1111110111100101111000100011000101101011010000010100110110100011111011011010" +
            "0110110011111101100111100100111110010110000100110110111010011110100011000111" +
            "1110110011010011110010111111110101000101101011101010010011110111010110001111" +
            "1101011111001011111000101111011001111010000011100111001111101111000101001010" +
            "010100100101110101001101100";

        // go a bit longer

        private static double[] ParseBinaryString (string text) {

            List<double> parts = new List<double>();

            double value = 1.0;
            double part = 0.0;
            int bits = 0;

            for (int i = 0; i < text.Length; i++) {
                value = value / 2.0;
                if (text[i] == '1') {
                    part += value;
                }
                if ((part != 0.0) || (parts.Count > 0)) bits++;
                if (bits == bitsPerPart) {
                    parts.Add(part);
                    part = 0.0;
                    bits = 0;
                }
            }
            if (part != 0.0) parts.Add(part);

            return (parts.ToArray());

        }

        private static readonly double[] twoOverPiParts = ParseBinaryString(twoOverPi);

        public static void Add (double a, double b, out double s, out double r) {
            if (Math.Abs(a) < Math.Abs(b)) Global.Swap(ref a, ref b);
            s = a + b;
            double z = s - a;
            r = b - z;
        }

        public static double[] Decompose (double x) {
            double hi, lo;
            Decompose(x, out hi, out lo);
            return (new double[] { hi, lo });
        }

        public static void Decompose (double x, out double hi, out double lo) {
            double p = K * x;
            double q = x - p;
            hi = p + q;
            lo = x - hi;
        }

        // This is the Dekker multiplication algorithm.

        // Multiply x by y. x is a double split into high and low parts xHi and xLo. y is a high-accuracy float whoose
        // parts are stored in yParts. The output is z0 + z1, where z0 is an integer and -0.5 < z1 < 0.5.

        public static void Multiply (double xHi, double xLo, double[] yParts, out long z0, out double z1) {

            z0 = 0L;
            z1 = 0.0;

            for (int i = 0; i < yParts.Length; i++) {

                long z0_old = z0;
                double z1_old = z1;

                double t = xHi * yParts[i];
                if (i > 0) t += xLo * yParts[i - 1];

                double tr = Math.Round(t);
                if (tr > K * K) continue;
                z0 += (long) tr;
                double dz = t - tr;
                z1 += dz;

                if (Math.Abs(z1) > 0.5) {
                    double zr = Math.Round(z1);
                    z0 += (long) zr;
                    z1 -= zr;
                }

                if ((z0 == z0_old) && (z1 == z1_old)) return;

            }

            throw new NonconvergenceException();

        }


        public static void Multiply (double[] a, double[] b, out long z0, out double z1) {

            z0 = 0L;
            z1 = 0.0;

            for (int i = 0; i < a.Length + b.Length - 1; i++) {

                long z0_old = z0; double z1_old = z1;

                double ab = 0.0; 
                for (int j = Math.Max(0, i - b.Length + 1); j < Math.Min(a.Length, i + 1); j++) {
                    ab += a[j] * b[i - j];
                }

                double zr = Math.Round(ab);

                if (Math.Abs(zr) > Int64.MaxValue / 4) continue;

                z0 += (long) zr;
                double dz = ab - zr;
                z1 += dz;

                if ((z0 == z0_old) && (z1 == z1_old)) return;
            }

            throw new NonconvergenceException();

        }

        private static long Mod (long n, long m) {
            long r = n % m;
            if (r < 0) r += m;
            return (r);
        }

        public static double Sin (long z0, double z1) {
            switch (Mod(z0, 4L)) {
                case 0L:
                    return (Math.Sin(z1 * Math.PI / 2.0));
                case 1L:
                    return (Math.Cos(z1 * Math.PI / 2.0));
                case 2L:
                    return (-Math.Sin(z1 * Math.PI / 2.0));
                case 3L:
                    return (-Math.Cos(z1 * Math.PI / 2.0));
                default:
                    throw new InvalidOperationException();
            }

        }

        public static double Cos (long z0, double z1) {
            switch (Mod(z0, 4L)) {
                case 0L:
                    return (Math.Cos(z1 * Math.PI / 2.0));
                case 1L:
                    return (-Math.Sin(z1 * Math.PI / 2.0));
                case 2L:
                    return (-Math.Cos(z1 * Math.PI / 2.0));
                case 3L:
                    return (Math.Sin(z1 * Math.PI / 2.0));
                default:
                    throw new InvalidOperationException();
            }
        }

        public static double Sin (double x) {
            long z0; double z1;
            ReduceByPiHalves(x, out z0, out z1);
            return (Sin(z0, z1));
        }

        public static double Cos (double x) {
            long z0; double z1;
            ReduceByPiHalves(x, out z0, out z1);
            return (Cos(z0, z1));
        }

        public static void ReduceByPiHalves (double x, out long x0, out double x1) {

            PrintFrExp(x);
            double[] xParts = Decompose(x);

            for (int i = 0; i < xParts.Length; i++) PrintFrExp(xParts[i]);
            for (int i = 0; i < 4; i++) PrintFrExp(twoOverPiParts[i]);

            Multiply(xParts, twoOverPiParts, out x0, out x1);
        }

        public static void ReduceByOnes (double y, out long y0, out double y1) {
            double yr = Math.Round(y);
            if (Math.Abs(yr) < intMax) {
                y0 = (long) yr;
            } else {
                y0 = Math.Sign(yr) * intMax;
            }
            y1 = y - yr;
        }

        // There are 54 bits of accuracy in a double, so any double larger than 2^54 is effectively an integer,
        // and any double larger than 2^56 is effectively an integer divisible by 4.
        // Since there are 63 bits of accuracy in a long, the integer part of any double lower than this can
        // fit in a long. Any double larger than this value we simply peg at this value.
        private const long intMax = (1L << 56);

        public static double Cos (double x, double y) {
            long x0; double x1;
            ReduceByPiHalves(x, out x0, out x1);

            long y0; double y1;
            ReduceByOnes(y, out y0, out y1);

            long z0 = x0 + y0;
            double z1 = x1 + y1;

            return (Cos(z0, z1));

        }

        // Extensive comparisons of the result of our range-reduction algorithm to the result of the built-in algorithm
        // indicates that loss of significant figures tends to occur near zeros. We have found y-values above which
        // agreement between the built-in and our algorithm is always better than 1 part in 10^{14}.
        //  1 < x < 2      requires     |y| < 2^{-22}       search granularity 2^{-24}
        //  2 < x < 4                   |y| < 2^{-20}
        //  4 < x < 8                   |y| < 2^{-18}
        //  8 < x < 16                  |y| < 2^{-18}
        //  16 < x < 32                 |y| < 2^{-16}
        //  32 < x < 64                 |y| < 2^{-16}       search granularity 2^{-20}
        //  64 < x < 128                |y| < 2^{-14}
        //  128 < x < 256               |y| < 2^{-14}
        //  256 < x < 512               |y| < 2^{-12}
        //  512 < x < 1024              |y| < 2^{-12}
        //  1024 < x < 2048             |y| < 2^{-10}       search granularity 2^{-16}
        //  2048 < x < 4096             |y| < 2^{-10}
        //  4096 < x < 8192             |y| < 2^{-8}
        //  8192 < x < 16384            |y| < 2^{-8}
        //  16384 < x < 32768           |y| < 2^{-6}
        // Pattern of significance loss appears quite consistent: we need to drop one power of 2 from our y-limit for every
        // power of 2 added to x. In summary, yMax = x / 2^{20}.

        public const double xLimit = (1 << 20);

        public static void FrExp (double x, out bool negative, out long mantissa, out int exponent) {

            long bits = BitConverter.ToInt64(BitConverter.GetBytes(x), 0);
            //long bits = BitConverter.DoubleToInt64Bits(x);

            negative = (bits < 0);
            exponent = (int) ((bits >> 52) & 0x7ffL);
            mantissa = bits & 0xfffffffffffffL;

            exponent -= 1075;

            mantissa = mantissa | (1L << 52);

            while ((mantissa & 1) == 0) {
                mantissa >>= 1;
                exponent++;
            }
        }

        public static void PrintFrExp (double x) {
            bool negative; long mantissa; int exponent;
            FrExp(x, out negative, out mantissa, out exponent);
            Debug.WriteLine("{0:R} = {1}{2} X 2^({3})", x, negative ? "-" : "+", mantissa, exponent);
        }

    }


}