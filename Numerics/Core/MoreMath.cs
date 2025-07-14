using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics {

    /// <summary>
    /// Contains additional basic math operations.
    /// </summary>
    /// <remarks>
    /// <para>The <see cref="System.Math"/> class defines many basic math operations, but a few that are important for optimal numerical
    /// practice are missing. They are defined by this class.</para>
    /// </remarks>
    public static class MoreMath {

        /// <summary>
        /// Raises an argument to an integer power.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <param name="n">The power.</param>
        /// <returns>The value of x<sup>n</sup>.</returns>
        /// <remarks>
        /// <para>Low integer powers can be computed by optimized algorithms much faster than the general
        /// algorithm for an arbitrary real power employed by <see cref="System.Math.Pow"/>.</para>
        /// </remarks>
        public static double Pow (double x, int n) {

            if (n < 0) return 1.0 / Pow(x, -n);

            switch (n) {
                case 0:
                    // We follow convention that 0^0 = 1.
                    return 1.0;
                case 1:
                    return x;
                case 2:
                    // 1 multiply
                    return x * x;
                case 3:
                    // 2 multiplies
                    return x * x * x;
                case 4: {
                        // 2 multiplies
                        double x2 = x * x;
                        return x2 * x2;
                    }
                case 5: {
                        // 3 multiplies
                        double x2 = x * x;
                        return x2 * x2 * x;
                    }
                case 6: {
                        // 3 multiplies
                        double x2 = x * x;
                        return x2 * x2 * x2;
                    }
                case 7: {
                        // 4 multiplies
                        double x3 = x * x * x;
                        return x3 * x3 * x;
                    }
                case 8: {
                        // 3 multiplies
                        double x2 = x * x;
                        double x4 = x2 * x2;
                        return x4 * x4;
                    }
                case 9: {
                        // 4 multiplies
                        double x3 = x * x * x;
                        return x3 * x3 * x3;
                    }
                case 10: {
                        // 4 multiplies
                        double x2 = x * x;
                        double x4 = x2 * x2;
                        return x4 * x4 * x2;
                    }
                case 12: {
                        // 4 multiplies
                        double x3 = x * x * x;
                        double x6 = x3 * x3;
                        return x6 * x6;
                    }
                case 16: {
                        // 4 multiplies
                        double x2 = x * x;
                        double x4 = x2 * x2;
                        double x8 = x4 * x4;
                        return x8 * x8;
                    }
                // Those are all the cases do-able in 4 or fewer multiplies.
                default:
                    return Math.Pow(x, n);
            }

            // I verified that this routine is measurably faster than Math.Pow for the
            // cases that require 4 or fewer hand-coded multiplies.

            // I also tried doing higher powers via exponentiation-by-squaring, but found
            // that was actually slightly slower than Math.Pow.

        }

        /// <summary>
        /// Computes the length of a right triangle's hypotenuse.
        /// </summary>
        /// <param name="x">The length of one side.</param>
        /// <param name="y">The length of another side.</param>
        /// <returns>The value of sqrt(x<sup>2</sup> + y<sup>2</sup>).</returns>
        /// <remarks>
        /// <para>The length is computed accurately, even in cases where
        /// x<sup>2</sup> or y<sup>2</sup> would overflow or underflow.</para>
        /// </remarks>
        public static double Hypot (double x, double y) {

            double ax = Math.Abs(x);
            double ay = Math.Abs(y);

            double small, big;
            if (ax < ay) {
                small = ax;
                big = ay;
            } else {
                small = ay;
                big = ax;
            }

            if (small == 0.0) {
                return big;
            } else if (Double.IsPositiveInfinity(big) && !Double.IsNaN(small)) {
                return Double.PositiveInfinity;
            } else {
                double ratio = small / big;
                return big * Math.Sqrt(1.0 + ratio * ratio);
            }

        }

        // Beebe, "Computation of expm1(x) = exp(x)  - 1", 2002 (http://www.math.utah.edu/~beebe/reports/expm1.pdf)
        // makes some good points about e^x - 1.
        //   * He shows that the points e^x = 1/2 and e^x = 3/2 are the relevant limits where Math.Exp(x) - 1.0
        //     looses one bit of accuracy.
        //   * He measures that the maximum number of terms in the Taylor series required in this region is 17.
        //   * He measures that the RMS error of the Taylor series in this region is ~0.8 bits and it's maximum
        //     relative error is ~ 2.7 bits.
        //   * He points out that the doubling formula expm1(x) = expm1(x/2) * (expm1(x/2) + 2) allows one
        //     to reduce this by 2 terms and save 3 flops. But this is hardly worth the complication
        //     and he finds that it actually looses a smidgeon of accuracy.
        //   * He reviews several much more complicated schemes, e.g. minimax rational approximations,
        //     which allow more significant efficiency gains and error reduction.
        // For now, I use the Taylor series with his limits.

        private const double expm1SeriesLowerLimit = -0.693147180;
        private const double expm1SeriesUpperLimit = 0.4054651081;

        // expm1(x) - 1 = 1 + x + x^2 / 2 + x^3 / 3! + \cdots - 1
        //              = x + x^2 / 2 + x^3 / 3! + \cdots
        //              = x (1 + x / 2 + x^2 / 3! + \cdots)

        private static double ReducedExpm1Series (double x) {
            double df = 0.5 * x;
            double f = 1.0 + df;
            for (int k = 3; k < 20; k++) {
                double f_old = f;
                df *= x / k;
                f += df;
                if (f == f_old) return f;
            }
            throw new NonconvergenceException();
        }

        /// <summary>
        /// Computes e<sup>x</sup>-1.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of e<sup>x</sup>-1.</returns>
        /// <remarks>
        /// <para>If x is close to 0, then e<sup>x</sup> is close to 1, and computing e<sup>x</sup>-1 by
        /// <c>Math.Exp(x) - 1.0</c> will be subject to severe loss of significance due to cancelation.
        /// This method maintains full precision for all values of x by switching to a series expansion
        /// when x is near zero.</para>
        /// </remarks>
        public static double ExpMinusOne (double x) {
            if ((expm1SeriesLowerLimit < x) && (x < expm1SeriesUpperLimit)) {
                return x * ReducedExpm1Series(x);
            } else {
                return (Math.Exp(x) - 1.0);
            }
        }

        internal static double ReducedExpMinusOne (double x) {
            if ((expm1SeriesLowerLimit < x) && (x < expm1SeriesUpperLimit)) {
                return ReducedExpm1Series(x);
            } else {
                return (Math.Exp(x) - 1.0) / x;
            }
        }

        // Computes [ e^(\epsilon * x) - 1 ] / \epsilon, i.e. the rate at which ExpMinusOne grows with x

        internal static double ReducedExpMinusOne (double x, double e) {
            if (x == Double.NegativeInfinity) return (0.0);
            double y = e * x;
            if (Math.Abs(y) < 0.125) {
                // For small x, use the series e^{x} = \sum_{k=0}^{\infty} \frac{x^k}{k!} = 1 + x + x^2 / 2! + x^3 / 3! + \cdots
                double dr = x;
                double r = dr;
                for (int k = 2; k < Global.SeriesMax; k++) {
                    double r_old = r;
                    dr *= y / k;
                    r += dr;
                    if (r == r_old) return (r);
                }
                throw new NonconvergenceException();
            } else {
                return ((Math.Exp(y) - 1.0) / e);
            }


        }

        // Theorem 4 of Goldberg's classic "What Every Computer Scientist Should Know About Floating-Point Arithmetic"
        // shows that, if log is 1/2-ulp accurate and arithmetic is performed with a guard digit and the compiler
        // respects parenthesis even when optimizing, then
        //   log1p(x) = x log(1 + x) / ((1 + x) - 1)
        // is accurate to 5 ulp.

        // It looks like the GO implementation of log1p uses a more complex algorithm involving polynomial fits
        // that they claim is accurate to 1 ulp, so that might be good to look into.

        // Previously, I used the series development, but not over a wide enough range of arguments. So this trick
        // is likely better and faster than my previous implementation.

        /// <summary>
        /// Computes log(1+x).
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of log(1+x).</returns>
        /// <remarks>
        /// <para>If x is close to 0, computing log(1+x) by <c>Math.Log(1.0 + x)</c> can result in a
        /// significant loss of accuracy.
        /// This function maintains full precision of all values of x by switching to a series expansion
        /// for values of x near zero.</para>
        /// </remarks>
        public static double LogOnePlus (double x) {
            double z = 1.0 + x;
            if (z == 1.0) {
                return x;
            } else {
                return (Math.Log(z) / (z - 1.0) * x);
            }
        }

        // Computes \log (1 + e * x) / e, i.e. the rate at which LogOnePlus changes with e

        internal static double ReducedLogOnePlus (double x) {
            double z = 1.0 + x;
            if (z == 1.0) {
                return z;
            } else {
                return Math.Log(z) / (z - 1.0);
            }
        }

        internal static double ReducedLogOnePlus (double x, double e) {
            double y = e * x;
            if (Math.Abs(y) < 0.125) {
                // For small x, use the series \log(1-x) = - \sum_{k=1}^{\infty} \frac{x^k}{k}. 
                double xk = x;
                double f = xk;
                for (int k = 2; k < Global.SeriesMax; k++) {
                    double f_old = f;
                    xk *= -y;
                    f += xk / k;
                    if (f == f_old) return (f);
                }
                throw new NonconvergenceException();
            } else {
                return (Math.Log(1.0 + y) / e);
            }

        }

        /// <summary>
        /// Computes x<sup>2</sup>.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The square of the argument.</returns>
        /// <remarks>
        /// <para>There is nothing numerically sophisticated inside this function; it exists simply for programmers' convenience. Given a complicated expression
        /// that needs to be squared, it is nice to be able to wrap it in a simple call to this function instead of explicitly assigning its value to a new variable and then,
        /// in a separate statement, multiplying that variable by itself. Even if you are hyper-vigilant about function call overhead, you should not worry
        /// about using this function, because even the most basic optimizing compiler will optimize away the call.</para>
        /// </remarks>
        public static double Sqr (double x) {
            return x * x;
        }

        /// <summary>
        /// The conversion factor from degrees to radians.
        /// </summary>
        /// <remarks>
        /// <para>This conversion factor makes it easier to compute trigonometric functions if arguments are
        /// given in degrees.
        /// Since trigonometric methods such as <see cref="Math.Sin"/> and <see cref="MoreMath.Cos"/> take arguments
        /// in dimensionless radians, you must convert a degree input to radians before passing it into one of these
        /// functions. This field makes it easy to do so via a simple and visually mnemonic multiplication. If x is
        /// in degrees and you wish to take its sine, just write: Math.Sin(x * MoreMath.Degrees).</para>
        /// </remarks>
        public static readonly double Degrees = Math.PI / 180.0;

        /// <summary>
        /// Computes the sine of the given value to full significance over the full range of arguments.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of sin(x).</returns>
        /// <remarks>
        /// <para>This method addresses several subtle shortcomings of the <see cref="System.Math.Sin" /> method.
        /// One shortcoming, quite striking but rarely encountered, is that <see cref="Math.Sin"/> returns entirely wrong results very large
        /// arguments -- for x larger than about 10<sup>20</sup>, it simply returns the argument as the function value!
        /// (I have no idea
        /// why the base class library designers did not at least choose to return <see cref="Double.NaN"/>, so as to signal to the caller
        /// that the result should not be trusted. No floating point standard specifies this crazy behavior.)
        /// Another shortcoming, more commonly encountered but often unnoticed, is that for large but not necessarily very large arguments,
        /// function values loose precision, particularly near zeros of the function.</para>
        /// <para>
        /// One way to view these shortcomings is that they are justified by the uncertainty inherent in floating point representations. In this
        /// view, any <see cref="System.Double"/> should be seen as an uncertain value with a relative error of ~10<sup>-16</sup>. If the
        /// argument is large enough, then the absolute size of this error can be as large or larger than 2&#x3C0;; in this circumstance we should not
        /// expect to be able to say anything about the value (except, of course, that it is between -1 and +1, which is violated by the designers' crazy
        /// choice to return the argument as the value). Even if the absolute error is just a non-negligible faction of 2&#x3C0;, there is a non-negligible fraction
        /// of the values between -1 and +1 in the corresponding range of function values; any of these values is as possible as any other as a value
        /// for the sine of our uncertain argument, so we should be satisfied with any returned value in this non-negligible range.
        /// </para>
        /// <para>A different view is that it is better to regard every representable floating point value as some exact rational number, and
        /// when computing functions of floating point numbers, we should strive to return the representable floating point value nearest
        /// to the actual function value for that exact rational.
        /// Callers are unlikely to complain if we are careful in this regard, and this behavior is particularly
        /// useful when the argument is an intermediate result that the programmer may not even realize has become large.
        /// Thus is the view that we adopt, and therefore we provide this improved trigonometric function.</para>
        /// <para>For typical arguments, say between -10<sup>4</sup> and 10<sup>4</sup>, the extra cost of calling this function instead of
        /// <see cref="Math.Sin"/> is just a couple of comparisons and a single floating point operation; less than 0.1% of
        /// arguments in this range are then routed to our much slower, high-accuracy algorithm. We therefore suggest that,
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
            } else if (x < Double.PositiveInfinity) {
                // If x is beyond the range of the built-in function entirely, use our range-reduction algorithm
                return (RangeReduction.Sin(x));
            } else {
                // If x is infinite or NaN, return NaN.
                return (Double.NaN);
            }
        }

        /// <summary>
        /// Computes the cosine of the given value to full significance over the full range of arguments.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of sin(x).</returns>
        /// <remarks>
        /// <para>For an explanation of this method, see the remarks for the <see cref="MoreMath.Sin"/> method.</para>
        /// </remarks>
        public static double Cos (double x) {

            // Ensure x is non-negative; this is easy because cosine is an even function.
            if (x < 0.0) x = -x;

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
            } else if (x < Double.PositiveInfinity) {
                // If x is beyond the range of the built-in function entirely, use our range-reduction algorithm
                return (RangeReduction.Cos(x));
            } else {
                return (Double.NaN);
            }
        }

        /// <summary>
        /// Computes the sine of the given multiple of &#x3C0;.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of sin(<paramref name="x"/>&#x3C0;).</returns>
        /// <remarks>
        /// <para>Formulas involving sin(&#x3C0;x) appear in many contexts. By using this method
        /// instead of <c>Math.Sin(Math.PI * x)</c>, you will increase performance and avoid inaccuracies
        /// that arise from the finite precision of the stored constant <see cref="Math.PI"/>.</para>
        /// <para>Suppose, for example, x = 1.0E6. Since x is an integer, sin(&#x3C0;x) = 0.0.
        /// However, due to the finite accuracy of Math.PI, Math.PI * x is not a perfect multiple of &#x3C0;, and
        /// <c>Math.Sin(Math.PI * x)</c> = -2.2318717360358953E-10. But <c>MoreMath.SinPi(x)</c> = 0.0 exactly. Even for
        /// arguments that are not exact integers, the accurary of MoreMath.SinPi will be better.
        /// </para>
        /// </remarks>
        /// <seealso cref="CosPi(double)"/>
        /// <seealso cref="TanPi(double)"/>
        public static double SinPi (double x) {
            RangeReduction.ReduceByOnes(2.0 * x, out long y0, out double y1);
            return (RangeReduction.Sin(y0, y1));
        }

        /// <summary>
        /// Computes the cosine of the given multiple of &#x3C0;.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of cos(<paramref name="x"/>&#x3C0;).</returns>
        /// <remarks>
        /// <para>For an explanation of why and when to use this function,
        /// see <see cref="SinPi(double)"/>.</para>
        /// </remarks>
        /// <seealso cref="SinPi(double)"/>
        /// <seealso cref="TanPi(double)"/>
        public static double CosPi (double x) {
            RangeReduction.ReduceByOnes(2.0 * x, out long y0, out double y1);
            return (RangeReduction.Cos(y0, y1));
        }

        /// <summary>
        /// Computes the tangent of the given multiple of &#x3C0;.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of tan(&#x3C0; <paramref name="x"/>).</returns>
        /// <remarks>
        /// <para>For an explanation of why and when to use this function,
        /// see <see cref="SinPi(double)"/>.</para>
        /// </remarks>
        /// <seealso cref="SinPi(double)"/>
        /// <seealso cref="CosPi(double)"/>
        public static double TanPi (double x) {
            RangeReduction.ReduceByOnes(2.0 * x, out long y0, out double y1);
            if (y0 % 2L == 0L) {
                return (Math.Tan(Math.PI / 2.0 * y1));
            } else {
                return (-1.0 / Math.Tan(Math.PI / 2.0 * y1));
            }
            // Should be possible to do even better, by reducing wrt pi / 4, but this is good enough for now.
        }

        /// <summary>
        /// Computes the sinc function.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of sinc(x) = sin(x) / x.</returns>
        /// <remarks>
        /// <para>The ratio sin(x) / x, understood to take its limiting value of 1 at x = 0, appears in many applications.
        /// This method allows you to compute it straightforwardly without having to implement the limit.</para>
        /// <para>Note that in signal processing applications, the function sin(&#x3C0;x)/(&#x3C0;x) is usually
        /// used instead, and is also called the sinc function. Meta.Numerics implements it as <see cref="SincPi(double)"/>.</para>
        /// <para>This function is also sometimes called the cardinal since function or the sampling function.</para>
        /// </remarks>
        /// <seealso href="http://mathworld.wolfram.com/SincFunction.html"/>
        /// <seealso href="https://en.wikipedia.org/wiki/Sinc_function"/>
        public static double Sinc (double x) {
            // Perhaps counter-intuitively, naive evaluation of sin(x) / x looses no accuracy,
            // even very close to zero, so long as x \ne 0. So we only branch on that one value.
            if (x == 0.0) {
                return 1.0;
            } else {
                return Sin(x) / x;
            }
        }

        /// <summary>
        /// Computes the sinc of the given multiple of &#x3C0;.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of sin(&#x3C0;x)/(&#x3C0;x).</returns>
        /// <remarks>
        /// <para>This definition of the sinc function is commonly used in signal processing applications.
        /// For the more common definition without the factor &#x3C0;, see <see cref="Sinc(double)"/>.
        /// </para>
        /// </remarks>
        public static double SincPi (double x) {
            if (x == 0.0) {
                return 1.0;
            } else {
                return SinPi(x) / (Math.PI * x);
            }
        }


        /// <summary>
        /// Returns the value of n mod m.
        /// </summary>
        /// <param name="n">The argument.</param>
        /// <param name="m">The modulus.</param>
        /// <returns>The value of n mod m.</returns>
        /// <remarks>The modulus operator in .NET languages (n % m in C#, n Mod m in VB) returns
        /// a negative remainder for negative n.
        /// This is not consistent with the mathematical conventions of modular arithmetic,
        /// which require that 0 &#x2264; n mod m &lt; m. For example, in C# -7 % 4 returns -3
        /// (which is simply the negative of 7 % 4), while mathematically -7 mod 4 = 1, as can
        /// seen by extending the repeating 0, 1, 2, 3 pattern to negative integers. In cases
        /// where it is important to obtain results consistent mathematical conventions,
        /// you should call this method, which fixes the .NET result.</remarks>
        public static long Mod (long n, long m) {
            long r = n % m;
            if (r < 0) r += m;
            return r;
            // ((n % m) + m) % m would probably be faster (mod being faster than test-and-branch),
            // but it would not only be less clear, but also suffer overflow problems for m near integer limits
        }

        /// <summary>
        /// Computes the midpoint of two integers.
        /// </summary>
        /// <param name="x">One integer.</param>
        /// <param name="y">Another integer.</param>
        /// <returns>The midpoint of the two integers (rounded down).</returns>
        /// <remarks>
        /// <para>The simple expression (x + y)/2 can overflow. The implementation of this
        /// method returns the accurate midpoint for all values of x and y, including all
        /// combinations of their extreme values.</para>
        /// </remarks>
        public static int Midpoint(int x, int y) {
            // This expression is from Warren, "Hacker's Delight", Section 2-5
            return (x & y) + ((x ^ y) >> 1);
        }

    }


    internal static class RangeReduction {

        // Start with Dekker arithmetic. Suppose x and y are floating point numbers with a limited mantissa, say 32 bits; then
        // x * y will have up to 64 bits, the last 32 bits of which will be lost when put the result into another 32-bit register.
        // But suppose we decompose x = xHi + xLo and y = yHi + yLo, with each part using only 16 bits but still being stored
        // in a 32-bit register. Then we can form products like xHi * yHi and xHi * yLo and be sure the product will fit fully
        // into a 32-bit register, so we can form their sum, product, etc. without any rounding error. This is Dekker arithmetic.
        // It's basically the idea as representing a long using two ints, but applied to floating point types.

        // Start with the decomposition, which is called Veltkamp splitting:

        public static void Decompose (double x, out double hi, out double lo) {
            double p = K * x;
            if (Double.IsInfinity(p) && !Double.IsInfinity(x)) {
                Debug.Assert(Math.Abs(x) >= Double.MaxValue / K);
                const int scaleFactor = 1 << (bitsPerPart + 1);
                double xPrime = x / scaleFactor;
                Decompose(xPrime, out double hiPrime, out double loPrime);
                hi = hiPrime * scaleFactor;
                lo = loPrime * scaleFactor;
                return;
            }
            double q = x - p;
            hi = p + q;
            lo = x - hi;
        }

        public static double[] Decompose (double x) {
            Decompose(x, out double hi, out double lo);
            return (new double[] { hi, lo });
        }

        // The number of bits to keep in each part when doing Veltkamp splitting

        private const int bitsPerPart = 26;

        // The corresponding constant used in the Veltkamp splitting algorithm

        private const double K = (1 << bitsPerPart) + 1;


        // Next comes arithmetic. All we need is multiplication.

        // The binary representation of 2/\pi. We will use this to find z = (2\pi) x.

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

        /*
        public static void Add (double a, double b, out double s, out double r) {
            if (Math.Abs(a) < Math.Abs(b)) Global.Swap(ref a, ref b);
            s = a + b;
            double z = s - a;
            r = b - z;
        }





        // This is the Dekker multiplication algorithm.

        // Multiply x by y. x is a double split into high and low parts xHi and xLo. y is a high-accuracy float whose
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
        */

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

        public static double Sin (long z0, double z1) {
            switch (MoreMath.Mod(z0, 4L)) {
                case 0L:
                    return (Math.Sin(Math.PI / 2.0 * z1));
                case 1L:
                    return (Math.Cos(Math.PI / 2.0 * z1));
                case 2L:
                    return (-Math.Sin(Math.PI / 2.0 * z1));
                case 3L:
                    return (-Math.Cos(Math.PI / 2.0 * z1));
                default:
                    throw new InvalidOperationException();
            }

        }

        public static double Cos (long z0, double z1) {
            switch (MoreMath.Mod(z0, 4L)) {
                case 0L:
                    return (Math.Cos(Math.PI / 2.0 * z1));
                case 1L:
                    return (-Math.Sin(Math.PI / 2.0 * z1));
                case 2L:
                    return (-Math.Cos(Math.PI / 2.0 * z1));
                case 3L:
                    return (Math.Sin(Math.PI / 2.0 * z1));
                default:
                    throw new InvalidOperationException();
            }
        }

        public static double Sin (double x) {
            ReduceByPiHalves(x, out long z0, out double z1);
            return (Sin(z0, z1));
        }

        public static double Cos (double x) {
            ReduceByPiHalves(x, out long z0, out double z1);
            return (Cos(z0, z1));
        }

        public static void ReduceByPiHalves (double x, out long x0, out double x1) {
            double[] xParts = Decompose(x);
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

        /*
        public static double Cos (double x, double y) {
            long x0; double x1;
            ReduceByPiHalves(x, out x0, out x1);

            long y0; double y1;
            ReduceByOnes(y, out y0, out y1);

            long z0 = x0 + y0;
            double z1 = x1 + y1;

            return (Cos(z0, z1));

        }
        */

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

        /*
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
        */

    }


}