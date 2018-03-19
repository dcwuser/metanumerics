using System;
using System.Diagnostics;
using System.Text;

namespace Meta.Numerics.Extended {

    /// <summary>
    /// Represents a floating point number with extended precision.
    /// </summary>
    /// <remarks>
    /// <para>The double double format uses two <see cref="Double"/> values to effectively
    /// double the precision with which a number can be stored and manipulated as compared to
    /// to <see cref="Double"/>, to approximately 31 decimal digits of precision.</para>
    /// <para>Of all the extended precision floating point systems, double double is the
    /// fastest when implemented in software. A typical floating point operation on
    /// double doubles is 3-4 times slower than on <see cref="Double"/>s.</para>
    /// </remarks>
    public struct DoubleDouble : IEquatable<DoubleDouble>, IComparable<DoubleDouble> {
        internal DoubleDouble (double hi, double lo) {
            Debug.Assert(Math.Abs(lo) <= (1.0E-10) * Math.Abs(hi));
            this.hi = hi;
            this.lo = lo;
        }

        /// <summary>
        /// Initializes a new double double number from the given string.
        /// </summary>
        /// <param name="s">The base-10 representation of the number.</param>
        public DoubleDouble (string s) {
            if (s == null) throw new ArgumentNullException(nameof(s));
            DoubleDouble r = DoubleDouble.Parse(s);
            this.hi = r.hi;
            this.lo = r.lo;
        }

        private readonly double hi;
        private readonly double lo;

        // Constant values

        /// <summary>
        /// The double double number zero.
        /// </summary>
        public static readonly DoubleDouble Zero = new DoubleDouble(0.0, 0.0);

        /// <summary>
        /// The double double number one.
        /// </summary>
        public static readonly DoubleDouble One = new DoubleDouble(1.0, 0.0);

        private static readonly DoubleDouble ten = new DoubleDouble(10.0, 0.0);

        private static readonly DoubleDouble log2 = new DoubleDouble("0.693147180559945309417232121458176568");

        /// <summary>
        /// The double double value of pi.
        /// </summary>
        public static readonly DoubleDouble Pi = DoubleDouble.Parse("3.141592653589793238462643383279502884");

        /// <summary>
        /// The double double vale of the base of natural logarithms.
        /// </summary>
        public static readonly DoubleDouble E = DoubleDouble.Parse("2.71828182845904523536028747135266249");

        private static readonly DoubleDouble Euler = DoubleDouble.Parse("0.57721566490153286060651209008240243");



        // Arithmetic

        /// <summary>
        /// Negates a double double number. 
        /// </summary>
        /// <param name="a">The number to negate.</param>
        /// <returns>The additive inverse of <paramref name="a"/>.</returns>
        public static DoubleDouble operator - (DoubleDouble a) {
            return new DoubleDouble(-a.hi, -a.lo);
        }

        /// <summary>
        /// Computes the sum of two double double numbers.
        /// </summary>
        /// <param name="a">The first number.</param>
        /// <param name="b">The second number.</param>
        /// <returns>The value of <paramref name="a"/> + <paramref name="b"/>.</returns>
        public static DoubleDouble operator + (DoubleDouble a, DoubleDouble b) {
            // Add high components
            double sHi, sLo;
            ExtendedMath.TwoSum(a.hi, b.hi, out sHi, out sLo);

            // Add low components
            double tHi, tLo;
            ExtendedMath.TwoSum(a.lo, b.lo, out tHi, out tLo);

            double vHi, vLo;
            ExtendedMath.TwoSum(sHi, sLo + tHi, out vHi, out vLo);

            double zHi, zLo;
            ExtendedMath.FastTwoSum(vHi, tLo + vLo, out zHi, out zLo);

            return new DoubleDouble(zHi, zLo);
        }

        /// <summary>
        /// Computes the difference of two double double numbers.
        /// </summary>
        /// <param name="a">The first number.</param>
        /// <param name="b">The second number.</param>
        /// <returns>The value of <paramref name="a"/> - <paramref name="b"/>.</returns>
        public static DoubleDouble operator - (DoubleDouble a, DoubleDouble b) {
            return (a + (-b));
        }

        /// <summary>
        /// Computes the product of two double double numbers.
        /// </summary>
        /// <param name="a">The first number.</param>
        /// <param name="b">The second number.</param>
        /// <returns>The product of <paramref name="a"/> and <paramref name="b"/>.</returns>
        public static DoubleDouble operator * (DoubleDouble a, DoubleDouble b) {
            double p0, p1;
            ExtendedMath.TwoProduct(a.hi, b.hi, out p0, out p1);

            double p2, p4;
            ExtendedMath.TwoProduct(a.hi, b.lo, out p2, out p4);

            double p3, p5;
            ExtendedMath.TwoProduct(a.lo, b.hi, out p3, out p5);

            double p6 = a.lo * b.lo;

            double t1, t2;
            ExtendedMath.ThreeSum(p1, p2, p3, out t1, out t2);

            t2 += p4 + p5 + p6;

            double pHi, pLo;
            ExtendedMath.ThreeSum(p0, t1, t2, out pHi, out pLo);

            return new DoubleDouble(pHi, pLo);
            /*
            double zHi, zLo;
            ExtendedMath.TwoProduct(a.hi, b.hi, out zHi, out zLo);
            zLo = (zLo + a.hi * b.lo) + a.lo * b.hi;
            //zLo = (a.hi * b.lo + a.lo * b.hi) + zLo;

            return new DoubleDouble(zHi, zLo);
            */
        }

        /// <summary>
        /// Computes the quotient of two double double numbers.
        /// </summary>
        /// <param name="a">The dividend.</param>
        /// <param name="b">The divisor.</param>
        /// <returns>The value of <paramref name="a"/> / <paramref name="b"/>.</returns>
        public static DoubleDouble operator / (DoubleDouble a, DoubleDouble b) {
            double q1 = a.hi / b.hi;
            DoubleDouble r = a - q1 * b;
            double q2 = r.hi / b.hi;
            r = r - q2 * b;
            double q3 = r.hi / b.hi;

            double qHi, qLo;

            double t1, t2, t3;
            ExtendedMath.TwoSum(q1, q2, out t1, out t2);
            ExtendedMath.TwoSum(q3, t1, out qHi, out t3);
            ExtendedMath.TwoSum(t2, t3, out qLo, out t1);

            return new DoubleDouble(qHi, qLo);

            /*
            // Compute an initial approximation by dividing hi parts.
            double qHi = a.hi / b.hi;

            // Compute the product of qHi and bHi to full precision.
            // Ideally this would be exactly aHi but it won't be.
            double uHi, uLo;
            ExtendedMath.TwoProduct(qHi, b.hi, out uHi, out uLo);

            // Now compute the correction and but it in qLo.
            double qLo = ((((a.hi - uHi) - uLo) + a.lo) - qHi * b.lo) / b.hi;

            return new DoubleDouble(qHi, qLo);
            */
        }

        // Casts

        /// <summary>
        /// Converts a double double value to a double value.
        /// </summary>
        /// <param name="a">The value to be converted.</param>
        /// <returns>The <see cref="Double"/> value closest to the original double double value.</returns>
        /// <remarks><para>Note that this is a narrowing operator; the extra precision of the double double type
        /// is lost in the conversion.</para></remarks>
        public static explicit operator double (DoubleDouble a) {
            return (a.hi + a.lo);
        }

        /// <summary>
        /// Converts a double value to a double double value.
        /// </summary>
        /// <param name="x">The value to be converted.</param>
        /// <returns>The double double value equal to the original <see cref="Double" /> value.</returns>
        public static implicit operator DoubleDouble (double x) {
            return new DoubleDouble(x, 0.0);
        }

        // Basic functions

        /// <summary>
        /// Computes the absolute value of a double double number.
        /// </summary>
        /// <param name="x">The number.</param>
        /// <returns>The absolute value of <paramref name="x"/>.</returns>
        public static DoubleDouble Abs (DoubleDouble x) {
            // Don't just take abs of both components, since they may have different signs.
            if (x.hi < 0.0) {
                return new DoubleDouble(-x.hi, -x.lo);
            } else {
                return (x);
            }
        }

        /// <summary>
        /// Raises a double double number to an integer power.
        /// </summary>
        /// <param name="x">The number.</param>
        /// <param name="n">The exponent.</param>
        /// <returns>The value of x<sup>n</sup>.</returns>
        public static DoubleDouble Pow (DoubleDouble x, int n) {
            if (n == 0) return DoubleDouble.One;

            if (n < 0) {
                x = DoubleDouble.One / x;
                n = -n;
            }

            if (n == 1) return x;

            Debug.Assert(n > 1);
            /*
            // Exponentiation by squaring
            DoubleDouble y = DoubleDouble.One;
            while (n > 0)
            {
                if (n % 2 != 0)
                {
                    n -= 1;
                    y *= x;
                }
                n /= 2;
                x *= x;
            }
            return (y);

            DoubleDouble y = DoubleDouble.One;
            while (true)
            {
                if (n % 2 != 0)
                {
                    y *= x;
                }
                n /= 2;
                if (n == 0)
                {
                    return (y);
                }
                x *= x;
            }
            */

            // Exponentiation by squaring
            DoubleDouble y = DoubleDouble.One;
            while (n > 1) {
                if (n % 2 != 0) y = x * y;
                x = Sqr(x); /* Replace with optimized square function */
                n /= 2; /* Integer division effectively deals with odd n case by dropping the remainder */
            }
            return x * y;
        }

        private static DoubleDouble Sqr (DoubleDouble x) {
            return (x * x);
        }

        /// <summary>
        /// Computes the square root of a double double value.
        /// </summary>
        /// <param name="x">The value of which the square root will be computed.</param>
        /// <returns>The value of the square root of x.</returns>
        public static DoubleDouble Sqrt (DoubleDouble x) {
            if (Double.IsNaN(x.hi)) return (x);

            if (x.hi < 0.0) return (Double.NaN);

            if (x.hi == 0.0) return (0.0);

            double yHi = Math.Sqrt(x.hi);

            double uHi, uLo;
            ExtendedMath.TwoProduct(yHi, yHi, out uHi, out uLo);

            double yLo = (((x.hi - uHi) - uLo) + x.lo) / (2.0 * yHi);

            return new DoubleDouble(yHi, yLo);
        }

        private static DoubleDouble Log1P (DoubleDouble x) {
            if (DoubleDouble.One + x == DoubleDouble.One) {
                return (x);
            } else {
                DoubleDouble mx = -x;
                DoubleDouble t = x;
                DoubleDouble f = t;
                for (int k = 2; k < 100; k++) {
                    DoubleDouble f_old = f;
                    t *= mx;
                    f += t / k;
                    if (f == f_old) {
                        return (f);
                    }
                }
                throw new InvalidOperationException();
            }
        }

        /// <summary>
        /// Computes the natural logarithm of a double double value.
        /// </summary>
        /// <param name="x">The argument of the logarithm.</param>
        /// <returns>The value of ln(x).</returns>
        public static DoubleDouble Log (DoubleDouble x) {
            if (Double.IsNaN(x.hi)) return (x);
            if (x.hi < 0.0) return (Double.NaN);
            if (x.hi == 0.0) return (Double.NegativeInfinity);

            DoubleDouble r = x;
            int e = (int) Math.Round(Math.Log(x.hi) / Math.Log(2.0));
            if (e < 0) {
                r *= DoubleDouble.Pow(2.0, -e);
            } else if (e > 0) {
                r /= DoubleDouble.Pow(2.0, e);
            }
            // At this point 1/\sqrt{2} <= r <= \sqrt{2},
            // i.e.  0.707 <= r <= 1.414
            r -= DoubleDouble.One;
            // Now -0.293 <= r - 1 <= 0.414.
            // We have lost some accuracy if r ~ 1, i.e. x was very close
            // to an exact power of 2.

            return e * log2 + Log1P(r);
        }

        private static DoubleDouble Exp_Series (DoubleDouble x) {
            DoubleDouble t = x;
            DoubleDouble f = DoubleDouble.One + t;
            for (int k = 2; k < 100; k++) {
                DoubleDouble f_old = f;
                t *= x / k;
                f += t;
                if (f == f_old) {
                    return (f);
                }
            }
            throw new InvalidOperationException();
        }

        /// <summary>
        /// Computes the exponential of a double double value.
        /// </summary>
        /// <param name="x">The argument of the exponential.</param>
        /// <returns>The value of e<sup>x</sup>.</returns>
        public static DoubleDouble Exp (DoubleDouble x) {
            if (x.hi == 0.0) return (DoubleDouble.One);

            DoubleDouble m = x / log2;
            int n = (int) Math.Round(m.hi);
            x -= n * log2;

            return DoubleDouble.Pow(2.0, n) * Exp_Series(x);

        }

        // Printing and Parsing

        /// <summary>
        /// Parses a string representation of a double double value.
        /// </summary>
        /// <param name="s">The string representation of the value.</param>
        /// <returns>The corresponding double double value.</returns>
        /// <remarks>
        /// <para>Double double supports the same string representations as
        /// <see cref="Double"/>.</para>
        /// </remarks>
        public static DoubleDouble Parse (string s) {
            if (s == null) throw new ArgumentNullException(nameof(s));

            s = s.Trim();

            // Parse negative sign.
            bool isNegative = false;
            if (s[0] == '-') {
                isNegative = true;
                s = s.Substring(1);
            }

            // Parse exponent.
            int exponent = 0;
            int eIndex = s.IndexOfAny(new char[] { 'e', 'E' });
            if (eIndex >= 0) {
                exponent = Int32.Parse(s.Substring(eIndex + 1));
                s = s.Substring(0, eIndex);
            }

            // Remove trailing and leading 0s. Otherwise we do unnecessary multiplications
            // (and later divisions) by 10, which also introduces floating point wobble.
            s = s.Trim('0');

            // Parse the mantissa.
            int digitCount = 0;
            int decimalPoint = -1;
            DoubleDouble result = 0.0;
            foreach (char c in s) {
                if (c == '.') {
                    if (decimalPoint > 0) throw new FormatException();
                    decimalPoint = digitCount;
                } else {
                    double value = Char.GetNumericValue(c);
                    result *= ten;
                    result += value;
                    digitCount++;
                }
            }

            if (decimalPoint >= 0) {
                exponent += decimalPoint - digitCount;
            }

            if (exponent != 0) {
                result *= DoubleDouble.Pow(ten, exponent);
            }

            if (isNegative) {
                result = -result;
            }

            return (result);
        }

        private string Write () {
            if (Double.IsNaN(hi) || Double.IsInfinity(hi)) {
                return (hi.ToString());
            } else if (hi == 0.0) {
                return ((0.0).ToString());
            } else {
                // Multiply by a power of 10 to put leading digit in the ones place
                int e = (int) Math.Floor(Math.Log10(Math.Abs(hi)));
                DoubleDouble r = Abs(this);
                if (e > 0) {
                    r /= Pow(ten, e);
                } else if (e < 0) {
                    r *= Pow(ten, -e);
                }

                StringBuilder s = new StringBuilder();
                if (hi < 0.0) s.Append('-');

                // This passes unit tests, but it needs to be cleaned up.
                //   1. Don't print unnecessary trailing 0s.
                //   2. For -2 <= e <= 2, don't use E notation.
                //   3. Be more careful with rounding.
                for (int i = 0; i < 32; i++) {
                    double t = Math.Floor(r.hi);
                    if (t < 0) t = 0.0; /* if we subtract to zero, sometimes the result can be a very tiny negative number that floors to -1 */
                    int d = (int) t;
                    Debug.Assert((0 <= d) && (d < 10));
                    s.Append(d);
                    if (i == 0) s.Append(".");
                    r -= t;
                    r *= ten;
                }

                if (e != 0) s.Append(String.Format("E{0}", e));
                return (s.ToString());
            }
        }

        /// <summary>
        /// Writes a string representation of the value.
        /// </summary>
        /// <returns>A string representation of the value.</returns>
        public override string ToString () {
            return (Write());
        }

        // Equality

        /// <summary>
        /// Determines whether two double double values are equal.
        /// </summary>
        /// <param name="a">The first value.</param>
        /// <param name="b">The second value.</param>
        /// <returns><see langword="true"/> if <paramref name="a"/> and <paramref name="b"/> are equal, otherwise <see langword="false"/>.</returns>
        private static bool Equals (DoubleDouble a, DoubleDouble b) {
            return ((a.hi == b.hi) && (a.lo == b.lo));
        }

        /// <summary>
        /// Determines whether two double double values are equal.
        /// </summary>
        /// <param name="a">The first value.</param>
        /// <param name="b">The second value.</param>
        /// <returns><see langword="true"/> if <paramref name="a"/> and <paramref name="b"/> are equal, otherwise <see langword="false"/>.</returns>
        public static bool operator == (DoubleDouble a, DoubleDouble b) {
            return (Equals(a, b));
        }

        /// <summary>
        /// Determines whether two double double values are unequal.
        /// </summary>
        /// <param name="a">The first value.</param>
        /// <param name="b">The second value.</param>
        /// <returns><see langword="false"/> if <paramref name="a"/> and <paramref name="b"/> are equal, otherwise <see langword="true"/>.</returns>
        public static bool operator != (DoubleDouble a, DoubleDouble b) {
            return (!Equals(a, b));
        }

        /// <summary>
        /// Determines whether the current value is equal to another value.
        /// </summary>
        /// <param name="other">The other value to compare.</param>
        /// <returns><see langword="true"/> if <paramref name="other"/> is equal to the current value, otherwise <see langword="false"/>.</returns>
        public bool Equals (DoubleDouble other) {
            return (Equals(this, other));
        }

        /// <summary>
        /// Determines whether the current value is equal to another object.
        /// </summary>
        /// <param name="obj">The other value to compare.</param>
        /// <returns><see langword="true"/> if <paramref name="obj"/> is equal to the current value, otherwise <see langword="false"/>.</returns>
        public override bool Equals (object obj) {
            return ((obj is DoubleDouble) && Equals(this, (DoubleDouble) obj));
        }

        /// <summary>
        /// Gets a hash code for the current value.
        /// </summary>
        /// <returns>A hash of the current value.</returns>
        public override int GetHashCode () {
            unchecked {
                return (hi.GetHashCode() + 31 * lo.GetHashCode());
            }
        }

        // Comparison

        /// <summary>
        /// Determines whether the first value is less than the second value.
        /// </summary>
        /// <param name="a">The first value.</param>
        /// <param name="b">The second value.</param>
        /// <returns><see langword="true"/> if <paramref name="a"/> is less than <paramref name="b"/>, otherwise <see langword="false"/>.</returns>
        public static bool operator < (DoubleDouble a, DoubleDouble b) {
            return ((a.hi < b.hi) || ((a.hi == b.hi) && (a.lo < b.lo)));
        }

        /// <summary>
        /// Determines whether the first value is greater than the second value.
        /// </summary>
        /// <param name="a">The first value.</param>
        /// <param name="b">The second value.</param>
        /// <returns><see langword="true"/> if <paramref name="a"/> is greater than <paramref name="b"/>, otherwise <see langword="false"/>.</returns>
        public static bool operator > (DoubleDouble a, DoubleDouble b) {
            return ((a.hi > b.hi) || ((a.hi == b.hi) && (a.lo > b.lo)));
        }

        /// <summary>
        /// Determines whether the first value is less than or equal to the second value.
        /// </summary>
        /// <param name="a">The first value.</param>
        /// <param name="b">The second value.</param>
        /// <returns><see langword="true"/> if <paramref name="a"/> is less than or equal to <paramref name="b"/>, otherwise <see langword="false"/>.</returns>
        public static bool operator <= (DoubleDouble a, DoubleDouble b) {
            return ((a.hi < b.hi) || ((a.hi == b.hi) && (a.lo <= b.lo)));
        }

        /// <summary>
        /// Determines whether the first value is greater than or equal to the second value.
        /// </summary>
        /// <param name="a">The first value.</param>
        /// <param name="b">The second value.</param>
        /// <returns><see langword="true"/> if <paramref name="a"/> is greater than or equal to <paramref name="b"/>, otherwise <see langword="false"/>.</returns>
        public static bool operator >= (DoubleDouble a, DoubleDouble b) {
            return ((a.hi > b.hi) || ((a.hi == b.hi) && (a.lo >= b.lo)));
        }

        /// <summary>
        /// Compares the current value to another.
        /// </summary>
        /// <param name="other">The other value to compare.</param>
        /// <returns>-1 if this value is less than the other value, +1 if it is greater than the other value, 0 if they are equal.</returns>
        public int CompareTo (DoubleDouble other) {
            if (this < other) {
                return (-1);
            } else if (this > other) {
                return (+1);
            } else {
                return (0);
            }
        }

        // random

        /// <summary>
        /// Gets a random double double value.
        /// </summary>
        /// <param name="rng">A random number generator.</param>
        /// <returns>A random value in [0, 1).</returns>
        public static DoubleDouble GetRandomValue (Random rng) {
            if (rng == null) throw new ArgumentNullException(nameof(rng));
            double hi = rng.NextDouble();
            double lo = rng.NextDouble() * hi * eps;
            return (new DoubleDouble(hi, lo));
        }

        private static readonly double eps = 1.0 / (1L << 54);

    }

}
