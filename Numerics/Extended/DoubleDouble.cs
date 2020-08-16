using System;
using System.Diagnostics;
using System.Globalization;
using System.Text;

namespace Meta.Numerics.Extended {

    /// <summary>
    /// Represents a floating point number with quadruple precision.
    /// </summary>
    /// <remarks>
    /// <para>The <see cref="DoubleDouble"/> structure uses two <see cref="Double"/> values to achieve
    /// twice the precision with which a floating point number can be stored and manipulated as compared to
    /// to the <see cref="Double"/> structure, approximately 31 decimal digits.</para>
    /// <para>Of all the extended precision floating point systems, double double is the
    /// fastest when implemented in software. A typical floating point operation using
    /// <see cref="DoubleDouble"/>s is just 3-4 times slower than using <see cref="Double"/>s.</para>
    /// <para>To instantiate a <see cref="DoubleDouble"/>, you can use <see cref="DoubleDouble.TryParse(string, out DoubleDouble)"/>,
    /// or <see cref="DoubleDouble.Parse(string)"/>, or the constructor <see cref="DoubleDouble.DoubleDouble(string)"/>
    /// to parse the text representation of the decimal value you want. If the value you want can be represented as a <see cref="Double"/>
    /// or <see cref="Int32"/> or other built in numeric type, you can cast that value to a <see cref="DoubleDouble"/>.</para>
    /// <para>When casting a <see cref="Double"/> to a <see cref="DoubleDouble"/>, there is a gotcha that you must be careful
    /// to avoid. Suppose you write <c>DoubleDouble x = 0.2</c> or <c>DoubleDouble x = 1.0 / 5.0</c>. You might think that this produces the 
    /// <see cref="DoubleDouble"/> representation of 1/5, but you would be wrong. The problem is that the compiler intreprets 0.2 or 1.0/5.0
    /// as <see cref="Double"/>s, and 1/5th is not exactly representable as a double, since it is not a rational number with a power-of-two denominator.
    /// Double's best attempt at 1/5th is 3602879701896397 X 2^<sup>-54</sup> = 0.20000000000000001110223024625157..., which
    /// is accurate to 16 decimal digits, but not to 32. Therefore when it is cast to a <see cref="DoubleDouble"/> it is much
    /// farther away from 1/5th than <see cref="DoubleDouble"/> can achieve. To obtain 1/5th to the accuracy of a <see cref="DoubleDouble"/>,
    /// you must write <c>DoubleDouble x = new DoubleDouble("0.2")</c> or <c>DoubleDouble x = (DoubleDouble) 1 / 5</c>. (The latter works
    /// because 1 and 5 are exactly representable and the division is performed as <see cref="DoubleDouble"/> division. All integers in range,
    /// indeed all rational numbers with in-range numerators and power-of-two denominators, are exactly representable. So, for example,
    /// <c>DoubleDouble x = 0.25</c> <i>does</i> work as expected, because 1/4 is exactly representable. But to avoid the gotcha
    /// it's best to simply train yourself to avoid assigning <see cref="DoubleDouble"/> variables from factional <see cref="Double"/> values.)</para>
    /// <para>Many of the mathematical functions which are implemented for <see cref="Double"/> arguments by static methods of the <see cref="Math"/> class
    /// are implemented for <see cref="DoubleDouble"/> arguments by static methods of the <see cref="DoubleDouble"/> type itself,
    /// for example <see cref="DoubleDouble.Sqrt(DoubleDouble)"/> and <see cref="DoubleDouble.Log(DoubleDouble)"/>.
    /// Some of the advanced functions which are implemented for <see cref="Double"/> arguments by static methods of the <see cref="Meta.Numerics.Functions.AdvancedMath"/> class
    /// are implemented for <see cref="DoubleDouble"/> arguments by static methods of the <see cref="AdvancedDoubleDoubleMath"/> class.
    /// </para>
    /// <para>You may wonder why <see cref="DoubleDouble"/> is not simply named "Quad". The reason is that "Quad" would propertly refer to an
    /// implementation of the <see href="https://en.wikipedia.org/wiki/Quadruple-precision_floating-point_format">IEEE 754 quadruple-precision binary floating
    /// point format</see>, which would have not only the extended presion of <see cref="DoubleDouble"/>, but also an extended range (up to 10<sup>4932</sup>).
    /// </para>
    /// </remarks>
    public struct DoubleDouble : IEquatable<DoubleDouble>, IComparable<DoubleDouble> {

        internal DoubleDouble (double hi, double lo) {
            // Unless number is zero or NaN, lo parts should be lower than hi part by 2^52
            Debug.Assert((hi == 0.0 && lo == 0.0) || Math.Abs(lo) < (1.0E-14) * Math.Abs(hi) || Double.IsNaN(hi));
            this.hi = hi;
            this.lo = lo;
        }

        /// <summary>
        /// Initializes a new double double number from the given string.
        /// </summary>
        /// <param name="s">The base-10 representation of the number.</param>
        public DoubleDouble (string s) {
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

        /// <summary>
        /// The double double negative infinite value.
        /// </summary>
        public static readonly DoubleDouble NegativeInfinity = new DoubleDouble(Double.NegativeInfinity, 0.0);

        /// <summary>
        /// The double double positive infinite value.
        /// </summary>
        public static readonly DoubleDouble PositiveInfinity = new DoubleDouble(Double.PositiveInfinity, 0.0);

        /// <summary>
        /// The double double not-a-number value.
        /// </summary>
        public static readonly DoubleDouble NaN = new DoubleDouble(Double.NaN, 0.0);

        private static readonly DoubleDouble ten = new DoubleDouble(10.0, 0.0);

        private static readonly DoubleDouble log2 = new DoubleDouble("0.693147180559945309417232121458176568");

        /// <summary>
        /// The double double value of pi.
        /// </summary>
        /// <seealso href="https://en.wikipedia.org/wiki/Pi"/>
        public static readonly DoubleDouble Pi = DoubleDouble.Parse("3.141592653589793238462643383279502884");

        /// <summary>
        /// The double double value of the base of natural logarithms.
        /// </summary>
        /// <seealso href="https://en.wikipedia.org/wiki/E_(mathematical_constant)"/>
        public static readonly DoubleDouble E = DoubleDouble.Parse("2.71828182845904523536028747135266249");

        // Arithmetic

        /// <summary>
        /// Negates a double double number. 
        /// </summary>
        /// <param name="x">The number to negate.</param>
        /// <returns>The additive inverse of <paramref name="x"/>.</returns>
        public static DoubleDouble operator - (DoubleDouble x) {
            return new DoubleDouble(-x.hi, -x.lo);
        }

        /// <summary>
        /// Computes the sum of two double double numbers.
        /// </summary>
        /// <param name="x">The first number.</param>
        /// <param name="y">The second number.</param>
        /// <returns>The value of <paramref name="x"/> + <paramref name="y"/>.</returns>
        public static DoubleDouble operator + (DoubleDouble x, DoubleDouble y) {

            // Add high components
            ExtendedMath.TwoSum(x.hi, y.hi, out double sHi, out double sLo);

            if (ExtendedMath.IsNotFinite(sHi)) return (DoubleDouble) sHi;

            // Add low components
            ExtendedMath.TwoSum(x.lo, y.lo, out double tHi, out double tLo);
            ExtendedMath.TwoSum(sHi, sLo + tHi, out double vHi, out double vLo);
            ExtendedMath.FastTwoSum(vHi, tLo + vLo, out double zHi, out double zLo);

            return new DoubleDouble(zHi, zLo);
        }

        /// <summary>
        /// Computes the difference of two double double numbers.
        /// </summary>
        /// <param name="x">The first number.</param>
        /// <param name="y">The second number.</param>
        /// <returns>The value of <paramref name="x"/> - <paramref name="y"/>.</returns>
        public static DoubleDouble operator - (DoubleDouble x, DoubleDouble y) {
            return (x + (-y));
        }

        /// <summary>
        /// Computes the product of two double double numbers.
        /// </summary>
        /// <param name="x">The first number.</param>
        /// <param name="y">The second number.</param>
        /// <returns>The product of <paramref name="x"/> and <paramref name="y"/>.</returns>
        public static DoubleDouble operator * (DoubleDouble x, DoubleDouble y) {

            ExtendedMath.TwoProduct(x.hi, y.hi, out double p0, out double p1);

            if (p0 == 0.0 || ExtendedMath.IsNotFinite(p0)) return (DoubleDouble) p0;

            ExtendedMath.TwoProduct(x.hi, y.lo, out double p2, out double p4);

            ExtendedMath.TwoProduct(x.lo, y.hi, out double p3, out double p5);

            double p6 = x.lo * y.lo;

            ExtendedMath.ThreeSum(p1, p2, p3, out double t1, out double t2);

            t2 += p4 + p5 + p6;

            ExtendedMath.ThreeSum(p0, t1, t2, out double pHi, out double pLo);

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
        /// <param name="x">The dividend.</param>
        /// <param name="y">The divisor.</param>
        /// <returns>The value of <paramref name="x"/> / <paramref name="y"/>.</returns>
        public static DoubleDouble operator / (DoubleDouble x, DoubleDouble y) {

            double q1 = x.hi / y.hi;

            // If leading order result is NaN or infinity or zero, we are done.
            // To continue would introduce NaNs even if result is infinite, so this early return is necessary.
            if (q1 == 0.0 || ExtendedMath.IsNotFinite(q1)) return (DoubleDouble) q1;

            DoubleDouble r = x - q1 * y;
            double q2 = r.hi / y.hi;
            r = r - q2 * y;
            double q3 = r.hi / y.hi;

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
        /// <param name="x">The value to be converted.</param>
        /// <returns>The <see cref="Double"/> value closest to the original double double value.</returns>
        /// <remarks><para>Note that this is a narrowing operator; the extra precision of the double double type
        /// is lost in the conversion.</para></remarks>
        public static explicit operator double (DoubleDouble x) {
            return (x.hi + x.lo);
        }

        /// <summary>
        /// Converts a double value to a double double value.
        /// </summary>
        /// <param name="x">The value to be converted.</param>
        /// <returns>The double double value equal to the original <see cref="Double" /> value.</returns>
        /// <remarks><para>This cast preserves the value of the <see cref="Double"/> (as can be verified
        /// by round-tripping), but keep in mind that that value may only be an approximation of less
        /// than the desired precision.</para></remarks>
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
                x = Sqr(x);
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

            if (x.hi == 0.0) return (0.0);

            double yHi = Math.Sqrt(x.hi);

            if (ExtendedMath.IsNotFinite(yHi)) return (DoubleDouble) yHi;

            ExtendedMath.TwoProduct(yHi, yHi, out double uHi, out double uLo);

            double yLo = (((x.hi - uHi) - uLo) + x.lo) / (2.0 * yHi);

            return new DoubleDouble(yHi, yLo);

        }

        internal static DoubleDouble Log1P (DoubleDouble x) {
            if (DoubleDouble.One + x == DoubleDouble.One) {
                return (x);
            } else {
                DoubleDouble mx = -x;
                DoubleDouble t = x;
                DoubleDouble f = t;
                for (int k = 2; k < Global.SeriesMax; k++) {
                    DoubleDouble f_old = f;
                    t *= mx;
                    f += t / k;
                    if (f == f_old) {
                        return (f);
                    }
                }
                throw new NonconvergenceException();
            }
        }

        /// <summary>
        /// Computes the natural logarithm of a double double value.
        /// </summary>
        /// <param name="x">The argument of the logarithm.</param>
        /// <returns>The value of ln(x).</returns>
        public static DoubleDouble Log (DoubleDouble x) {

            double logHi = Math.Log(x.hi);

            if (ExtendedMath.IsNotFinite(logHi)) return (DoubleDouble) logHi;

            int e = (int) Math.Round(logHi / Global.LogTwo);
            if (e < 0) {
                x *= DoubleDouble.Pow(2.0, -e);
            } else if (e > 0) {
                x /= DoubleDouble.Pow(2.0, e);
            }
            // At this point 1/\sqrt{2} <= r <= \sqrt{2},
            // i.e.  0.707 <= r <= 1.414
            x -= DoubleDouble.One;
            // Now -0.293 <= r - 1 <= 0.414.
            // We have lost some accuracy if r ~ 1, i.e. x was very close
            // to an exact power of 2.

            return e * log2 + Log1P(x);
        }

        private static DoubleDouble Exp_Series (DoubleDouble x) {
            DoubleDouble t = x;
            DoubleDouble f = DoubleDouble.One + t;
            for (int k = 2; k < Global.SeriesMax; k++) {
                DoubleDouble f_old = f;
                t *= x / k;
                f += t;
                if (f == f_old) {
                    return (f);
                }
            }
            throw new NonconvergenceException();
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

        /// <summary>
        /// Computes the sine of a double double value.
        /// </summary>
        /// <param name="x">The argument of the sine.</param>
        /// <returns>The value of sin(x).</returns>
        /// <remarks><para>Unlike <see cref="MoreMath.Sin(double)"/>, this function does not currenly do range reduction with a value
        /// of pi even more accurate than <see cref="DoubleDouble.Pi"/>. Therefore, beyond the first few periods, results can slowly
        /// loose accuracy, particularly near zeros, at the rate of about one digit per order of magnitude of the argument.</para></remarks>
        public static DoubleDouble Sin (DoubleDouble x) {
            return DoubleDoubleMath.Sin(x);
        }

        /// <summary>
        /// Computes the cosine of a double double value.
        /// </summary>
        /// <param name="x">The argument of the cosine.</param>
        /// <returns>The value of cos(x).</returns>
        /// <remarks><para>Unlike <see cref="MoreMath.Cos(double)"/>, this function does not currenly do range reduction with a value
        /// of pi even more accurate than <see cref="DoubleDouble.Pi"/>. Therefore, beyond the first few periods, results can slowly
        /// loose accuracy, particularly near zeros, at the rate of about one digit per order of magnitude of the argument.</para></remarks>
        public static DoubleDouble Cos (DoubleDouble x) {
            return DoubleDoubleMath.Cos(x);
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
            ParseResult result = TryParseInternal(s, out DoubleDouble x);
            if (result == ParseResult.Null) throw new ArgumentNullException(nameof(s));
            if (result != ParseResult.Success) throw new FormatException();
            return x;
        }

        /// <summary>
        /// Attempts to parse a string representation of a double double value.
        /// </summary>
        /// <param name="s">The string representation of the value.</param>
        /// <param name="x">The value, if the parse was successful.</param>
        /// <returns>True if the string was successfully parsed, false otherwise.</returns>
        /// <remarks>
        /// <para>Double double supports the same string representations as
        /// <see cref="Double"/>.</para>
        /// </remarks>
        public static bool TryParse (string s, out DoubleDouble x) {
            ParseResult result = TryParseInternal(s, out x);
            return (result == ParseResult.Success);
        }

        private static ParseResult TryParseInternal (string s, out DoubleDouble x) {

            x = DoubleDouble.Zero;

            if (s is null) return ParseResult.Null;

            s = s.Trim();

            if (s.Length == 0) return ParseResult.Empty;

            if (s == "NaN") {
                x = DoubleDouble.NaN;
                return ParseResult.Success;
            }

            // Parse sign.
            bool negative = false;
            if (s[0] == '-') {
                negative = true;
                s = s.Substring(1);
            } else if (s[0] == '+') {
                s = s.Substring(1);
            }

            // Parse exponent.
            int exponent = 0;
            int eIndex = s.IndexOfAny(new char[] { 'e', 'E' });
            if (eIndex >= 0) {
                bool exponentResult = Int32.TryParse(s.Substring(eIndex + 1), out exponent);
                if (!exponentResult) return ParseResult.Format;
                s = s.Substring(0, eIndex);
            }

            // Convert decimal value into a straight-up integer with an adjusted exponent
            int dIndex = s.IndexOf('.');
            if (dIndex >= 0) {
                string s0 = s.Substring(0, dIndex);
                string s1 = s.Substring(dIndex + 1);
                exponent -= s1.Length;
                s = s0 + s1;
            }

            // We have done a lot of trimming. If nothing remains, that was not a valid format.
            if (s.Length == 0) return ParseResult.Format;

            // Drop leading zeros
            s = s.TrimStart('0');

            // Absorb trailing zeros into the exponent
            int l0 = s.Length;
            s = s.TrimEnd('0');
            int l1 = s.Length;
            exponent += (l0 - l1);

            // Compute the integer
            foreach (char c in s) {
                double d = Char.GetNumericValue(c);
                if (d < 0.0) return ParseResult.Format;
                x *= ten;
                x += d;
            }

            // Apply the exponent
            if (exponent != 0) {
                x *= DoubleDouble.Pow(ten, exponent);
            }

            if (negative) x = -x;

            return ParseResult.Success;

        }

        /// <summary>
        /// Produces a string representation of the double double value.
        /// </summary>
        /// <returns>A string representation of the value.</returns>
        public override string ToString () {
            return ToString(CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// Produces a text representation of the double double value using the given format provider.
        /// </summary>
        /// <param name="format">The format provider.</param>
        /// <returns>A text representation of the value.</returns>
        private string ToString (IFormatProvider format) {

            // This algorithm doesn't handle infinities, NaNs, and zeros, but writing them is trivial. 
            if (ExtendedMath.IsNotFinite(hi) || hi == 0.0) {
                return hi.ToString(format);
            }

            // Multiply by a power of 10 to put leading digit in the ones place
            Debug.Assert(hi != 0.0);
            int e = (int) Math.Floor(Math.Log10(Math.Abs(hi)));
            DoubleDouble r = Abs(this);
            if (e > 0) {
                r /= Pow(ten, e);
            } else if (e < 0) {
                r *= Pow(ten, -e);
            }

            // If scaling has screwed us up because of overflow, exit. Fix this.
            if (Double.IsNaN(r.hi) || Double.IsNaN(r.lo)) {
                return "X";
            }

            // Spit out up to 32 digits
            StringBuilder s = new StringBuilder();
            int z = 0;
            for (int i = 0; i < 32; i++) {
                if (r.hi == 0.0) break;
                double t = Math.Floor(r.hi); /* If nothing left, we can stop early */
                if (t < 0) t = 0.0; /* if we subtract to zero, sometimes the result can be a very tiny negative number that floors to -1 */
                int d = (int) t;
                Debug.Assert((0 <= d) && (d < 10));
                s.Append(d);
                if (d == 0) { z++; } else { z = 0; } /* Keep track of number of trailing zeros for later truncation */
                r -= t;
                r *= ten;
            }

            // Remove trailing zeros (which are all after the decimal point)
            if (z > 0) s.Remove(s.Length - z, z);

            // Express in scientific notation for large exponents,
            // or by filling in zeros for small exponents.
            if (e >= 6) {
                if (s.Length > 1) s.Insert(1, '.');
                s.Append($"E{e}");
            } else if (e >= 0) {
                if (s.Length > e + 1) {
                    s.Insert(e + 1, '.');
                } else {
                    s.Append("000000", 0, e - s.Length + 1);
                }
            } else if (e >= -4) {
                s.Insert(0, "0.00000".Substring(0, 1 - e)); /* StringBuilder.Append has a substring overload, but StringBuilder.Insert doesn't */
            } else {
                if (s.Length > 1) s.Insert(1, '.');
                s.Append($"E{e}");
            }

            // Add a negative sign if necessary
            if (hi < 0.0) s.Insert(0, '-');

            return s.ToString();

        }

        // Equality

        /// <summary>
        /// Determines whether two double double values are equal.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns><see langword="true"/> if <paramref name="x"/> and <paramref name="y"/> are equal, otherwise <see langword="false"/>.</returns>
        public static bool Equals (DoubleDouble x, DoubleDouble y) {
            // Is this really enough? I see no bugs in our tests, but isn't it possible for (hi+1, -lo) to represent same value as (hi,lo)?
            return ((x.hi == y.hi) && (x.lo == y.lo));
        }

        /// <summary>
        /// Determines whether two double double values are equal.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns><see langword="true"/> if <paramref name="x"/> and <paramref name="y"/> are equal, otherwise <see langword="false"/>.</returns>
        public static bool operator == (DoubleDouble x, DoubleDouble y) {
            return (Equals(x, y));
        }

        /// <summary>
        /// Determines whether two double double values are unequal.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns><see langword="false"/> if <paramref name="x"/> and <paramref name="y"/> are equal, otherwise <see langword="true"/>.</returns>
        public static bool operator != (DoubleDouble x, DoubleDouble y) {
            return (!Equals(x, y));
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
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns><see langword="true"/> if <paramref name="x"/> is less than <paramref name="y"/>, otherwise <see langword="false"/>.</returns>
        public static bool operator < (DoubleDouble x, DoubleDouble y) {
            return ((x.hi < y.hi) || ((x.hi == y.hi) && (x.lo < y.lo)));
        }

        /// <summary>
        /// Determines whether the first value is greater than the second value.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns><see langword="true"/> if <paramref name="x"/> is greater than <paramref name="y"/>, otherwise <see langword="false"/>.</returns>
        public static bool operator > (DoubleDouble x, DoubleDouble y) {
            return ((x.hi > y.hi) || ((x.hi == y.hi) && (x.lo > y.lo)));
        }

        /// <summary>
        /// Determines whether the first value is less than or equal to the second value.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns><see langword="true"/> if <paramref name="x"/> is less than or equal to <paramref name="y"/>, otherwise <see langword="false"/>.</returns>
        public static bool operator <= (DoubleDouble x, DoubleDouble y) {
            return ((x.hi < y.hi) || ((x.hi == y.hi) && (x.lo <= y.lo)));
        }

        /// <summary>
        /// Determines whether the first value is greater than or equal to the second value.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns><see langword="true"/> if <paramref name="x"/> is greater than or equal to <paramref name="y"/>, otherwise <see langword="false"/>.</returns>
        public static bool operator >= (DoubleDouble x, DoubleDouble y) {
            return ((x.hi > y.hi) || ((x.hi == y.hi) && (x.lo >= y.lo)));
        }

        /// <summary>
        /// Compares two values.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns>-1 if <paramref name="x"/> is less than <paramref name="y"/>, +1 if <paramref name="x"/> is greater than
        /// <paramref name="y"/>, 0 if <paramref name="x"/> and <paramref name="y"/> are equal.</returns>
        public static int Compare (DoubleDouble x, DoubleDouble y) {
            if (x < y) {
                return -1;
            } else if (x > y) {
                return +1;
            } else {
                return 0;
            }
        }

        /// <summary>
        /// Compares the current value to another.
        /// </summary>
        /// <param name="other">The other value to compare.</param>
        /// <returns>-1 if this value is less than the other value, +1 if it is greater than the other value, 0 if they are equal.</returns>
        public int CompareTo (DoubleDouble other) {
            return Compare(this, other);
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

        private const double eps = 1.0 / (1L << 54);

        /// <summary>
        /// Determines whether a double double value is not-a-number.
        /// </summary>
        /// <param name="x">The value.</param>
        /// <returns>True is <paramref name="x"/> is not-a-number, otherwise false.</returns>
        /// <remarks><para>By the floating point standard respect by <see cref="Double"/> and floating point types in
        /// essentially all languages and frameworks, equality testing NaN always returns false. Therefore it is
        /// necessary to have a specific method to test for not-a-number. <see cref="Double.IsNaN(double)"/> is that
        /// method for <see cref="Double"/>, and this is that method for <see cref="DoubleDouble"/>.</para></remarks>
        public static bool IsNaN (DoubleDouble x) {
            return Double.IsNaN(x.hi);
        }

    }

}
