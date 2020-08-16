using System;
using System.Diagnostics;
using System.Numerics;

namespace Meta.Numerics.Extended
{

    /// <summary>
    /// Represents a signed integer with a 128 bit register width.
    /// </summary>
    /// <remarks>
    /// <para>The built-in integer types short, int, and long are fixed-width registers with the characteristics shown below. <see cref="Int128"/> is
    /// a 128-bit signed integer register with analogous behavior and greatly extended range.</para>
    /// <table>
    ///     <tr>
    ///         <th>Type</th>
    ///         <th>C# Name</th>
    ///         <th>Width</th>
    ///         <th>MaxValue</th>
    ///     </tr>
    ///     <tr>
    ///         <td><see cref="Int16"/></td>
    ///         <td>short</td>
    ///         <td>16b (2B)</td>
    ///         <td>~32 X 10<sup>3</sup></td>
    ///     </tr>
    ///     <tr>
    ///         <td><see cref="Int32"/></td>
    ///         <td>int</td>
    ///         <td>32b (4B)</td>
    ///         <td>~2.1 X 10<sup>9</sup></td>
    ///     </tr>
    ///     <tr>
    ///         <td><see cref="Int64"/></td>
    ///         <td>long</td>
    ///         <td>64b (8B)</td>
    ///         <td><see cref="Int64.MaxValue">~9.2 X 10<sup>18</sup></see></td>
    ///     </tr>
    ///     <tr>
    ///         <td><see cref="Int128"/></td>
    ///         <td><see cref="Int128"/></td>
    ///         <td>128b (16B)</td>
    ///         <td><see cref="Int128.MaxValue">~1.7 X 10<sup>38</sup></see></td>
    ///     </tr>
    /// </table>
    /// <para>To instantiate a 128-bit signed integer, you can use the constructor <see cref="Int128.Int128(string)"/>, or the parsing methods <see cref="Int128.Parse(string)"/>
    /// or <see cref="Int128.TryParse(string, out Int128)"/> to parse a decimal representation provided as a string. The parsing then occurs at run-time.
    /// If the value you want is representable by a built-in integer type (e.g. <see cref="Int64"/>), you can
    /// simply assign a <see cref="Int128"/> variable from a built-in integer type. This is particularly useful in source code, since the compiler will parse fixed
    /// values of these types at compile-time.</para>
    /// <para>If you know that the numbers you need to work with all fit in a 128 bit register, arithmetic with <see cref="Int128"/> is typically significantly faster than with <see cref="BigInteger"/>.
    /// Addition and subtraction are about six times faster, multiplication is about twice as fast, and division is about 50% faster. (Meta.Numerics itself uses <see cref="Int128"/>
    /// internally in the implementation of some discrete distributions that require very large integers.)</para>
    /// <para>Operations that overflow and underflow <see cref="Int128"/> behave like they do for the native unsigned integer types, with results equal to the lower 128 bits of the
    /// full result, or wrapping around from <see cref="Int128.MaxValue"/> to <see cref="Int128.MinValue"/>.</para>
    /// <para>Explicit casts between <see cref="Int128"/> and the built-in integer types behave like unchecked explicit casts between the built-in integer types:
    /// the low-order bits of the binary representation are preserved and the high-order bits are, if necessary, discarded. This preserves values that are representable
    /// by the target type, but values that are not representable by the target type may be transformed in ways that, while correct in terms of the bit-level rules,
    /// are unexpected in terms of numerical values. While this is like the behavior of the built-in integer types, it is different than the behavior of <see cref="BigInteger"/>,
    /// which performs casts as if checked, throwing an exception if the value is not representable in the target type. If you want checked behavior, your code must explicitly
    /// check whether the value is in the range of the target before casting.</para>
    /// <para><see cref="Int128"/> does not support bit-wise logical and shift operations, but the unsigned 128-bit integer type <see cref="UInt128"/>
    /// does. If you want to perform bit-wise operations on 128-bit registers, use <see cref="UInt128"/>.</para>
    /// </remarks>
    public struct Int128 : IEquatable<Int128>, IComparable<Int128> {

        /// <summary>
        /// Initializes a new 128-bit integer with the given decimal representation.
        /// </summary>
        /// <param name="s">The decimal representation of th 128-bit integer.</param>
        public Int128 (string s) {
            Int128 x = Int128.Parse(s);
            this.u = x.u;
        }

        private readonly UInt128 u;

        private Int128 (ulong s1, ulong s0) {
            this.u = new UInt128(s1, s0);
        }

        private Int128 (UInt128 u) {
            this.u = u;
        }

        // Static values

        /// <summary>
        /// The zero value of the signed 128-bit integer.
        /// </summary>
        public static readonly Int128 Zero = new Int128(0L, 0UL);

        /// <summary>
        /// The unit value of the signed 128-bit integer.
        /// </summary>
        public static readonly Int128 One = new Int128(0L, 1UL);

        /// <summary>
        /// The maximum value of the signed 128-bit integer.
        /// </summary>
        /// <remarks>This has the value 170,141,183,460,469,231,731,687,303,715,884,105,728, or about 1.7 X 10<sup>38</sup>.</remarks>
        public static readonly Int128 MaxValue = new Int128(ulong.MaxValue >> 1, ulong.MaxValue);

        /// <summary>
        /// The minimum value of the signed 128-bit integer.
        /// </summary>
        /// <remarks><para>This has the value -170,141,183,460,469,231,731,687,303,715,884,105,729, or about -1.7 X 10<sup>38</sup>.</para>
        /// <para>As is true for the native signed integer registers and any signed integer register that uses two's complement representation
        /// of negative numbers, this is one larger in absolute value than <see cref="Int128.MaxValue"/>, so the corresponding
        /// positive value is not within the representable range of the type.</para></remarks>
        public static readonly Int128 MinValue = new Int128(1UL << 63, 0UL);

        // Equality

        /// <summary>
        /// Tests the equality of two 128-bit integers.
        /// </summary>
        /// <param name="x">The first integer.</param>
        /// <param name="y">The second integer.</param>
        /// <returns><see langword="true"/> if <paramref name="x"/> equals <paramref name="y"/>, <see langword="false"/> otherwise.</returns>
        public static bool Equals (Int128 x, Int128 y) {
            return (x.u == y.u);
        }

        /// <summary>
        /// Tests whether the current instance equals another 128-bit integer.
        /// </summary>
        /// <param name="other">The other integer.</param>
        /// <returns><see langword="true"/> if the current instance equals <paramref name="other"/>, <see langword="false"/> otherwise.</returns>
        public bool Equals (Int128 other) {
            return Equals(this, other);
        }

        /// <summary>
        /// Tests the equality of two 128-bit integers.
        /// </summary>
        /// <param name="x">The first integer.</param>
        /// <param name="y">The second integer.</param>
        /// <returns><see langword="true"/> if <paramref name="x"/> equals <paramref name="y"/>, <see langword="false"/> otherwise.</returns>
        public static bool operator == (Int128 x, Int128 y) {
            return Equals(x, y);
        }

        /// <summary>
        /// Tests the inequality of two 128-bit integers.
        /// </summary>
        /// <param name="x">The first integer.</param>
        /// <param name="y">The second integer.</param>
        /// <returns><see langword="false"/> if <paramref name="x"/> equals <paramref name="y"/>, <see langword="true"/> otherwise.</returns>
        public static bool operator != (Int128 x, Int128 y) {
            return !Equals(x, y);
        }

        /// <summary>
        /// Tests whether the current instance equals another object.
        /// </summary>
        /// <param name="obj">The other object.</param>
        /// <returns><see langword="true"/> if the current instance equals <paramref name="obj"/>, <see langword="false"/> otherwise.</returns>
        public override bool Equals (object obj) {
            if (obj is Int128) {
                return Equals(this, (Int128) obj);
            } else {
                return false;
            }
        }

        /// <summary>
        /// Return a hash code for the 128-bit integer.
        /// </summary>
        /// <returns>A hash code for the value.</returns>
        public override int GetHashCode () {
            return u.GetHashCode() + 31;
        }

        // Comparison

        private static readonly UInt128 mask = UInt128.One << 127;

        /// <summary>
        /// Indicates whether the first value is less than the second value.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns><see langword="true" /> if <paramref name="x"/> is less than <paramref name="y"/>,
        /// otherwise <see langword="false"/>.</returns>
        public static bool operator < (Int128 x, Int128 y) {
            return (x.u ^ mask) < (y.u ^ mask);
        }

        /// <summary>
        /// Indicates whether the first value is greater than the second value.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns><see langword="true" /> if <paramref name="x"/> is greater than <paramref name="y"/>,
        /// otherwise <see langword="false"/>.</returns>
        public static bool operator > (Int128 x, Int128 y) {
            return (x.u ^ mask) > (y.u ^ mask);
        }

        /// <summary>
        /// Indicates whether the first value is greater than or equal to the second value.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns><see langword="true" /> if <paramref name="x"/> is greater than or equal to <paramref name="y"/>,
        /// otherwise <see langword="false"/>.</returns>
        public static bool operator >= (Int128 x, Int128 y) {
            return (x.u ^ mask) >= (y.u ^ mask);
        }

        /// <summary>
        /// Indicates whether the first value is less than or equal to the second value.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns><see langword="true" /> if <paramref name="x"/> is less than or equal to <paramref name="y"/>,
        /// otherwise <see langword="false"/>.</returns>
        public static bool operator <= (Int128 x, Int128 y) {
            return (x.u ^ mask) <= (y.u ^ mask);
        }

        /// <summary>
        /// Compares two 128-bit integers.
        /// </summary>
        /// <param name="x">The first integer.</param>
        /// <param name="y">The second integer.</param>
        /// <returns>1 if <paramref name="x"/> is greater, -1 if <paramref name="y"/> is greater,
        /// and 0 if <paramref name="x"/> and <paramref name="y"/> are equal.</returns>
        public static int Compare (Int128 x, Int128 y) {
            return UInt128.Compare(x.u ^ mask, y.u ^ mask);
        }

        /// <summary>
        /// Compares the current instance to another 128-bit integer.
        /// </summary>
        /// <param name="x">The other integer.</param>
        /// <returns>1 if the current instance is greater, -1 if <paramref name="x"/> is greater,
        /// and 0 if the current instance is equal to <paramref name="x"/>.</returns>
        public int CompareTo (Int128 x) {
            return Compare(this, x);
        }

        // Implicit conversions

        /// <summary>
        /// Converts a 64-bit integer into a 128-bit integer.
        /// </summary>
        /// <param name="s">The 64-bit integer to convert.</param>
        /// <returns>THe corresponding 128-bit integer.</returns>
        public static implicit operator Int128 (long s) {
            if (s < 0L) {
                return -(new Int128(0UL, (ulong) (-s)));
            } else {
                return new Int128(0UL, (ulong) s);
            }

        }

        /// <summary>
        /// Converts a 128-bit integer into a floating point value.
        /// </summary>
        /// <param name="s">The 128-bit integer to convert.</param>
        /// <returns>The correspoding floating point value.</returns>
        public static implicit operator double (Int128 s) {
            if (s.u.IsNegative) {
                return -((double) s.u.Negate());
            } else {
                return (double) s.u;
            }
        }

        /// <summary>
        /// Converts a 128-bit integer into an arbitrary-size big integer.
        /// </summary>
        /// <param name="s">The 128-bit integer to convert.</param>
        /// <returns>The corresponding big integer.</returns>
        public static implicit operator BigInteger (Int128 s) {
            byte[] bytes = s.u.GetBytesForBigInteger(true);
            return new BigInteger(bytes);
        }

        // Explicit conversions

        /// <summary>
        /// Converts a 128-bit integer into a 64-bit integer.
        /// </summary>
        /// <param name="s">The 128-bit integer to convert.</param>
        /// <returns>The 64-bit integer with the same lower 64 bits.</returns>
        public static explicit operator long (Int128 s) {
            return unchecked((long) (ulong) s.u);
        }

        /// <summary>
        /// Converts an arbitrary-size big integer into a 128-bit integer.
        /// </summary>
        /// <param name="b">The big integer to convert.</param>
        /// <returns>The 128-bit integer with the same lower 128 bits.</returns>
        public static explicit operator Int128 (BigInteger b) {
            UInt128 u = (UInt128) b;
            return new Int128(u);
        }

        /// <summary>
        /// Converts a floating-point value into a 128-bit integer.
        /// </summary>
        /// <param name="x">The floating-point number to convert.</param>
        /// <returns>The lower 128 bits of the integer part of the floating-point number.</returns>
        /// <exception cref="InvalidCastException"><paramref name="x"/> is NaN, or infinite.</exception>
        public static explicit operator Int128 (double x) {
            DoubleInfo s = new DoubleInfo(Math.Truncate(x));
            if (!s.IsFinite) throw new InvalidCastException();
            int e = s.Exponent;
            if (e < 0) return Int128.Zero;
            UInt128 u = (ulong) s.Mantissa;
            u = u << e;
            if (s.IsNegative) u = u.Negate();
            return new Int128(u);
        }

        // Serialization and deserialization

        /// <summary>
        /// Produces a string representation of the 128-bit integer.
        /// </summary>
        /// <returns>A string containing the base-10 representation of the 128-bit integer value.</returns>
        public override string ToString () {
            if (u.IsNegative) {
                UInt128 t = u.Negate();
                return "-" + t.ToString();
            } else {
                return u.ToString();
            }
        }

        /// <summary>
        /// Produces a 128-bit integer from its string representation.
        /// </summary>
        /// <param name="text">A string containing the base-10 representation of a 128-bit integer value.</param>
        /// <returns>The 128-bit integer represented by <paramref name="text"/>.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="text"/> was <see langword="null"/>.</exception>
        /// <exception cref="FormatException"><paramref name="text"/> did not contain a valid representation of an integer.</exception>
        /// <exception cref="OverflowException"><paramref name="text"/> represented an integer outside the range of <see cref="Int128"/>.</exception>
        public static Int128 Parse (string text) {
            ParseResult result = UInt128.TryParse(text, true, out UInt128 u);
            if (result == ParseResult.Null) throw new ArgumentNullException(nameof(text));
            if (result == ParseResult.Overflow) throw new OverflowException();
            if (result != ParseResult.Success) throw new FormatException();
            return new Int128(u);
        }

        /// <summary>
        /// Attempts to produce a 128-bit integer from a string representation.
        /// </summary>
        /// <param name="text">A string containing the base-10 representation of a 128-bit integer value.</param>
        /// <param name="value">If the parse was successful, the 128-bit integer represented by <paramref name="text"/>.</param>
        /// <returns><see langword="true"/> if the parse was successful, otherwise <see langword="false"/>.</returns>
        public static bool TryParse (string text, out Int128 value) {
            ParseResult result = UInt128.TryParse(text, true, out UInt128 u);
            if (result == ParseResult.Success) {
                value = new Int128(u);
                return true;
            } else {
                value = Int128.Zero;
                return false;
            }
        }

        // Negation

        /// <summary>
        /// Negates a 128 bit integer.
        /// </summary>
        /// <param name="x">The integer.</param>
        /// <returns>The addative inverse of <paramref name="x"/>.</returns>
        public static Int128 operator - (Int128 x) {
            return new Int128(x.u.Negate());
        }

        // Arithmetic

        /// <summary>
        /// Computes the sum of two 128-bit integers.
        /// </summary>
        /// <param name="x">The first integer.</param>
        /// <param name="y">The second integer.</param>
        /// <returns>Tne sum <paramref name="x"/> + <paramref name="y"/>.</returns>
        public static Int128 operator + (Int128 x, Int128 y) {
            return new Int128(x.u + y.u);
        }

        /// <summary>
        /// Computes the difference of two 128-bit integers.
        /// </summary>
        /// <param name="x">The first integer.</param>
        /// <param name="y">The second integer.</param>
        /// <returns>Tne difference <paramref name="x"/> - <paramref name="y"/>.</returns>
        public static Int128 operator - (Int128 x, Int128 y) {
            return new Int128(x.u - y.u);
        }

        /// <summary>
        /// Increments a 128-bit integer.
        /// </summary>
        /// <param name="x">The integer.</param>
        /// <returns>One more than <paramref name="x"/>.</returns>
        public static Int128 operator ++ (Int128 x) {
            UInt128 v = x.u;
            v++;
            return new Int128(v);
        }

        /// <summary>
        /// Decrements a 128-bit unsigned integer.
        /// </summary>
        /// <param name="x">The integer.</param>
        /// <returns>One less than <paramref name="x"/>.</returns>
        public static Int128 operator -- (Int128 x) {
            UInt128 v = x.u;
            v--;
            return new Int128(v);
        }

        // Multiplication

        /// <summary>
        /// Computes the product of two 128-bit integers.
        /// </summary>
        /// <param name="x">The first integer.</param>
        /// <param name="y">The second integer.</param>
        /// <returns>The product <paramref name="x"/> X <paramref name="y"/>.</returns>
        public static Int128 operator * (Int128 x, Int128 y) {
            return new Int128(x.u * y.u);
        }

        /// <summary>
        /// Computer the quotient and remainder of two 128-bit integers.
        /// </summary>
        /// <param name="x">The dividend.</param>
        /// <param name="y">The divisor.</param>
        /// <param name="r">The remainder <paramref name="x"/> % <paramref name="y"/>.</param>
        /// <returns>The quotient <paramref name="x"/> / <paramref name="y"/>.</returns>
        public static Int128 DivRem (Int128 x, Int128 y, out Int128 r) {

            bool xNegative = false;
            UInt128 xu = x.u;
            if (xu.IsNegative) {
                xu = xu.Negate();
                xNegative = true;
            }

            bool yNegative = false;
            UInt128 yu = y.u;
            if (yu.IsNegative) {
                yu = yu.Negate();
                yNegative = true;
            }

            UInt128 qu = UInt128.DivRem(xu, yu, out UInt128 ru);
            if (xNegative ^ yNegative) {
                qu = qu.Negate();
            }
            if (xNegative) {
                ru = ru.Negate();
            }

            r = new Int128(ru);
            return new Int128(qu);

        }

        /// <summary>
        /// Computes the quotient of two 128 bit integers.
        /// </summary>
        /// <param name="x">The dividend.</param>
        /// <param name="y">The divisor.</param>
        /// <returns>The quotient <paramref name="x"/> / <paramref name="y"/>.</returns>
        public static Int128 operator / (Int128 x, Int128 y) {

            bool negate = false;

            UInt128 xu = x.u;
            if (xu.IsNegative) {
                xu = xu.Negate();
                negate = !negate;
            }

            UInt128 yu = y.u;
            if (yu.IsNegative) {
                yu = yu.Negate();
                negate = !negate;
            }

            UInt128 qu = xu / yu;
            if (negate) qu = qu.Negate();

            return new Int128(qu);

        }

        /// <summary>
        /// Computes the remainder of two 128 bit integers.
        /// </summary>
        /// <param name="x">The dividend.</param>
        /// <param name="y">The divisor.</param>
        /// <returns>The remainder of <paramref name="x"/> divided by <paramref name="y"/>.</returns>
        public static Int128 operator % (Int128 x, Int128 y) {
            Int128.DivRem(x, y, out Int128 r);
            return r;
        }

        // Absolute Value

        /// <summary>
        /// Gets the absolute value of a 128 bit integer.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The absolute value of the argument.</returns>
        /// <remarks><para>Because <see cref="Int128.MinValue"/> is one unit smaller in absolute value than <see cref="Int128.MaxValue"/>, this
        /// method throws an <see cref="OverflowException"/> when passed <see cref="Int128.MinValue"/>. Built in types such as <see cref="Int32"/>
        /// have analogous behavior. All other values are supported.</para></remarks>
        /// <exception cref="OverflowException"><paramref name="x"/> was <see cref="Int128.MinValue"/>, which has no corresponding positive
        /// value in the range of the type.</exception>
        public static Int128 Abs (Int128 x) {
            UInt128 xu = x.u;
            if (xu.IsNegative) {
                xu = xu.Negate();
                if (xu.IsNegative) throw new OverflowException();
            }
            Debug.Assert(!xu.IsNegative);
            return new Int128(xu);
        }

    }


}
