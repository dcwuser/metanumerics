using System;
using System.Diagnostics;
using System.Numerics;

namespace Meta.Numerics.Extended
{

    /// <summary>
    /// Represents a signed integer with a 128 big register width.
    /// </summary>
    /// <remarks>
    /// <para>To instantiate a 128-bit signed integer, you can use the <see cref="Int128.Parse(string)"/> or <see cref="Int128.TryParse(string, out Int128)"/>
    /// methods to parse a decimal string representation. If the value you want is representable by a built-in integer type (e.g. <see cref="Int64"/>), you can
    /// simply assign a <see cref="Int128"/>-typed variable from a built-in integer type.</para>
    /// <para>If you know that the numbers you need fit in a 128 bit register, arithmetic with <see cref="Int128"/> is typically faster than with <see cref="BigInteger"/>.
    /// Addition and subtraction are about four times faster, multiplication is about twice as fast, and division is about equally fast.</para>
    /// <para>Explicit casts between <see cref="Int128"/> and the built-in integer types behave like unchecked explicit casts between the built-in integer types:
    /// the low-order bits of the binary representation are preserved and the high-order bits are, if necessary, discarded. This preserves values that are representable
    /// by the target type, but values that are not representable by the target type may be transformed in ways that, while correct in terms of the bit-level rules,
    /// are unexpected in terms of numerical values. This is different than the behavior of <see cref="BigInteger"/>, which performs casts as if checked, throwing
    /// an exception if the value is not representable in the target type. If you want checked behavior, your code must explicitly check whether the value is in
    /// the range of the target before casting.</para>
    /// <para>The signed 128-bit integer type does not support bit-wise logical and shift operations, but the unsigned 128-bit integer type <see cref="UInt128"/>
    /// does. If you want to perform bit-wise operations on 128-bit registers, use <see cref="UInt128"/>.</para>
    /// </remarks>
    [CLSCompliant(false)]
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
        public static readonly Int128 MaxValue = new Int128(ulong.MaxValue >> 1, ulong.MaxValue);

        /// <summary>
        /// The minimum value of the signed 128-bit integer.
        /// </summary>
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
        /// Converts an aribitrary-sized big integer into a 128-bit integer.
        /// </summary>
        /// <param name="b">Te big integer to be converted.</param>
        /// <returns>The 128-bit integer with the same lower 128 bits.</returns>
        public static explicit operator Int128 (BigInteger b) {
            UInt128 u = (UInt128) b;
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

        // Absolute Value

        /// <summary>
        /// Gets the absolute value of a 128 bit integer.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The absolute value of the argument.</returns>
        /// <remarks><para>Because <see cref="Int128.MinValue"/> is one unit smaller in absolute value than <see cref="Int128.MaxValue"/>, this
        /// method throws an <see cref="OverflowException"/> when passed <see cref="Int128.MinValue"/>. Built in types such as <see cref="Int32"/>
        /// have analogous behavior. All other values are supported.</para></remarks>
        /// <exception cref="OverflowException"><paramref name="x"/> was <see cref="Int128.MinValue"/>.</exception>
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
