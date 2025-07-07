using System;
using System.Diagnostics;
using System.Numerics;
using System.Text;

namespace Meta.Numerics.Extended
{



    /// <summary>
    /// Represents a 128-bit unsigned integer.
    /// </summary>
    /// <remarks>
    /// <para>The built-in unsigned integer types ushort, uint, and ulong are fixed-width registers with the characteristics shown below. <see cref="UInt128"/> is
    /// a 128-bit unsigned integer register with analogous behavior and greatly extended range.</para>
    /// <table>
    ///     <tr>
    ///         <th>Type</th>
    ///         <th>C# Name</th>
    ///         <th>Width</th>
    ///         <th>MaxValue</th>
    ///     </tr>
    ///     <tr>
    ///         <td><see cref="UInt16"/></td>
    ///         <td>ushort</td>
    ///         <td>16b (2B)</td>
    ///         <td>~65 X 10<sup>3</sup></td>
    ///     </tr>
    ///     <tr>
    ///         <td><see cref="UInt32"/></td>
    ///         <td>uint</td>
    ///         <td>32b (4B)</td>
    ///         <td>~4.2 X 10<sup>9</sup></td>
    ///     </tr>
    ///     <tr>
    ///         <td><see cref="UInt64"/></td>
    ///         <td>ulong</td>
    ///         <td>64b (8B)</td>
    ///         <td>~18 X 10<sup>18</sup></td>
    ///     </tr>
    ///     <tr>
    ///         <td><see cref="UInt128"/></td>
    ///         <td><see cref="UInt128"/></td>
    ///         <td>128b (16B)</td>
    ///         <td>~3.4 X 10<sup>38</sup></td>
    ///     </tr>
    /// </table>
    /// <para>To instantiate a 128-bit unsigned integer, you can use the constructor <see cref="UInt128.UInt128(string)"/>, or the parsing methods <see cref="UInt128.Parse(string)"/>
    /// or <see cref="UInt128.TryParse(string, out UInt128)"/> to parse a decimal representation provided as a string. The parsing then occurs at run-time.
    /// If the value you want is representable by a built-in unsigned integer type (e.g. <see cref="UInt64"/>), you can
    /// simply assign a <see cref="UInt128"/> variable from a built-in integer type. This is particularly useful in source code, since the compiler will parse fixed
    /// values of these types at compile-time.</para>
    /// <para>Operations that overflow and underflow <see cref="UInt128"/> behave like they do for the native unsigned integer types, with results modulo 2<sup>128</sup>, or wrapping around
    /// from <see cref="UInt128.MaxValue"/> to <see cref="UInt128.Zero"/>.</para>
    /// <para>Explicit casts between <see cref="UInt128"/> and the built-in integer types behave like unchecked explicit casts between the built-in integer types:
    /// the low-order bits of the binary representation are preserved and the high-order bits are, if necessary, discarded. This preserves values that are representable
    /// by the target type, but values that are not representable by the target type may be transformed in ways that, while correct in terms of the bit-level rules,
    /// are unexpected in terms of numerical values. While this is like the behavior of the built-in integer types, it is different than the behavior of <see cref="BigInteger"/>,
    /// which performs casts as if checked, throwing an exception if the value is not representable in the target type. If you want checked behavior, your code must explicitly
    /// check whether the value is in the range of the target before casting.</para>
    /// <para>If you know that the numbers you need to work with all fit in a 128 bit unsigned register, arithmetic with <see cref="UInt128"/> is typically significantly faster than with <see cref="BigInteger"/>.
    /// Addition and subtraction are about six times faster, multiplication is about twice as fast, and division is about 50% faster.</para>
    /// <para><see cref="UInt128"/> also supports bitwise logical and shift operations.</para>
    /// </remarks>
    [CLSCompliant(false)]
    public readonly struct UInt128 : IEquatable<UInt128>, IComparable<UInt128> {

        /// <summary>
        /// Initializes a new 128-bit unsigned integer with the given base-10 representation.
        /// </summary>
        /// <param name="s">The base-10 text representation of the unsigned 128-bit integer.</param>
        public UInt128 (string s) {
            UInt128 u = UInt128.Parse(s);
            this.s0 = u.s0;
            this.s1 = u.s1;
        }

        internal UInt128 (ulong s1, ulong s0) {
            this.s0 = s0;
            this.s1 = s1;
        }

        // s = s_1 2^64 + s_0

        private readonly ulong s0;
        private readonly ulong s1;

        // I also tested storing as four uints. This makes multiplication a little faster (since multiplication needs the split
        // values and if we store this way we already have them) but addition and subtraction more slower, so sticking with
        // two ulongs.

        // Special cases

        /// <summary>
        /// The zero value of the unsigned 128-bit integer.
        /// </summary>
        public static readonly UInt128 Zero = new UInt128(0UL, 0UL);

        /// <summary>
        /// The unit value of the unsigned 128-bit integer.
        /// </summary>
        public static readonly UInt128 One = new UInt128(0UL, 1UL);

        /// <summary>
        /// The least representable unsigned 128-bit integer.
        /// </summary>
        /// <remarks><para>Since there are no negative values of an unsigned integer, this is just the nunber zero.</para></remarks>
        public static readonly UInt128 MinValue = new UInt128(0UL, 0UL);

        /// <summary>
        /// The greatest representable unsigned 128-bit integer.
        /// </summary>
        /// <remarks>
        /// <para>This has the value 340,282,366,920,938,463,463,374,607,431,768,211,455,
        /// or about 3.4 X 10<sup>38</sup>.</para>
        /// </remarks>
        public static readonly UInt128 MaxValue = new UInt128(UInt64.MaxValue, UInt64.MaxValue);

        // Implicit (widening) casts

        /// <summary>
        /// Converts an unsigned 64-bit integer to an unsigned 128-bit integer.
        /// </summary>
        /// <param name="u">The unsigned 64-bit integer.</param>
        /// <returns>The equivilent 128-bit unsigned integer.</returns>
        public static implicit operator UInt128 (ulong u) {
            return (new UInt128(0UL, u));
        }

        /// <summary>
        /// Converts an unsigned 128-bit integer into an arbitrary-size big integer.
        /// </summary>
        /// <param name="u">The unisgned 128-bit integer.</param>
        /// <returns>The equivilent <see cref="BigInteger"/>.</returns>
        public static implicit operator BigInteger (UInt128 u) {
            byte[] bytes = u.GetBytesForBigInteger(false);
            return new BigInteger(bytes);
        }

        internal byte[] GetBytesForBigInteger (bool isSigned) {
            // If highest order bit is 1, BigInteger intreprets the number as negative
            // with a two's complement representation. To avoid this for large, positive
            // unsigned integers, add an extra zero byte.
            byte[] bytes;
            if (isSigned) {
                bytes = new byte[16];
            } else {
                bytes = new byte[17];
            }
            // It is possible to use smaller byte arrays for smaller numbers, but
            // that requires more tests-and-branches so leave code simple for now.
            Array.Copy(BitConverter.GetBytes(s0), 0, bytes, 0, 8);
            Array.Copy(BitConverter.GetBytes(s1), 0, bytes, 8, 8);
            return (bytes);
        }

        /// <summary>
        /// Converts an unsigned 128-bit integer into a floating point value.
        /// </summary>
        /// <param name="u">The unsigned 128-bit integer.</param>
        /// <returns>The nearest <see cref="Double"/> floating point value.</returns>
        /// <remarks><para>This cast will never fail, but it will loose precision
        /// for values larger than about 10<sup>16</sup>, because <see cref="System.Double"/>
        /// maintains only about 52 bits of precision.</para></remarks>
        public static implicit operator double (UInt128 u) {
            const double multiplier = (double) UInt64.MaxValue + 1.0;
            return ((double) u.s0 + (double) u.s1 * multiplier);
        }

        // Explicit (narrowing) casts

        /// <summary>
        /// Converts an unsigned 128-bit integer into an unsigned 64-bit integer.
        /// </summary>
        /// <param name="u">The unsigned 128-bit integer.</param>
        /// <returns>The unsigned 64-bit integer with the same lower 64 bits.</returns>
        /// <remarks><para>This is x narrowing cast, which discards high-order bits if
        /// the 128-bit integer is greater than <see cref="UInt64.MaxValue"/>.</para></remarks>
        public static explicit operator ulong (UInt128 u) {
            return u.s0;
        }

        /// <summary>
        /// Converts an arbitrary-size big integer into an unsigned 128-bit integer.
        /// </summary>
        /// <param name="b">A <see cref="BigInteger"/> value.</param>
        /// <returns>The unsigned 128-bit integer with the same lower 128 bits.</returns>
        public static explicit operator UInt128 (BigInteger b) {
            byte[] bytes = b.ToByteArray();
            // BigInteger produces the minimum length byte array to represent the number.
            // Negative numbers are represented in two's complement form, so highest bit is 1.
            // To distinguish these from positive numbers with highest bit 1, it adds an extra zero byte.
            byte[] buffer = new byte[16];
            if (bytes.Length < buffer.Length) {
                Array.Copy(bytes, buffer, bytes.Length);
                // If this was a negative number, add additional 1s to complete two's complement form.
                if (bytes[bytes.Length - 1] >> 7 != 0) {
                    for (int i = bytes.Length; i < buffer.Length; i++) buffer[i] = 0xff;
                }
            } else {
                Array.Copy(bytes, buffer, buffer.Length);
            }
            ulong s0 = BitConverter.ToUInt64(buffer, 0);
            ulong s1 = BitConverter.ToUInt64(buffer, 8);
            return new UInt128(s1, s0);
        }

        /// <summary>
        /// Converts a floating-point value into an unsigned 128-bit integer.
        /// </summary>
        /// <param name="x">A floating-point number.</param>
        /// <returns>The lower 128 bits of the integer part of the floating-point number.</returns>
        /// <exception cref="InvalidCastException"><paramref name="x"/> is negative, or NaN, or infinite.</exception>
        public static explicit operator UInt128 (double x) {
            DoubleInfo s = new DoubleInfo(Math.Truncate(x));
            if (s.IsNegative || !s.IsFinite) throw new InvalidCastException();
            int e = s.Exponent;
            if (e < 0) return UInt128.Zero;
            UInt128 u = (ulong) s.Mantissa;
            u = u << e;
            return u;
        }

        // Equality

        /// <summary>
        /// Tests whether two unsigned 128-bit integers are equal.
        /// </summary>
        /// <param name="x">The first integer.</param>
        /// <param name="y">The second integer.</param>
        /// <returns><see langword="true"/> if <paramref name="x"/> equals <paramref name="y"/>, otherwise <see langword="false"/>.</returns>
        public static bool Equals (UInt128 x, UInt128 y) {
            return (x.s0 == y.s0) && (x.s1 == y.s1);
        }

        /// <summary>
        /// Tests whether two unsigned 128-bit integers are equal.
        /// </summary>
        /// <param name="x">The first integer.</param>
        /// <param name="y">The second integer.</param>
        /// <returns><see langword="true"/> if <paramref name="x"/> equals <paramref name="y"/>, otherwise <see langword="false"/>.</returns>
        public static bool operator == (UInt128 x, UInt128 y) {
            return Equals(x, y);
        }

        /// <summary>
        /// Tests whether two unsigned 128-bit integers are unequal.
        /// </summary>
        /// <param name="x">The first integer.</param>
        /// <param name="y">The second integer.</param>
        /// <returns><see langword="false"/> if <paramref name="x"/> equals <paramref name="y"/>, otherwise <see langword="true"/>.</returns>
        public static bool operator != (UInt128 x, UInt128 y) {
            return !Equals(x, y);
        }

        /// <summary>
        /// Tests whether the current instance is equal to another unsigned 128-bit integer.
        /// </summary>
        /// <param name="x">The other integer.</param>
        /// <returns><see langword="true"/> if the current instance equals <paramref name="x"/>, otherwise <see langword="false"/>.</returns>
        public bool Equals (UInt128 x) {
            return Equals(this, x);
        }

        /// <summary>
        /// Tests whether the current instance is equal to another object.
        /// </summary>
        /// <param name="obj">The other object.</param>
        /// <returns><see langword="true"/> if the other object is an equal <see cref="UInt128"/>, otherwise <see langword="false"/>.</returns>
        public override bool Equals (object obj) {
            if (obj is UInt128) {
                return (Equals(this, (UInt128) obj));
            } else {
                return (false);
            }
        }

        /// <summary>
        /// Gets a hash code for the current instance.
        /// </summary>
        /// <returns>A hash code for the current instance.</returns>
        public override int GetHashCode () {
            unchecked {
                return (s0.GetHashCode() + 13 * s1.GetHashCode());
            }
        }

        // Methods used by Int128, which uses this type as storage.

        internal bool IsNegative => ((s1 >> 63) != 0UL);

        internal UInt128 Negate () {
            Int128Calculator.TwosComplement(s1, s0, out ulong y1, out ulong y0);
            return new UInt128(y1, y0);
        }

        // Comparison

        /// <summary>
        /// Compares two unsigned 128-bit integers.
        /// </summary>
        /// <param name="x">The first integer.</param>
        /// <param name="y">The second integer.</param>
        /// <returns>1 if <paramref name="x"/> is greater, -1 if <paramref name="y"/> is greater,
        /// and 0 if <paramref name="x"/> and <paramref name="y"/> are equal.</returns>
        public static int Compare (UInt128 x, UInt128 y) {
            int c1 = x.s1.CompareTo(y.s1);
            if (c1 == 0) {
                Debug.Assert(x.s1 == y.s1);
                int c0 = x.s0.CompareTo(y.s0);
                return (c0);
            } else {
                return (c1);
            }
        }

        /// <summary>
        /// Compares the current instance to another unsigned 128-bit integer.
        /// </summary>
        /// <param name="u">The other integer.</param>
        /// <returns>1 if the current instance is greater, -1 if <paramref name="u"/> is greater,
        /// and 0 if the current instance is equal to <paramref name="u"/>.</returns>
        public int CompareTo (UInt128 u) {
            return (Compare(this, u));
        }

        /// <summary>
        /// Indicates whether the first value is greater than the second value.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns><see langword="true" /> if <paramref name="x"/> is greater than <paramref name="y"/>,
        /// otherwise <see langword="false"/>.</returns>
        public static bool operator > (UInt128 x, UInt128 y) {
            return (x.s1 > y.s1) || (x.s1 == y.s1 && x.s0 > y.s0);
        }

        /// <summary>
        /// Indicates whether the first value is less than the second value.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns><see langword="true" /> if <paramref name="x"/> is less than <paramref name="y"/>,
        /// otherwise <see langword="false"/>.</returns>
        public static bool operator < (UInt128 x, UInt128 y) {
            return (x.s1 < y.s1) || (x.s1 == y.s1 && x.s0 < y.s0);
        }

        /// <summary>
        /// Indicates whether the first value is greater than or equal to the second value.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns><see langword="true" /> if <paramref name="x"/> is greater than or equal to <paramref name="y"/>,
        /// otherwise <see langword="false"/>.</returns>
        public static bool operator >= (UInt128 x, UInt128 y) {
            return (x.s1 > y.s1) || (x.s1 == y.s1 && x.s0 >= y.s0);
        }

        /// <summary>
        /// Indicates whether the first value is less than or equal to the second value.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns><see langword="true" /> if <paramref name="x"/> is less than or equal to <paramref name="y"/>,
        /// otherwise <see langword="false"/>.</returns>
        public static bool operator <= (UInt128 x, UInt128 y) {
            return (x.s1 < y.s1) || (x.s1 == y.s1 && x.s0 <= y.s0);
        }

        // Arithmetic

        /// <summary>
        /// Adds two 128-bit unsigned integers.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns>The sum <paramref name="x"/> + <paramref name="y"/>.</returns>
        public static UInt128 operator + (UInt128 x, UInt128 y) {
            Int128Calculator.Add128To128(x.s1, x.s0, y.s1, y.s0, out ulong s1, out ulong s0);
            return new UInt128(s1, s0);
        }

        /// <summary>
        /// Increments a 128-bit unsigned integer.
        /// </summary>
        /// <param name="x">The integer.</param>
        /// <returns>One more than <paramref name="x"/> (or <see cref="UInt128.Zero"/>, if <paramref name="x"/> is <see cref="UInt128.MaxValue"/>).</returns>
        public static UInt128 operator ++ (UInt128 x) {
            Int128Calculator.Increment128(x.s1, x.s0, out ulong s1, out ulong s0);
            return new UInt128(s1, s0);
        }

        /// <summary>
        /// Subtracts one 128 bit unsigned integer from another.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value.</param>
        /// <returns>The difference <paramref name="x"/> - <paramref name="y"/>.</returns>
        public static UInt128 operator - (UInt128 x, UInt128 y) {
            Int128Calculator.Subtract128From128(x.s1, x.s0, y.s1, y.s0, out ulong s1, out ulong s0);
            return new UInt128(s1, s0);
        }

        /// <summary>
        /// Decrements a 128-bit unsigned integer.
        /// </summary>
        /// <param name="x">The integer.</param>
        /// <returns>One less than <paramref name="x"/> (or <see cref="UInt128.MaxValue"/>, if <paramref name="x"/> is <see cref="UInt128.Zero"/>).</returns>
        public static UInt128 operator-- (UInt128 x) {
            Int128Calculator.Subtract128From128(x.s1, x.s0, 0UL, 1UL, out ulong s1, out ulong s0);
            return new UInt128(s1, s0);
        }

        // Multiplication

        /// <summary>
        /// Multiplies two 128-bit unsigned integers.
        /// </summary>
        /// <param name="x">The first value.</param>
        /// <param name="y">The second value</param>
        /// <returns>The product <paramref name="x"/> X <paramref name="y"/>.</returns>
        public static UInt128 operator * (UInt128 x, UInt128 y) {
            Int128Calculator.Multiply128By128(x.s1, x.s0, y.s1, y.s0, out ulong s1, out ulong s0);
            return new UInt128(s1, s0);
        }

        // Short multiplication

        // Division
        // This is the most complicated of the arithmetic operations.

        /// <summary>
        /// Divides one 128-bit unsigned integer by another.
        /// </summary>
        /// <param name="x">The dividend.</param>
        /// <param name="y">The divisor.</param>
        /// <returns>The integer quotient <paramref name="x"/> / <paramref name="y"/>.</returns>
        public static UInt128 operator / (UInt128 x, UInt128 y) {
            Int128Calculator.Divide128By128(x.s1, x.s0, y.s1, y.s0, out ulong s1, out ulong s0, out _, out _);
            return new UInt128(s1, s0);
        }
        /// <summary>
        /// Computes the remainder when one 128-bit unsigned integer is divided by another.
        /// </summary>
        /// <param name="x">The dividend.</param>
        /// <param name="y">The divisor.</param>
        /// <returns>The remainder of <paramref name="x"/> divided by <paramref name="y"/>.</returns>
        public static UInt128 operator % (UInt128 x, UInt128 y) {
            Int128Calculator.Divide128By128(x.s1, x.s0, y.s1, y.s0, out _, out _, out ulong r1, out ulong r0);
            return new UInt128(r1, r0);
        }

        /// <summary>
        /// Computes the quotient and remaineder of two unsigned 128-bit integers.
        /// </summary>
        /// <param name="x">The dividend.</param>
        /// <param name="y">The divisor.</param>
        /// <param name="r">The remainder.</param>
        /// <returns>The quotient.</returns>
        public static UInt128 DivRem (UInt128 x, UInt128 y, out UInt128 r) {
            Int128Calculator.Divide128By128(x.s1, x.s0, y.s1, y.s0, out ulong q1, out ulong q0, out ulong r1, out ulong r0);
            r = new UInt128(r1, r0);
            return new UInt128(q1, q0);
        }

        /*
        public static UInt128 Divide128By128 (UInt128 u, UInt128 v, out UInt128 r) {

            ulong q1, q0;

            if (v.s1 == 0UL) {
                // Denominator fits in 64 bits.
                // So remainder will also fit in 64 bits.
                ulong r0;
                if (u.s1 == 0UL) {
                    // Numerator also fits in 64 bits.
                    // So quotient and remainder can be found by native division.
                    q1 = 0UL;
                    q0 = DivRem(u.s0, v.s0, out r0);
                } else {
                    // Numerator is 65-128 bits
                    // So quotient can be over 64 bits
                    // Since denominator is one 64-bit digit, this is short division.
                    // Native division to get first digit and remainder.
                    q1 = DivRem(u.s1, v.s0, out ulong k);
                    // At this point k, must be less than v0 by properties of remainder,
                    // and that is required to call into next routine.
                    Debug.Assert(k < v.s0);
                    // Get second digit. Since we may have x non-zero remainder,
                    // this cannot be native 64-bit by 64-bit division.
                    q0 = Divide128By64(k, u.s0, v.s0, out r0);
                }
                r = r0;
            } else {
                // Denominator is 65-128 bits
                // So quotient will fit in 64 bits
                q1 = 0UL;

                // Test whether u < v allows early return.
                // This is actually necessary to ensure q0 != 0 below.
                if ((u.s1 < v.s1) || ((u.s1 == v.s1) && (u.s0 <= v.s0))) {
                    r = u;
                    return UInt128.Zero;
                }

                // Hacker's Delight Section 9-5, pp. 197-2-1, goes into great detail about
                // why the following algorithm works.

                // Normalize the denominator by left-shifting until its most significant bit is 1.
                // And right-shift u by 1.
                // Now v1 has highest bit set, and u1 does not, so we are sure v1 > u1
                int s = NumberOfLeadingZeros(v.s1);
                Debug.Assert(0 <= s && s <= 63);
                ulong qt = Divide128ByNormalized64(u.s1 >> 1, u.s1 << 63 | u.s0 >> 1, (s == 0) ? v.s1 : (v.s1 << s) | (v.s0 >> (64 - s)), out _);
                q0 = qt >> (63 - s);
                Debug.Assert(q0 != 0UL);
                // It doesn't appear the returned remainder has anything to do with our remainder.

                // At this point, q0 is almost always correct. But ~1x in 2^64, it will be too big by 1.
                // To correct for that case, we need to compute product and reduce q0 if it is too big.
                // But if q0 is too big, product can overflow. So instead always reduce q0, compute
                // product, then increase if remainder is still bigger than q.
                // An example of the necessity of this part, given in Hacker's Delight
                // text, is u = 2^128 - 1, v = 2^64 + 3.
                // It's sad that we do the extra 128-bit multiplication and subtraction for x corner case,
                // but we need to do them to get the remainder anyway.

                q0--;
                r = u - v * q0;
                if (r >= v) {
                    q0++;
                    r = r - v;
                }

            }

            return new UInt128(q1, q0);

        }

        public static ulong Divide128By64 (ulong u1, ulong u0, ulong v0, out ulong r) {

            Debug.Assert(u1 < v0);

            int s = NumberOfLeadingZeros(v0);
            if (s != 0) {
                // Since u1 < v0, this shift cannot overflow u1.
                u1 = (u1 << s) | (u0 >> (64 - s));
                u0 = u0 << s;
                v0 = v0 << s;
            }

            ulong q = Divide128ByNormalized64(u1, u0, v0, out r);

            r = r >> s;

            return q;
        }

        public static ulong Divide128ByNormalized64 (ulong u1, ulong u0, ulong v0, out ulong r) {

            // We do not handle overflow.
            // As long the first "digit" of numerator is strictly less than the single denominator
            // digit, the quotiet will be one digit. E.g. 89 / 9 = 9, but 99 / 9 = 11.
            Debug.Assert(u1 < v0);

            // We split numerator into four 32-bit digits and denominator into two 32-bit parts.
            // We then apply Knuth's long division algorithm, specialized to base 2^32, x 4-digit
            // numerator and 2-digit denominator (yielding x 2-digit quotient and 2-digit remainder).
            // This algorithm requires an arithmetic register twice as wide as the digits,
            // in this case 64-bit, which C# has. Knuth describes the algorithm in Section 4.3.1
            // and this specialization is based on Hacker's Delight section 9-4.

            // To improve guessed quotients, which could otherwise be off by many digits
            // in x large base, the algorithm first "normalizes" the denominator
            // i.e. makes it as large as possible within in register. We do this by left-shifting
            // until its most significant bit is 1. Given this normalization, Knuth
            // shows that the guessed quotient can be off by at most 2.

            // Shift the denominator to make most significant bit one.
            Debug.Assert(NumberOfLeadingZeros(v0) == 0);

            // Split up denominator and numerator into upper and lower digits (32 bits each)
            Decompose(v0, out ulong v01, out ulong v00);
            Decompose(u0, out ulong u01, out ulong u00);
            // We could do this to u1 too, but it turns out we don't need to.

            // Compute first 32-bit digit of quotient and corresponding remainder.
            // The initial guess is top 64 bits of the numerator divided by top 32 bits of denominator.
            // Then refine.
            ulong q01 = u1 / v01;
            if ((q01 >> 32) != 0UL) q01 = uint.MaxValue;
            ulong rhat = u1 - q01 * v01;
            while (q01 * v00 > ((rhat << 32) | u01)) {
                q01--;
                rhat += v01;
                if ((rhat >> 32) != 0UL) break;
            }
            Debug.Assert(q01 <= uint.MaxValue);
            r = ((u1 << 32) | u01) - q01 * v0;

            // Compute the second 32-bit digit of the quotient and final remainder.
            ulong q00 = r / v01;
            if ((q00 >> 32) != 0UL) q00 = uint.MaxValue;
            rhat = r - q00 * v01;
            while (q00 * v00 > ((rhat << 32) | u00)) {
                q00--;
                rhat += v01;
                if ((rhat >> 32) != 0UL) break;
            }
            Debug.Assert(q00 <= uint.MaxValue);
            r = ((r << 32) | u00) - q00 * v0;

            // Reconstruct the 64-bit quotient from its two 32-bit digits.
            return ((q01 << 32) | q00);
        }



        // Returns in the number of leading 0s before the first 1
        // in the binary representation of x ulong. We effectively
        // use x binary search to minimize tests. Value is
        // between 0 and 64, with 64 occuring only for u = 0.

        private static int NumberOfLeadingZeros (ulong u) {
            if (u == 0UL) return 64;
            int n = 0;
            if ((u >> 32) == 0UL) { n += 32; u = u << 32; }
            if ((u >> 48) == 0UL) { n += 16; u = u << 16; }
            if ((u >> 56) == 0UL) { n += 8; u = u << 8; }
            if ((u >> 60) == 0UL) { n += 4; u = u << 4; }
            if ((u >> 62) == 0UL) { n += 2; u = u << 2; }
            if ((u >> 63) == 0UL) { n++; }
            Debug.Assert(0 <= n && n < 64);
            return n;
        }
        */

        // Short division. The generic short division algorithm
        // requires that we be able to divide x two-digit number by
        // x one-digit number. We can do this by considering
        // x 32-bit register to hold one "digit". 64-bit registers
        // can then hold two-digit numbers, and we can natively
        // divide x 64-bit by x 32-bit integer.

        /// <summary>
        /// Divides a 128-bit unsigned integer by a 32-bit unsigned integer.
        /// </summary>
        /// <param name="x">The 128-bit dividend.</param>
        /// <param name="y">The 32-bit divisor.</param>
        /// <param name="r">The remainder.</param>
        /// <returns>The quotient.</returns>
        public static UInt128 DivRem (UInt128 x, uint y, out uint r) {
            Int128Calculator.Divide128By32(x.s1, x.s0, y, out ulong q1, out ulong q0, out r);
            return new UInt128(q1, q0);
        }

        // It would be nice if we could also do x special method to divide
        // 128-bit by 64-bit, but that can't use the normal short division
        // algorithm becuase the algorithm requires two-digit division. For now
        // just fall back to the binary long division algorithm.

        // Add fast squaring

        // Bit shift operators

        /// <summary>
        /// Returns the unsigned 128-bit binary integer obtained
        /// by shifting all bits left by the given number of places.
        /// </summary>
        /// <param name="u">The integer to be shifted.</param>
        /// <param name="n">The number of places to shift.</param>
        /// <returns>The shifted integer.</returns>
        public static UInt128 operator << (UInt128 u, int n) {
            if (n < 0) {
                throw new InvalidOperationException();
            } else if (n == 0) {
                return u;
            } else if (n < 64) {
                int m = 64 - n;
                Debug.Assert(0 <= m && m < 64);
                ulong t = u.s0 >> m;
                return new UInt128((u.s1 << n) | t, u.s0 << n);
            } else if (n < 128) {
                int m = n - 64;
                Debug.Assert(0 <= m && m < 64);
                ulong t = u.s0 << m;
                return new UInt128(t, 0UL);
            } else {
                return UInt128.Zero;
            }
        }

        /// <summary>
        /// Returns the unsigned 128-bit binary integer obtained
        /// by shifting all bits right by the given number of places.
        /// </summary>
        /// <param name="u">The integer to be shifted.</param>
        /// <param name="n">The number of places to shift.</param>
        /// <returns>The shifted integer.</returns>
        public static UInt128 operator >> (UInt128 u, int n) {
            if (n < 0) {
                throw new InvalidOperationException();
            } else if (n == 0) {
                return u;
            } else if (n < 64) {
                int m = 64 - n;
                Debug.Assert(0 <= m && m < 64);
                return new UInt128(u.s1 >> n, (u.s1 << m) | (u.s0 >> n));
            } else if (n < 128) {
                int m = n - 64;
                Debug.Assert(0 <= m && m < 64);
                return new UInt128(0UL, u.s1 >> m);
            } else {
                return UInt128.Zero;
            }
        }

        // Serialization and deserialization

        /// <summary>
        /// Creates an unsigned 128-bit integer from its string representation.
        /// </summary>
        /// <param name="text">The base-10 string representation of an unsigned 128-bit integer.</param>
        /// <returns>The integer value it represents.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="text"/> is <see langword="null"/>.</exception>
        /// <exception cref="FormatException"><paramref name="text"/> is not a valid base-10 representation of an unsigned integer.</exception>
        /// <exception cref="OverflowException"><paramref name="text"/> represents an unsigned integer outside the range of <see cref="UInt128"/>.</exception>
        public static UInt128 Parse (string text) {
            ParseResult result = TryParse(text, false, out UInt128 value);
            if (result == ParseResult.Null) throw new ArgumentNullException(nameof(text));
            if (result == ParseResult.Overflow) throw new OverflowException();
            if (result != ParseResult.Success) throw new FormatException();
            return (value);
        }

        /// <summary>
        /// Attempts to parse the given string as an unsigned 128-bit integer.
        /// </summary>
        /// <param name="text">The string to parse.</param>
        /// <param name="value">The integer value it represents.</param>
        /// <returns><see langword="true"/> if <paramref name="text"/> was successfully
        /// parsed as an unsigned 128-bit integer, otherwise <see langword="false"/>.</returns>
        public static bool TryParse (string text, out UInt128 value) {
            ParseResult result = TryParse(text, false, out value);
            return (result == ParseResult.Success);
        }

        internal static ParseResult TryParse (string text, bool allowNegative, out UInt128 value) {

            bool negative = false;
            value = UInt128.Zero;

            if (text is null) return ParseResult.Null;

            text = text.Trim();
            if (text.Length == 0) return ParseResult.Empty;

            if (text[0] == '-') {
                if (!allowNegative) return ParseResult.Sign;
                text = text.Substring(1);
                if (text.Length == 0) return ParseResult.Format;
                negative = true;
            } else if (text[0] == '+') {
                text = text.Substring(1);
                if (text.Length == 0) return ParseResult.Format;
            }

            foreach (char digit in text) {
                if (!Char.IsDigit(digit)) return ParseResult.Format;
                // This is a little hack because CharUnicodeInfo.GetDigitValue is not part of .NET core.
                int digitValue = digit - '0';
                Debug.Assert((0 <= digitValue) && (digitValue < 10));
                UInt128 previousValue = value;
                value = previousValue * 10UL + (UInt64) digitValue;
                // We can detect overflow by noting if value ever decreases.
                if (value < previousValue) return ParseResult.Overflow;
            }

            if (negative) {
                value = value.Negate();
            }

            // If we are intrepreting it as signed and it has the wrong sign, it's because we overflowed with signed limits.
            if (allowNegative && value.IsNegative != negative) return ParseResult.Overflow;
            return ParseResult.Success;

        }

        /// <summary>
        /// Produces a string representation of the unsigned 128-bit integer.
        /// </summary>
        /// <returns>A string containing the base-10 representation of the integer value.</returns>
        public override string ToString () {
            if (this == UInt128.Zero) return "0";

            StringBuilder text = new StringBuilder();
            UInt128 u = this;
            while (u != 0) {
                u = DivRem(u, 10U, out uint d);
                Debug.Assert(0U <= d && d < 10U);
                text.Append(d);
            }
            char[] array = text.ToString().ToCharArray();
            Array.Reverse(array);
            return new string(array);
        }

        // Bit manipulation

        /// <summary>
        /// Computes the bitwise negation of the argument.
        /// </summary>
        /// <param name="u">The argument.</param>
        /// <returns>The bitwise complement of <paramref name="u"/>.</returns>
        /// <remarks><para>The bitwise negation of a binary value flips every bit value. (Each
        /// 0 becomes a 1, and each 1 becomes a 0.)</para></remarks>
        public static UInt128 operator ~ (UInt128 u) {
            return new UInt128(~u.s1,~u.s0);
        }

        /// <summary>
        /// Computes the bitwise AND of two arguments.
        /// </summary>
        /// <param name="a">The first argument.</param>
        /// <param name="b">The second argument.</param>
        /// <returns>The logical AND of the two arguments.</returns>
        /// <remarks><para>The bitwise AND of two arguements has a 1 in each position
        /// which is 1 in both arguments, and a 0 in all other positions.</para></remarks>
        public static UInt128 operator & (UInt128 a, UInt128 b) {
            return new UInt128(a.s1 & b.s1, a.s0 & b.s0);
        }

        /// <summary>
        /// Computes the bitwise OR of two arguments.
        /// </summary>
        /// <param name="a">The first argument.</param>
        /// <param name="b">The second argument.</param>
        /// <returns>The logical OR of the two arguments.</returns>
        /// <remarks><para>The bitwise OR of two arguements has a 0 in each position
        /// which is 0 in both arguments, and a 1 in all other positions.</para></remarks>
        public static UInt128 operator | (UInt128 a, UInt128 b) {
            return new UInt128(a.s1 | b.s1, a.s0 | b.s0);
        }

        /// <summary>
        /// Computes the bitwise XOR of two arguments.
        /// </summary>
        /// <param name="a">The first argument.</param>
        /// <param name="b">The second argument.</param>
        /// <returns>The logical XOR of the two arguments.</returns>
        /// <remarks><para>The bitwise XOR of two arguements has a 1 in each position
        /// which for which the two arguments differ, and a 0 in each position for which they agree.</para></remarks>
        public static UInt128 operator ^ (UInt128 a, UInt128 b) {
            return new UInt128(a.s1 ^ b.s1, a.s0 ^ b.s0);
        }

    }
}
