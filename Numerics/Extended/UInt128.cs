using System;
using System.Diagnostics;
using System.Globalization;

namespace Meta.Numerics.Extended
{

    [CLSCompliant(false)]
    public struct UInt128B : IEquatable<UInt128B>, IComparable<UInt128B> {

        // Originally, I implemented this with two ulongs instead of four uints.
        // Addition was about 2x faster, but multiplication was about 10x slower,
        // because I needed to decompose ulongs into uint parts in order to do
        // multiplication within a ulong register. May want to revisit this
        // to see if I can use both uints and ulongs in optimal places. 

        public UInt128B (uint s3, uint s2, uint s1, uint s0) {
            this.s3 = s3;
            this.s2 = s2;
            this.s1 = s1;
            this.s0 = s0;
        }

        private readonly uint s3, s2, s1, s0;

        // Special instances

        public static readonly UInt128B Zero = new UInt128B(0, 0, 0, 0);

        public static readonly UInt128B One = new UInt128B(0, 0, 0, 1);

        public static readonly UInt128B MinValue = Zero;

        public static readonly UInt128B MaxValue = new UInt128B(UInt32.MaxValue, UInt32.MaxValue, UInt32.MaxValue, UInt32.MaxValue);

        // Equality operators

        public static bool operator == (UInt128B a, UInt128B b) {
            return (Equals(a, b));
        }

        public static bool operator != (UInt128B a, UInt128B b) {
            return (!Equals(a, b));
        }

        // Equality implementations

        private static bool Equals (UInt128B a, UInt128B b) {
            return (a.s3 == b.s3 && a.s2 == b.s2 && a.s1 == b.s1 && a.s0 == b.s0);
        }

        public bool Equals (UInt128B other) {
            return (Equals(this, other));
        }

        public override bool Equals (object other) {
            if (other is UInt128B) {
                return (Equals(this, (UInt128B) other));
            } else {
                return (false);
            }
        }

        public override int GetHashCode () {
            return (s0.GetHashCode() + 5 * s1.GetHashCode() + 7 * s2.GetHashCode() + 13 * s3.GetHashCode());
        }

        // Comparison operators

        public static bool operator < (UInt128B a, UInt128B b) {
            return (Compare(a, b) < 0);
        }

        public static bool operator > (UInt128B a, UInt128B b) {
            return (Compare(a, b) > 0);
        }

        public static bool operator <= (UInt128B a, UInt128B b) {
            return (Compare(a, b) <= 0);
        }

        public static bool operator >= (UInt128B a, UInt128B b) {
            return (Compare(a, b) >= 0);
        }

        // Comparison implementations

        private static int Compare (UInt128B a, UInt128B b) {
            int result = a.s3.CompareTo(b.s3);
            if (result != 0) return (result);
            result = a.s2.CompareTo(b.s2);
            if (result != 0) return (result);
            result = a.s1.CompareTo(b.s1);
            if (result != 0) return (result);
            result = a.s0.CompareTo(b.s0);
            return (result);
        }

        public int CompareTo (UInt128B other) {
            return (Compare(this, other));
        }

        // Arithmetic operators

        public static UInt128B operator + (UInt128B a, UInt128B b) {
            return (Add(a, b));
        }

        public static UInt128B operator - (UInt128B a, UInt128B b) {
            return (Subtract(a, b));
        }

        public static UInt128B operator * (UInt128B a, UInt128B b) {
            return (Multiply(a, b));
        }

        public static UInt128B operator / (UInt128B a, UInt128B b) {
            UInt128B r;
            UInt128B q = DivRem(a, b, out r);
            return (q);
        }

        public static UInt128B operator % (UInt128B a, UInt128B b) {
            UInt128B r;
            UInt128B q = DivRem(a, b, out r);
            return (r);
        }

        // Arithmetic operator implementations

        public static UInt128B Add (UInt128B a, UInt128B b) {
            uint carry = 0;
            uint s0 = AddWithCarry(a.s0, b.s0, ref carry);
            uint s1 = AddWithCarry(a.s1, b.s1, ref carry);
            uint s2 = AddWithCarry(a.s2, b.s2, ref carry);
            uint s3 = AddWithCarry(a.s3, b.s3, ref carry);
            return (new UInt128B(s3, s2, s1, s0));
        }

        private static uint AddWithCarry (uint a, uint b, ref uint c) {
            ulong sum = (ulong) a + b + c;
            c = (uint) (sum >> 32);
            return ((uint) sum);
            // I also implemented this with addition in a uint register
            // and a comparison test for overflow and the speed was
            // essentially the same.
        }

        public static UInt128B Subtract (UInt128B a, UInt128B b) {
            uint borrow = 0;
            uint s0 = SubtractWithBorrow(a.s0, b.s0, ref borrow);
            uint s1 = SubtractWithBorrow(a.s1, b.s1, ref borrow);
            uint s2 = SubtractWithBorrow(a.s2, b.s2, ref borrow);
            uint s3 = SubtractWithBorrow(a.s3, b.s3, ref borrow);
            return (new UInt128B(s3, s2, s1, s0));
        }

        private static uint SubtractWithBorrow(uint a, uint b, ref uint c) {
            ulong difference = (ulong) a - b - c;
            c = (uint) - (int) (difference >> 32);
            return ((uint) difference);
        }

        public static UInt128B Multiply (UInt128B a, UInt128B b) {

            uint carry = 0;
            uint s0 = MultiplyWithCarry(0, a.s0, b.s0, ref carry);
            uint s1 = MultiplyWithCarry(0, a.s1, b.s0, ref carry);
            uint s2 = MultiplyWithCarry(0, a.s2, b.s0, ref carry);
            uint s3 = MultiplyWithCarry(0, a.s3, b.s0, ref carry);

            carry = 0;
            s1 = MultiplyWithCarry(s1, a.s0, b.s1, ref carry);
            s2 = MultiplyWithCarry(s2, a.s1, b.s1, ref carry);
            s3 = MultiplyWithCarry(s3, a.s2, b.s1, ref carry);

            carry = 0;
            s2 = MultiplyWithCarry(s2, a.s0, b.s2, ref carry);
            s3 = MultiplyWithCarry(s3, a.s1, b.s2, ref carry);

            carry = 0;
            s3 = MultiplyWithCarry(s3, a.s0, b.s3, ref carry);

            return (new UInt128B(s3, s2, s1, s0));
        }

        private static UInt128B Multiply(UInt128B a, uint b) {

            uint carry = 0;
            uint s0 = MultiplyWithCarry(0, a.s0, b, ref carry);
            uint s1 = MultiplyWithCarry(0, a.s1, b, ref carry);
            uint s2 = MultiplyWithCarry(0, a.s2, b, ref carry);
            uint s3 = MultiplyWithCarry(0, a.s3, b, ref carry);

            return (new UInt128B(s3, s2, s1, s0));

        }

        private static uint MultiplyWithCarry (uint start, uint a, uint b, ref uint c) {
            ulong sum = (ulong) start + (ulong) a * b + (ulong) c;
            c = (uint) (sum >> 32);
            // Unexpectedly to me, this cast doesn't throw when sum > unit.MaxValue, but
            // instead throws away the high bits, which is what we want here. So its
            // equivilent to (sum << 32) >> 32, but without the bit shifting.
            uint d = (uint) sum;
            return (d);
            // Make a version with sum = 0, since that is the case in about half the calls
            // into this method.
        }

        public static UInt128B DivRem (UInt128B a, uint b, out uint r) {
            ulong remainder = 0UL;
            uint s3 = (uint) DivRem(a.s3, b, out remainder);
            uint s2 = (uint) DivRem((remainder << 32) + a.s2, b, out remainder);
            uint s1 = (uint) DivRem((remainder << 32) + a.s1, b, out remainder);
            uint s0 = (uint) DivRem((remainder << 32) + a.s0, b, out remainder);
            r = (uint) remainder;
            return (new UInt128B(s3, s2, s1, s0));
        }

        private static ulong DivRem (ulong a, uint b, out ulong r) {
            r = a % b;
            return ( a / b);
        }

        public static UInt128B AltDivRem (UInt128B a, UInt128B b, out UInt128B r) {

            // Find most significant part of denominator
            uint d = b.s3;
            if (d == 0) {
                d = b.s2;
                if (d == 0) {
                    d = b.s1;
                    if (d == 0) {
                        d = b.s0;
                        if (d == 0) throw new DivideByZeroException();
                    }
                }
            }

            // Guestimate each digit in turn
            UInt128B n = a;
            uint q3 = 0;
            if (n.s3 != 0) {
                ulong m = ((ulong) n.s3) | n.s2;
                uint c = (uint) DivRem(m, d, out _);
                UInt128B nh = Multiply(b, c);
                // Guaranteed to be almost right
                n = n - nh;
                q3 = c;
            }

            uint q2 = 0;
            if (n.s2 != 0) {
                ulong m = ((ulong) n.s2) | n.s1;
                uint c = (uint) DivRem(m, d, out _);
            }


            r = UInt128B.Zero;
            return (new UInt128B(q3, q2, 0, 0));
        }

        public static UInt128B DivRem (UInt128B a, UInt128B b, out UInt128B r) {

            if (a == UInt128B.Zero) throw new DivideByZeroException();

            if (a <= b) {
                if (a == b) {
                    r = UInt128B.Zero;
                    return (UInt128B.One);
                } else {
                    r = a;
                    return (UInt128B.Zero);
                }
            }

            // At this point, we know d >=1 and n > d, so also n > 1.

            //if (n.s1 == 0UL) {
            //    Debug.Assert(d.s1 == 0UL);
            //    r = new UInt128(0UL, n.s0 % d.s0);
            //    return (new UInt128(0UL, n.s0 / d.s0));
            //}

            // Find the highest bit
            int n = 127;
            while (!GetBit(a, n)) n--;

            // Binary long division algorithm
            UInt128B q = UInt128B.Zero;
            r = UInt128B.Zero;
            for (int i = n; i >= 0; i--) {
                r = LeftShiftOne(r);
                r = SetBit(r, 0, GetBit(a, i));
                if (r >= b) {
                    r = r - b;
                    q = SetBit(q, i, true);
                }
            }

            return (q);
        }

        private static UInt128B LeftShiftOne (UInt128B v) {
            uint s0 = v.s0 << 1;
            uint s1 = (v.s1 << 1) | (v.s0 >> 31);
            uint s2 = (v.s2 << 1) | (v.s1 >> 31);
            uint s3 = (v.s3 << 1) | (v.s2 >> 31);
            return (new UInt128B(s3, s2, s1, s0));
        }

        private static bool GetBit (UInt128B v, int n) {
            if (n < 32) {
                return ((v.s0 & (1 << n)) != 0);
            } else if (n < 2 * 32) {
                return ((v.s1 & (1 << (n - 32))) != 0);
            } else if (n < 3 * 32) {
                return ((v.s2 & (1 << (n - 2 * 32))) != 0);
            } else {
                return ((v.s3 & (1 << (n - 3 * 32))) != 0);
            }
        }

        private static UInt128B SetBit (UInt128B v, int n, bool bit) {
            Debug.Assert(0 <= n && n < 128);
            if (n < 32) {
                return (new UInt128B(v.s3, v.s2, v.s1, SetBit(v.s0, n, bit)));
            } else if (n < 2 * 32) {
                return (new UInt128B(v.s3, v.s2, SetBit(v.s1, n - 32, bit), v.s0));
            } else if (n < 3 * 32) {
                return (new UInt128B(v.s3, SetBit(v.s2, n - 2 * 32, bit), v.s1, v.s0));
            } else {
                return (new UInt128B(SetBit(v.s3, n - 3 * 32, bit), v.s2, v.s1, v.s0));
            }
        }

        private static uint SetBit (uint v, int n, bool value) {
            Debug.Assert(0 <= n && n < 32);
            uint mask = ~(1U << n);
            uint result = v & mask;
            if (value) result = result | (1U << n);
            return (result);
        }

        // Conversion casts

        public static implicit operator UInt128B (uint value) {
            return (new UInt128B(0, 0, 0, value));
        }

        public static explicit operator uint (UInt128B value) {
            return (value.s0);
        }

        public static implicit operator UInt128B (ulong value) {
            return (new UInt128B(0, 0, (uint) value >> 32, (uint) value));
        }

        public static explicit operator ulong (UInt128B value) {
            return ((((ulong) value.s1) << 32) + value.s0);
        }

        // Printing and parsing

        public static UInt128B Parse (string text) {
            UInt128B value;
            bool success = TryParse(text, out value);
            if (!success) throw new FormatException();
            return (value);
        }

        public static bool TryParse (string text, out UInt128B value) {
            value = UInt128B.Zero;
            if (String.IsNullOrEmpty(text)) return (false);
            foreach (char digit in text) {
                if (!Char.IsDigit(digit)) return (false);
                // This is a little hack because CharUnicodeInfo.GetDigitValue is not part of .NET core.
                int digitValue = digit - '0';
                Debug.Assert((0 <= digitValue) && (digitValue < 10));
                value = Multiply(value, 10) + (uint) digitValue;
            }
            return (true);
        }

        private static string ToString (UInt128B value) {
            char[] buffer = new char[42];
            int i = buffer.Length - 1;
            while (i >= 0) {
                uint digit;
                value = DivRem(value, 10, out digit);
                buffer[i] = (char) ((uint) '0' + digit);
                if (value == UInt128B.Zero) break;
                i--;
            }
            return (new String(buffer, i, buffer.Length - i));
        }

        public override string ToString () {
            return (ToString(this));
        }

    }

    /// <summary>
    /// Represents a 128-bit unsigned integer.
    /// </summary>
    /// <remarks>
    /// <para>The range of this structure is 0 to approximately 3.4E38.</para>
    /// <para>If you know that the integers you need to work with fit within this range, operations using
    /// this structure will be faster than using <see cref="System.Numerics.BigInteger"/>.
    /// (Addition and subtraction are about 5x faster. Multiplication is about 1x faster.)
    /// </para>
    /// </remarks>
    [CLSCompliant(false)]
    public struct UInt128 : IEquatable<UInt128>, IComparable<UInt128> {

        internal UInt128 (ulong s1, ulong s0) {
            this.s0 = s0;
            this.s1 = s1;
        }

        // s = s_1 2^64 + s_0

        private readonly ulong s0;
        private readonly ulong s1;

        // Special cases

        /// <summary>
        /// The zero value of the unsigned 128-bit integer.
        /// </summary>
        public static readonly UInt128 Zero = new UInt128(0UL, 0UL);

        public static readonly UInt128 One = new UInt128(0UL, 1UL);

        /// <summary>
        /// The least possible value of the unsigned 128-bit integer.
        /// </summary>
        public static readonly UInt128 MinValue = new UInt128(0UL, 0UL);

        /// <summary>
        /// The greatest possible value of the unsigned 128-bit integer.
        /// </summary>
        public static readonly UInt128 MaxValue = new UInt128(UInt64.MaxValue, UInt64.MaxValue);

        // Casts

        /// <summary>
        /// Creates an unsigned 128-bit integer from an unsigned 64-bit integer.
        /// </summary>
        /// <param name="c">The unsigned 64-bit integer.</param>
        /// <remarks><para>This is a widening cast, which cannot fail.</para></remarks>
        public static implicit operator UInt128 (UInt64 c) {
            return (new UInt128(0UL, c));
        }

        /// <summary>
        /// Creates an unsigned 64-bit integer from an unsigned 128-bit integer.
        /// </summary>
        /// <param name="c">The unsigned 128-bit integer.</param>
        /// <remarks><para>This is a narrowing cast, which will fail if the value of
        /// the 128-bit integer is greater than <see cref="UInt64.MaxValue"/>.</para></remarks>
        public static explicit operator UInt64 (UInt128 c) {
            if (c.s1 == 0UL) {
                return (c.s0);
            } else {
                throw new InvalidCastException();
            }
        }

        /// <summary>
        /// Creates an double from an unsigned 128-bit integer.
        /// </summary>
        /// <param name="c">The unsigned 128-bit integer.</param>
        /// <remarks><para>This cast will never fail, but it will loose precision
        /// for values larger than about 1.0E16, because <see cref="System.Double"/>
        /// maintains fewer than 128 bits of precision.</para></remarks>
        public static implicit operator double (UInt128 c) {
            const double multiplier = (double) UInt64.MaxValue + 1.0;
            return ((double) c.s0 + (double) c.s1 * multiplier);
        }

        // Equality

        public static bool Equals (UInt128 a, UInt128 b) {
            return ((a.s0 == b.s0) && (a.s1 == b.s1));
        }

        public static bool operator == (UInt128 a, UInt128 b) {
            return (Equals(a, b));
        }

        public static bool operator != (UInt128 a, UInt128 b) {
            return (!Equals(a, b));
        }

        public bool Equals (UInt128 b) {
            return (Equals(this, b));
        }

        public override bool Equals (object obj) {
            if (obj is UInt128) {
                return (Equals(this, (UInt128) obj));
            } else {
                return (false);
            }
        }

        public override int GetHashCode () {
            unchecked {
                return (s0.GetHashCode() + 13 * s1.GetHashCode());
            }
        }

        // Comparable

        private static int Compare (UInt128 a, UInt128 b) {
            int c1 = a.s1.CompareTo(b.s1);
            if (c1 == 0) {
                Debug.Assert(a.s1 == b.s1);
                int c0 = a.s0.CompareTo(b.s0);
                return (c0);
            } else {
                return (c1);
            }
        }

        public int CompareTo (UInt128 b) {
            return (Compare(this, b));
        }

        /// <summary>
        /// Returns a value that indicates whether the first value is greater than the second value.
        /// </summary>
        /// <param name="a">The first value.</param>
        /// <param name="b">The second value.</param>
        /// <returns><see langword="true" /> if <paramref name="a"/> is greater than <paramref name="b"/>,
        /// otherwise <see langword="false"/>.</returns>
        public static bool operator > (UInt128 a, UInt128 b) {
            int comparison = Compare(a, b);
            return (comparison > 0);
        }

        /// <summary>
        /// Returns a value that indicates whether the first value is less than the second value.
        /// </summary>
        /// <param name="a">The first value.</param>
        /// <param name="b">The second value.</param>
        /// <returns><see langword="true" /> if <paramref name="a"/> is less than <paramref name="b"/>,
        /// otherwise <see langword="false"/>.</returns>
        public static bool operator < (UInt128 a, UInt128 b) {
            int comparison = Compare(a, b);
            return (comparison < 0);
        }

        /// <summary>
        /// Returns a value that indicates whether the first value is greater than or equal to the second value.
        /// </summary>
        /// <param name="a">The first value.</param>
        /// <param name="b">The second value.</param>
        /// <returns><see langword="true" /> if <paramref name="a"/> is greater than or equal to <paramref name="b"/>,
        /// otherwise <see langword="false"/>.</returns>
        public static bool operator >= (UInt128 a, UInt128 b) {
            int comparison = Compare(a, b);
            return (comparison >= 0);
        }

        /// <summary>
        /// Returns a value that indicates whether the first value is less than or equal to the second value.
        /// </summary>
        /// <param name="a">The first value.</param>
        /// <param name="b">The second value.</param>
        /// <returns><see langword="true" /> if <paramref name="a"/> is less than or equal to <paramref name="b"/>,
        /// otherwise <see langword="false"/>.</returns>
        public static bool operator <= (UInt128 a, UInt128 b) {
            int comparison = Compare(a, b);
            return (comparison <= 0);
        }

        // Arithmetic

        /// <summary>
        /// Adds two values.
        /// </summary>
        /// <param name="a">The first value.</param>
        /// <param name="b">The second value.</param>
        /// <returns>The sum of <paramref name="a"/> and <paramref name="b"/>.</returns>
        public static UInt128 operator + (UInt128 a, UInt128 b) {
            return (Add(a, b));
        }

        public static UInt128 operator - (UInt128 a, UInt128 b) {
            return (Subtract(a, b));
        }

        public static UInt128 operator * (UInt128 a, UInt128 b) {
            return (Multiply(a, b));
        }

        private static UInt128 Multiply (UInt128 a, UInt128 b) {
            if (a.s1 == 0UL) {
                if (b.s1 == 0UL) {
                    return (Multiply64(a.s0, b.s0));
                } else {
                    return (Multiply128By64(b, a.s0));
                }
            } else {
                if (b.s1 == 0UL) {
                    return (Multiply128By64(a, b.s0));
                } else {
                    return (Multiply128(a, b));
                }
            }
        }

        private static UInt128 Multiply128By64 (UInt128 a, UInt64 b) {

            // (a1 2^64 + a0) * b0
            //   = a1 b0 2^64 + a0 b0 = q 2^64 + p
            //   = (q1 2^64 + q0) 2^64 + (p1 2^64 + p0)
            //   = q1 2^128 + (q0 + p1) 2^64 + p0

            UInt128 p = Multiply64(a.s0, b);
            UInt128 q = Multiply64(a.s1, b);

            // We ignore q1, which is overflow.

            UInt64 r = q.s0 + p.s1;

            return (new UInt128(r, p.s0));
        }

        private static UInt128 Multiply128 (UInt128 a, UInt128 b) {

            // (a1 2^64 + a0) * (b1 2^64 + b0) =
            //   = (a1 b1) 2^128 + (a1 b0 + a0 b1) 2^64 + a0 b0
            // Define p = a0 b0, q = a1 b0 + a0 b1, r = a1 b1, and decompose.
            //   = (r1 2^64 + r0) 2^128 + (q1 2^64 + q0) 2^64 + (p1 2^64 + p0)
            //   = [r1 2^64 + (r0 + q1)] 2^128 + [(q0 + p1) 2^64 + p0]

            UInt128 p = Multiply64(a.s0, b.s0);
            UInt128 qa = Multiply64(a.s1, b.s0);
            UInt128 qb = Multiply64(a.s0, b.s1);
            UInt128 q = Add(qa, qb);

            UInt64 s0 = p.s0;
            UInt64 s1 = Add64(q.s0, p.s1, out _);

            return (new UInt128(s1, s0));
        }

        // Multiply two 64-bit ints to get a 128-bit result.
        // Overflow should be impossible.
        private static UInt128 Multiply64 (UInt64 a, UInt64 b) {
            checked {

                // Decompose into 32-bit parts because we can
                // multiply two of them within an existing 64-bit
                // register without overflow.

                // a = a_1 2^{32} + a_0
                ulong a1, a0;
                Decompose(a, out a1, out a0);

                // b = b_1 2^{32} + b_1
                ulong b1, b0;
                Decompose(b, out b1, out b0);

                // ab = a_1 b_1 2^{64} + a_0 b_0 +
                //      (a_1 b_0 + a_0 b_1) 2^{32}
                //    = s_1 2^{64} + (t_a + t_b) 2^{32} + s_0

                // Guaranteed no overflow of these multiplies
                UInt64 s0 = Multiply32(a0, b0);
                UInt64 ta = Multiply32(a1, b0);
                UInt64 tb = Multiply32(a0, b1);
                UInt64 s1 = Multiply32(a1, b1);

                // Decompose the cross term, noting the possibility of a carry bit.
                //   t = t_a + t_b = c 2^{64} + t_1 2^{32} + t_0
                // Seperate its contributions into the upper and lower 64 bits
                // of the result.
                //   t 2^{32} = (c 2^{64} + t_1 2^{32} + t_0) 2^32
                //            = (c 2^{32} + t_1) 2^{64} + t_0 2^{32}

                bool carry;
                UInt64 t = Add64(ta, tb, out carry);
                if (carry) s1 = Add64(s1, 1UL << 32, out carry);
                Debug.Assert(!carry);

                ulong t1, t0;
                Decompose(t, out t1, out t0);

                s0 = Add64(s0, t0 << 32, out carry);
                if (carry) s1++;

                s1 = Add64(s1, t1, out carry);
                Debug.Assert(!carry);

                return (new UInt128(s1, s0));
            }
        }

        // Multiply 2 32-bit ints (stored in 64-bit registers for convenience)
        // to get a 64-bit result. Overflow is impossible.

        private static UInt64 Multiply32 (UInt64 a, UInt64 b) {
            // Overflow behavior: If inputs are actually 32-bit, this
            // multiplication will not overflow the 64-bit register
            // in which it occurs. We do it in a checked context only as 
            // a double-check. QED
            Debug.Assert(a <= UInt32.MaxValue);
            Debug.Assert(b <= UInt32.MaxValue);
            checked {
                return (a * b);
            }
        }

        private static UInt128 Add (UInt128 a, UInt128 b) {
            // Overflow behavior: We make the low-bit sum unchecked so that
            // we can carry overflow into the high-bit summand.
            UInt64 s0;
            unchecked {
                s0 = a.s0 + b.s0;
            }
            bool carry = (s0 < a.s0) && (s0 < b.s0);
            UInt64 s1 = a.s1 + b.s1;
            if (carry) s1++;
            return (new UInt128(s1, s0));          
        }

        private static UInt64 Add64 (UInt64 a, UInt64 b, out bool carry) {
            UInt64 sum;
            unchecked {
                sum = a + b;
            }
            carry = (sum < a) && (sum < b);
            return (sum);
        }

        private static UInt128 Subtract (UInt128 a, UInt128 b) {
            UInt64 s0;
            bool borrow;
            unchecked {
                borrow = (a.s0 < b.s0);
                s0 = a.s0 - b.s0;
            }
            UInt64 s1 = a.s1 - b.s1;
            if (borrow) s1--;
            return (new UInt128(s1, s0));
        }

        public static UInt128 DivRem (UInt128 n, UInt128 d, out UInt128 r) {

            if (d == UInt128.Zero) throw new DivideByZeroException();

            if (n <= d) {
                if (n == d) {
                    r = UInt128.Zero;
                    return (UInt128.One);
                } else {
                    r = n;
                    return (UInt128.Zero);
                }
            }

            // At this point, we know d >=1 and n > d, so also n > 1.

            if (n.s1 == 0UL) {
                Debug.Assert(d.s1 == 0UL);
                r = new UInt128(0UL, n.s0 % d.s0);
                return (new UInt128(0UL, n.s0 / d.s0));
            }

            // Find the highest bit
            int b = 127;
            while (!GetBit(n, b)) b--;

            // Binary long division algorithm
            UInt128 q = UInt128.Zero;
            r = UInt128.Zero;
            for (int i = b; i >= 0; i--) {
                r = LeftShiftOne(r);
                r = SetBit(r, 0, GetBit(n, i));
                if (r >= d) {
                    r = r - d;
                    q = SetBit(q, i, true);
                }
            }

            return (q);
        }

        public static UInt128 operator / (UInt128 a, UInt128 b) {
            UInt128 q = DivRem(a, b, out _);
            return (q);
        }

        public static UInt128 operator % (UInt128 a, UInt128 b) {
            UInt128 r;
            UInt128 q = DivRem(a, b, out r);
            return (r);
        }

        // Decomposes a 64-bit number into two 32-bit parts (stored in 64-bit
        // registers for convenience). v = hi 2^32 + lo.

        private static void Decompose (UInt64 v, out UInt64 hi, out UInt64 lo) {
            hi = v >> 32;
            lo = v & 0xFFFFFFFF;
            Debug.Assert(hi <= UInt32.MaxValue);
            Debug.Assert(lo <= UInt32.MaxValue);
        }

        // Add fast mixed operations
        // Add fast squaring

        // Shift operators

        private static UInt128 LeftShiftOne (UInt128 v) {
            // Get the leftmost bit of the lower register.
            UInt64 b = v.s0 >> 63;
            Debug.Assert(b == 0UL || b == 1UL);
            // Left shift the lower register (which looses the bit we read).
            UInt64 s0 = v.s0 << 1;
            // Left shift the higher register, and set the lowest bit to
            // the one we read from the lower register.
            UInt64 s1 = (v.s1 << 1) | b;
            // Return the result.
            return (new UInt128(s1, s0));
        }

        internal static bool GetBit (UInt128 v, int n) {
            Debug.Assert(0 <= n && n < 128);
            if (n < 64) {
                return (GetBit(v.s0, n));
            } else {
                return (GetBit(v.s1, n - 64));
            }
        }

        internal static bool GetBit(UInt64 v, int n) {
            Debug.Assert(0 <= n && n < 64);
            UInt64 mask = 1UL << n;
            UInt64 result = v & mask;
            return (result > 0UL);
        }

        internal static UInt128 SetBit (UInt128 v, int n, bool value) {
            Debug.Assert(0 <= n && n < 128);
            if (n < 64) {
                return (new UInt128(v.s1, SetBit(v.s0, n, value)));
            } else {
                return (new UInt128(SetBit(v.s1, n - 64, value), v.s0));
            }
        }

        internal static UInt64 SetBit (UInt64 v, int n, bool value) {
            Debug.Assert(0 <= n && n < 64);
            UInt64 mask = ~(1UL << n);
            UInt64 result = v & mask;
            if (value) result = result | (1UL << n);
            return (result);
        }

        public static UInt128 Parse (string text) {
            UInt128 value;
            bool parsed = TryParse(text, out value);
            if (!parsed) throw new FormatException();
            return (value);
        }

        public static bool TryParse (string text, out UInt128 value) {
            value = UInt128.Zero;
            if (String.IsNullOrEmpty(text)) return (false);
            foreach (char digit in text) {
                if (!Char.IsDigit(digit)) return (false);
                // This is a little hack because CharUnicodeInfo.GetDigitValue is not part of .NET core.
                int digitValue = digit - '0';
                Debug.Assert((0 <= digitValue) && (digitValue < 10));
                value = value * 10UL + (UInt64) digitValue;
            }
            return (true);
        }

        private static string ToString (UInt128 value) {
            UInt128 ten = 10UL;
            char[] buffer = new char[42];
            int i = buffer.Length - 1;
            while (i >= 0) {
                UInt128 digit;
                value = DivRem(value, ten, out digit);
                buffer[i] = (char) ((int) '0' + (int) digit);
                if (value == UInt128.Zero) break;
                i--;
            }
            return (new String(buffer, i, buffer.Length - i));
        }

        public override string ToString () {
            return (ToString(this));
        }

    }
}
