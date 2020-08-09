using System;
using System.Diagnostics;

namespace Meta.Numerics.Extended {


    // These routines are adapted from
    // Warren, Henry, "Hacker's Delight", 2nd edition


    internal static class Int128Calculator {

        // Adds two 128 bit numbers.

        public static void Add128To128 (ulong x1, ulong x0, ulong y1, ulong y0, out ulong z1, out ulong z0) {
            unchecked {
                // The lower sum should be unchecked, since it is okay if it overflows even in checked mode.
                z0 = x0 + y0;
            }
            // We need to add a carry bit to the higher sum, if the lower sum overflowed.
            // A fast, easy test for this is whether sum is less than either summand
            // i.e. if (z0 < x0) z1++, but tests indicate that this trick from HD 2-16
            // to produce a 0 or 1 without comparison or branching is even 10% faster than that.
            ulong c = ((x0 & y0) | ((x0 | y0) & ~z0)) >> 63;
            // The higher sum should inherit checked or unchecked, so it throws or not depending on caller's mode.
            z1 = x1 + y1 + c;
        }

        public static void Subtract128From128 (ulong x1, ulong x0, ulong y1, ulong y0, out ulong z1, out ulong z0) {
            unchecked {
                z0 = x0 - y0;
            }
            ulong b = ((~x0 & y0) | ((~(x0 ^ y0)) & z0)) >> 63;
            z1 = x1 - y1 - b;
        }

        public static void Increment128 (ulong x1, ulong x0, out ulong y1, out ulong y0) {
            y0 = x0 + 1UL;
            // This is a specialization of the trick above to avoid branching on (y0 == 0UL).
            ulong c = (x0 & ~y0) >> 63;
            y1 = x1 + c;
        }

        public static void TwosComplement (ulong x1, ulong x0, out ulong y1, out ulong y0) {
            Increment128(~x1, ~x0, out y1, out y0);
        }

        public static void Decompose (ulong u, out ulong u1, out ulong u0) {
            u1 = u >> 32;
            u0 = u & 0xffffffff;
            Debug.Assert(u1 <= uint.MaxValue);
            Debug.Assert(u0 <= uint.MaxValue);
        }

        // Multiplies two 64 bit numbers, and returns the resulting 128-bit number.
        // This is used in the multiplication of two 128-bit numbers.

        public static void Multiply64By64 (ulong u, ulong v, out ulong w1, out ulong w0) {

            Decompose(u, out ulong u1, out ulong u0);
            Decompose(v, out ulong v1, out ulong v0);

            // At fist I tried to do this by just multiplying out
            // (u1 2^32 + u0) * (v1 2^32 + v0) =
            //   (u1 v1) 2^64 + [(u0 v1) + (u1 v0)] 2^32 + (u0 v0) =
            //   [ (u1 v1) + (u0 v1)_hi + (u1 v0)_hi ] 2^64 +
            //   [ (u0 v0) + (u0 v1)_lo 2^32 + (u1 v0)_lo 2^32 ]
            // and then identifying the upper 64 bits with the first
            // bracked expression, and the lower 64 bits with the second.
            // This is _almost_ right, but the expression in the second brackets
            // can overflow, so you have to compute it to get the first.

            // That does work, but it's faster to turn the crank on 2 X 2
            // long multiplication. This gives the upper 64 bits in
            // 4 multiplies, 4 adds, and 4 bit operations. 

            ulong x = u0 * v0;
            ulong y = u1 * v0 + (x >> 32);
            ulong z = u0 * v1 + (y & 0xffffffff);
            w1 = u1 * v1 + (y >> 32) + (z >> 32);

            // The lower 64 bits can of course be constructed similiarly,
            // but that's not efficient, because they are also given
            // by just multiplying the numbers and discarding overflow.
            unchecked {
                w0 = u * v;
            }

        }

        // Multiplies two 128-bit numbers.
        // Overflow of product beyond 128 bits is ignored.

        public static void Multiply128By128 (ulong u1, ulong u0, ulong v1, ulong v0, out ulong p1, out ulong p0) {

            // [u1 2^64 + u0][v1 2^64 + v0] = u1 v1 2^128 + u1 v0 2^64 + u0 v1 2^64 + u0 v0 =
            // u1 v1 + [(u1 v0)_hi 2^64 + (u1 v0)_lo] 2^64 + [(u0 v1)_hi 2^64 + (u0 v1)_lo] 2^64 + (u0 v0)_hi 2^64 + (u0 v0)_lo =
            // [u1 v1 + (u1 v0)_hi + (u0 v1)_hi] 2^128 + [(u1 v0)_lo + (u0 v1)_lo + (u0 v0)_hi ] 2^64 + (u0 v0)_lo

            Multiply64By64(u0, v0, out ulong hi, out ulong low);

            p0 = low;
            p1 = hi + (u1 * v0) + (u0 * v1);

        }

        public static void Divide128By128 (ulong u1, ulong u0, ulong v1, ulong v0, out ulong q1, out ulong q0, out ulong r1, out ulong r0) {

            if (v1 == 0UL) {
                // Denominator fits in 64 bits.
                // So remainder will also fit in 64 bits.
                r1 = 0UL;
                if (u1 == 0UL) {
                    // Numerator fits in 64 bits.
                    // So quotient can be found by native division
                    q1 = 0UL;
                    if (u0 < v0) {
                        q0 = 0UL;
                        r0 = u0;
                    } else {
                        q0 = DivRem(u0, v0, out r0);
                    }
                } else {
                    // Numerator is 65-128 bits
                    // So quotient can be over 64 bits
                    // Since denominator is one 64-bit digit, this is short division.
                    // Native division to get first digit.
                    q1 = DivRem(u1, v0, out ulong k);
                    // At this point k, must be less than v0 by properties of remainder,
                    // and that is required to call into next routine.
                    Debug.Assert(k < v0);
                    // Get second digit. Since we may have a non-zero remainder,
                    // this cannot be native 64-bit by 64-bit division.
                    q0 = Divide128By64(k, u0, v0, out r0);
                }
            } else {
                // Denominator is 65-128 bits
                // So quotient will fit in 64 bits
                q1 = 0UL;

                // Test for early return
                // This is actually necessary to ensure q0 != 0 below.
                if ((u1 < v1) || ((u1 == v1) && (u0 <= v0))) {
                    q0 = 0UL;
                    r1 = u1;
                    r0 = u0;
                    return;
                }

                // Hacker's Delight Section 9-5, pp. 197-2-1, goes into great detail about
                // why the following algorithm works.

                // Normalize the denominator by left-shifting until its most significant bit is 1.
                // And right-shift u by 1.
                // Now v1 has highest bit set, and u1 does not, so we are sure v1 > u1
                int s = NumberOfLeadingZeros(v1);
                Debug.Assert(0 <= s && s <= 63);
                ulong qt = Divide128ByNormalized64(u1 >> 1, (u1 << 63) | (u0 >> 1), (s == 0) ? v1 : (v1 << s) | (v0 >> (64 - s)), out _);
                q0 = qt >> (63 - s);
                Debug.Assert(q0 != 0UL);
                // It doesn't appear the returned remainder has anything to do with our remainder.

                // At this point, q0 is almost always correct. But ~1x in 2^64, it
                // will be too big by 1.
                // To correct for that case, we need to compute product and
                // reduce q0 if it is too big.
                // But if q0 is too big, product can overflow. So instead always
                // reduce q0, compute product, then increase if remainder is
                // still bigger than q.
                // An example of the necessity of this part, given in Hacker's Delight
                // text, is u = 2^128 - 1, v = 2^64 + 3.

                //r1 = 0UL; r0 = 0UL;
                q0--;
                Multiply128By128(0UL, q0, v1, v0, out r1, out r0);
                Subtract128From128(u1, u0, r1, r0, out r1, out r0);
                if ((r1 > v1) || (r1 == v1 && r0 >= v0)) {
                    q0++;
                    Subtract128From128(r1, r0, v1, v0, out r1, out r0);
                }
            }

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
            // We then apply Knuth's long division algorithm, specialized to base 2^32, a 4-digit
            // numerator and 2-digit denominator (yielding a 2-digit quotient and 2-digit remainder).
            // This algorithm requires an arithmetic register twice as wide as the digits,
            // in this case 64-bit, which C# has. Knuth describes the algorithm in Section 4.3.1
            // and this specialization is based on Hacker's Delight section 9-4.

            // To improve guessed quotients, which could otherwise be off by many digits
            // in a large base, the algorithm first "normalizes" the denominator
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
        // in the binary representation of a ulong. We effectively
        // use a binary search to minimize tests. Value is
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

        public static void Divide128By32 (ulong u1, ulong u0, uint v, out ulong q1, out ulong q0, out uint r) {

            Decompose(u1, out ulong u11, out ulong u10);
            Decompose(u0, out ulong u01, out ulong u00);

            ulong q11 = DivRem(u11, v, out ulong s);
            ulong q10 = DivRem((s << 32) | u10, v, out s);
            ulong q01 = DivRem((s << 32) | u01, v, out s);
            ulong q00 = DivRem((s << 32) | u00, v, out s);

            Debug.Assert(s < v);
            r = (uint) s;

            q1 = (q11 << 32) | q10;
            q0 = (q01 << 32) | q00;
        }

        // Replace this by Math.DivRem if it is implemented via a single CPU instruction.

        private static ulong DivRem (ulong u, ulong v, out ulong r) {
            ulong q = u / v;
            r = u - q * v;
            Debug.Assert(r < v);
            return q;
        }

    }
}
