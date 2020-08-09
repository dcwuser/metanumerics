using System;
using System.Collections.Generic;

using System.Diagnostics;
using System.Linq;
using System.Numerics;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Extended;

namespace Test {

    [TestClass]
    public class UInt128Test {

        private static IEnumerable<ulong> GetRandomUInt64 (int count, Random rng = null) {
            if (rng == null) rng = new Random(1);
            for (int i = 0; i < count; i++) {
                ulong value = (ulong) Math.Round(Math.Exp(rng.NextDouble() * logMax64));
                //ulong value = (ulong) Math.Round(UInt64.MaxValue * rng.NextDouble());
                yield return value;
            }
        }

        private static double logMax64 = Math.Log(UInt64.MaxValue);

        private static IEnumerable<uint> GetRandomUInt32 (Random rng, int count) {
            for (int i = 0; i < count; i++) {
                uint value = (uint) Math.Round(Math.Exp(rng.NextDouble() * logMax32));
                //uint value = (uint) rng.Next();
                yield return value;
            }
        }

        private static double logMax32 = Math.Log(UInt32.MaxValue);

        private static IEnumerable<UInt128> GetRandomUInt128 (Random rng, int count) {
            for (int i = 0; i < count; i++) {
                UInt128 u = (UInt128) Math.Round(Math.Exp(rng.NextDouble() * logMax64));
                //UInt128 u = (UInt128) Math.Round(rng.NextDouble() * UInt64.MaxValue);
                if (rng.NextDouble() > 0.5) {
                    UInt128 u0 = (UInt128) Math.Round(rng.NextDouble() * UInt64.MaxValue);
                    u = u * UInt64.MaxValue + u0;
                }
                yield return u;
            }
        }

        [TestMethod]
        public void UInt128Equality () {

            UInt128 zero = UInt128.Zero;
            UInt128 one = UInt128.One;
            UInt128 m = UInt64.MaxValue;
            UInt128 n = m + 1;

            Assert.IsTrue(UInt128.MinValue == UInt128.Zero);
            Assert.IsTrue(UInt128.Equals(UInt128.MaxValue, UInt128.MaxValue));
            Assert.IsTrue(UInt128.One.Equals((object) UInt128.One));

            Assert.IsTrue(zero != one);
            Assert.IsTrue(!one.Equals(zero));

            Assert.IsTrue(m + one == n);
            Assert.IsTrue(n.Equals(m + one));
            Assert.IsTrue(!m.Equals(null));
            Assert.IsTrue((m + one).GetHashCode() == n.GetHashCode());
        }

        [TestMethod]
        public void UInt128Comparison () {
            UInt128[] orderedValues = new UInt128[] { UInt128.Zero, UInt128.One, Int64.MaxValue /* why without cast */, UInt64.MaxValue, /*Int128.MaxValue, */ UInt128.MaxValue };
            for (int i = 0; i < orderedValues.Length; i++) {
                for (int j = 0; j < orderedValues.Length; j++) {
                    if (i < j) {
                        Assert.IsTrue(UInt128.Compare(orderedValues[i], orderedValues[j]) == -1);
                        Assert.IsTrue(orderedValues[i].CompareTo(orderedValues[j]) == -1);
                        Assert.IsTrue(orderedValues[i] < orderedValues[j]);
                        Assert.IsFalse(orderedValues[i] >= orderedValues[j]);
                    } else if (i == j) {
                        Assert.IsTrue(UInt128.Compare(orderedValues[i], orderedValues[j]) == 0);
                        Assert.IsTrue(orderedValues[i].CompareTo(orderedValues[j]) == 0);
                        Assert.IsFalse(orderedValues[i] < orderedValues[j]);
                        Assert.IsFalse(orderedValues[i] > orderedValues[j]);
                        Assert.IsTrue(orderedValues[i] <= orderedValues[j]);
                    } else {
                        Assert.IsTrue(UInt128.Compare(orderedValues[i], orderedValues[j]) == +1);
                        Assert.IsTrue(orderedValues[i].CompareTo(orderedValues[j]) == +1);
                        Assert.IsTrue(orderedValues[i] > orderedValues[j]);
                        Assert.IsFalse(orderedValues[i] <= orderedValues[j]);
                    }
                }
            }

        }

        [TestMethod]
        public void UInt128Arithmetic () {

            UInt128 zero = UInt128.Zero;
            UInt128 one = UInt128.One;
            UInt128 m = UInt64.MaxValue;
            UInt128 n = m + 1;
            UInt128 p = UInt128.MaxValue - 1;
            UInt128 q = UInt128.MaxValue;

            // Comutativity
            Assert.IsTrue(m + n == n + m);

            // Associativity
            Assert.IsTrue((one + m) + n == one + (m + n));

            // Adding zero
            Assert.IsTrue(one + zero == one);
            Assert.IsTrue(zero + m == m);
            Assert.IsTrue(n == zero + n);
            Assert.IsTrue(q == q + zero);

            // Subtraction of self
            Assert.IsTrue(zero - zero == zero);
            Assert.IsTrue(one - one == zero);
            Assert.IsTrue(n - n == zero);
            Assert.IsTrue(q - q == zero);

            // Multiplication by one
            Assert.IsTrue(zero * one == zero);
            Assert.IsTrue(one * one == one);
            Assert.IsTrue(one * n == n);
            Assert.IsTrue(q == q * one);

            // Division by one
            Assert.IsTrue(zero / one == zero);
            Assert.IsTrue(one / one == one);
            Assert.IsTrue(n / one == n);
            Assert.IsTrue(p / one == p);

            // Basic arithmetic
            UInt128 two = one + one;
            UInt128 four = one + one + one + one;
            Assert.IsTrue(two * two == four);

        }

        [TestMethod]
        public void UInt128Increment () {
            foreach (UInt128 u in GetRandomUInt128(new Random(1), 4)) {
                UInt128 v = u;
                UInt128 p = v++;
                Assert.IsTrue(v == u + 1);
                Assert.IsTrue(p == u);

                v = u;
                p = v--;
                Assert.IsTrue(v == u - 1);
                Assert.IsTrue(p == u);

                v = u;
                p = ++v;
                Assert.IsTrue(v == u + 1);
                Assert.IsTrue(p == v);
            }
        }


        [TestMethod]
        public void UInt128Overflow () {
            unchecked {
                Assert.IsTrue(UInt128.MaxValue + UInt128.One == UInt128.Zero); // 99 + 1 = [1]00
                Assert.IsTrue(UInt128.MaxValue + UInt128.MaxValue == UInt128.MaxValue - UInt128.One); // 99 + 99 = [1]98
                Assert.IsTrue(UInt128.Zero - UInt128.One == UInt128.MaxValue); // [1]00 - 1 = 99
                Assert.IsTrue(UInt128.Zero - UInt128.MaxValue == UInt128.One); // [1]00 - 99 = 1
                Assert.IsTrue(UInt128.MaxValue * 2 == UInt128.MaxValue - UInt128.One); // 99 * 2 = [1]98
                Assert.IsTrue(UInt128.MaxValue * UInt128.MaxValue == UInt128.One); // 99 * 99 = [98]01

                UInt128 m = UInt128.MaxValue;
                m++;
                Assert.IsTrue(m == UInt128.Zero);
                m--;
                Assert.IsTrue(m == UInt128.MaxValue);
            }
        }

        [TestMethod]
        public void UInt128SpecialStrings () {
            Assert.IsTrue(UInt128.Zero.ToString() == "0");
            Assert.IsTrue(UInt128.One.ToString() == "1");

            Assert.IsTrue(new UInt128("0") == UInt128.Zero);
            Assert.IsTrue(new UInt128("1") == UInt128.One);
        }

        [TestMethod]
        public void UInt128ParseRoundtrip () {
            Random rng = new Random(1);
            foreach (UInt128 u in GetRandomUInt128(rng, 4)) {
                Assert.IsTrue(UInt128.Parse(u.ToString()) == u);
            }
        }

        [TestMethod]
        public void UInt128DoubleRoundtrip () {

            // For integers exactly representable by doubles (52 or fewer binary digits),
            // roundtrip through double should preserve value
            Random rng = new Random(8);
            foreach (UInt128 u in GetRandomUInt64(8, rng)) {
                if ((u >> 52) != UInt128.Zero) continue;
                double d = (double) u;
                UInt128 u1 = (UInt128) d;
                Assert.IsTrue(u1 == u);
            }

        }

        [TestMethod]
        public void UIn128CastFailure () {
            Assert.ThrowsException<InvalidCastException>(() => (UInt128) Double.NaN);
        }

        [TestMethod]
        public void UInt128ParseFailure () {

            Assert.IsFalse(UInt128.TryParse(null, out _));
            Assert.ThrowsException<ArgumentNullException>(() => UInt128.Parse(null));

            Assert.IsFalse(UInt128.TryParse(String.Empty, out _));
            Assert.ThrowsException<FormatException>(() => UInt128.Parse(String.Empty));

            Assert.IsFalse(UInt128.TryParse("-1", out _));
            Assert.ThrowsException<FormatException>(() => UInt128.Parse("-1"));

            Assert.IsFalse(UInt128.TryParse("+", out _));
            Assert.ThrowsException<FormatException>(() => UInt128.Parse("+"));

            // Max is parsable, but one more isn't
            Assert.IsTrue(UInt128.TryParse(UInt128.MaxValue.ToString(), out UInt128 max));
            Assert.IsTrue(max == UInt128.MaxValue);
            string tooBigByOne = MakeNumberOneBigger(UInt128.MaxValue.ToString());
            Assert.IsFalse(UInt128.TryParse(tooBigByOne, out _));
            Assert.ThrowsException<OverflowException>(() => UInt128.Parse(tooBigByOne));

        }

        // This is very fast-and-loose method increase the magnitude of a number represented as text by one.
        // It doesn't handle trim and will fail if the number is all 9s, but it's good enough for our test purposes.
        public static string MakeNumberOneBigger (string number) {
            char[] characters = number.ToCharArray();
            for (int i = characters.Length - 1; i > 0; i--) {
                if (characters[i] == '9') {
                    characters[i] = '0';
                } else {
                    characters[i] = (char) ((int) characters[i] + 1);
                    return new String(characters);
                }
            }
            throw new InvalidOperationException();
        }

        [TestMethod]
        public void UInt128UInt64Agreement () {
            Random rng = new Random(64);
            foreach (UInt64 a in GetRandomUInt64(8, rng)) {
                UInt128 aa = a;
                Assert.IsTrue(aa.ToString() == a.ToString());
                foreach (UInt64 b in GetRandomUInt64(8, rng)) {
                    UInt128 bb = b;

                    if (aa + bb < UInt64.MaxValue) Assert.IsTrue(a + b == aa + bb);

                    if (a >= b) {
                        Assert.IsTrue(a - b == aa - bb);
                    } else {
                        Assert.IsTrue(b - a == bb - aa);
                    }

                    if (aa * bb < UInt64.MaxValue) Assert.IsTrue(a * b == aa * bb);

                    Assert.IsTrue(a / b == aa / bb);
                    Assert.IsTrue(b / a == bb / aa);
                    Assert.IsTrue(a % b == aa % bb);
                    Assert.IsTrue(b % a == bb % aa);
                }
            }
        }

        [TestMethod]
        public void UInt128BigIntegerAgreement () {

            BigInteger b = BigInteger.Parse(UInt128.MaxValue.ToString()) + 1;

            foreach (UInt128 u in GetRandomUInt128(new Random(2), 16)) {
                BigInteger u1 = (BigInteger) u;
                foreach (UInt128 v in GetRandomUInt128(new Random(3), 16)) {
                    BigInteger v1 = (BigInteger) v;

                    // Comparison
                    Assert.IsTrue(UInt128.Compare(u, v) == BigInteger.Compare(u1, v1));
                    Assert.IsTrue((u < v) == (u1 < v1));
                    Assert.IsTrue((u > v) == (u1 > v1));

                    // Sum
                    UInt128 s = u + v;
                    BigInteger s1 = (u1 + v1) % b;
                    Assert.IsTrue((BigInteger) s == s1);

                    // Difference
                    if (u >= v) {
                        UInt128 d = u - v;
                        BigInteger d1 = u1 - v1;
                        Assert.IsTrue((BigInteger) d == d1);
                    }

                    // Product
                    UInt128 p = u * v;
                    BigInteger p1 = (u1 * v1) % b;
                    Assert.IsTrue((BigInteger) p == p1);

                    // Quotient
                    UInt128 q = u / v;
                    BigInteger q1 = u1 / v1;
                    Assert.IsTrue((BigInteger) q == q1);

                }
            }

        }

        [TestMethod]
        public void UInt128RandomArithmetic () {
            Random rng = new Random(1);
            foreach (UInt128 a in GetRandomUInt128(rng, 8)) {

                // ToString and Parse round-trip
                UInt128 rt = UInt128.Parse(a.ToString());
                Assert.IsTrue(rt == a);

                foreach (UInt128 b in GetRandomUInt128(rng, 8)) {

                    // Subtraction and addition are inverses
                    UInt128 d = a - b;
                    Assert.IsTrue(d + b == a);

                    UInt128 c = b - a;
                    Assert.IsTrue(a + c == b);

                    // Division and multiplication are inverses
                    UInt128 r;
                    UInt128 q = UInt128.DivRem(a, b, out r);
                    Assert.IsTrue(q * b + r == a);
                    Assert.IsTrue(a / b == q);
                    Assert.IsTrue(a % b == r);

                    UInt128 s = b / a;
                    UInt128 t = b % a;
                    Assert.IsTrue(s * a + t == b);
                }
            }
        }

        [TestMethod]
        public void PolynomialFactors () {

            // t^2 - 1 = (t+1)(t-1)
            UInt128 m = (UInt128) UInt64.MaxValue;
            UInt128 q = m * m - 1;
            Assert.IsTrue(q / (m - 1) == m + 1);
            Assert.IsTrue(q / (m + 1) == m - 1);

            Random rng = new Random(10);
            foreach (UInt128 t in GetRandomUInt64(8)) {
                UInt128 p = t * t - 1;
                Assert.IsTrue(p / (t + 1) == t - 1);
                Assert.IsTrue(p / (t - 1) == t + 1);
            }

        }


        [TestMethod]
        public void CornerCaseDivision () {

            // This division problem is chosen to exercise the rare case
            // when our initial estimated quotient is off-by-one.

            UInt128 u = UInt128.MaxValue;
            UInt128 v = (UInt128.One << 64) + 3;
            //UInt128 v = new UInt128(1UL, 3UL);

            UInt128 q = UInt128.DivRem(u, v, out UInt128 r);
            Assert.IsTrue(q * v + r == u);

        }

        [TestMethod]
        public void DivisionByZero () {

            Assert.ThrowsException<DivideByZeroException>(() => UInt128.MaxValue / UInt128.Zero);
            Assert.ThrowsException<DivideByZeroException>(() => UInt64.MaxValue / UInt128.Zero);
            Assert.ThrowsException<DivideByZeroException>(() => UInt128.Zero / UInt128.Zero);
            Assert.ThrowsException<DivideByZeroException>(() => UInt128.One / 0UL);
            Assert.ThrowsException<DivideByZeroException>(() => UInt128.One / 0U);

        }

        [TestMethod]
        public void CompareDivision () {
            Random rng = new Random(222);

            UInt128[] us = GetRandomUInt128(rng, 1000).ToArray();
            uint[] vs = new uint[1000];
            for (int i = 0; i < vs.Length; i++) vs[i] = (uint) rng.Next();

            Stopwatch t1 = Stopwatch.StartNew();
            foreach (UInt128 u in us) {
                foreach (uint v in vs) {
                    UInt128 q1 = UInt128.DivRem(u, (UInt128) v, out _);
                }
            }
            t1.Stop();
            Console.WriteLine(t1.ElapsedMilliseconds);

            Stopwatch t2 = Stopwatch.StartNew();
            foreach (UInt128 u in us) {
                foreach (uint v in vs) {
                    UInt128 q1 = UInt128.DivRem(u, v, out _);
                }
            }
            t2.Stop();
            Console.WriteLine(t2.ElapsedMilliseconds);



        }

        [TestMethod]
        public void UInt128TimingComparison () {

            UInt128[] values = GetRandomUInt128(new Random(1), 1000).ToArray();
            TimeOperations<UInt128>("Addition", (x, y) => { UInt128 z = x + y; }, values);
            TimeOperations<UInt128>("Subtraction", (x, y) => { UInt128 z = x - y; }, values);
            TimeOperations<UInt128>("Multiplication", (x, y) => { UInt128 z = x * y; }, values);
            TimeOperations<UInt128>("Division", (x, y) => { UInt128.DivRem(x, y, out _); }, values);

            BigInteger[] bigValues = values.Select(x => (BigInteger) x).ToArray();
            TimeOperations<BigInteger>("Big Addition", (x, y) => { BigInteger z = x + y; }, bigValues);
            TimeOperations<BigInteger>("Big Subtraction", (x, y) => { BigInteger z = x - y; }, bigValues);
            TimeOperations<BigInteger>("Big Multiplication", (x, y) => { BigInteger z = x * y; }, bigValues);
            TimeOperations<BigInteger>("Big Division", (x, y) => { BigInteger.DivRem(x, y, out _); }, bigValues);

            uint[] smallValues = GetRandomUInt32(new Random(2), 1000).ToArray();
            BigInteger[] smallAsBigValues = smallValues.Select(x => (BigInteger) x).ToArray();
            //TimeOperations2<UInt128, UInt128>("Division 2L", (x, y) => { UInt128 z = UInt128.DivRem(x, y, out _); }, values, smallValuesAsValues);
            TimeOperations2<UInt128, uint>("Division 128/32", (x, y) => { UInt128 z = UInt128.DivRem(x, y, out _); }, values, smallValues);
            TimeOperations2<BigInteger, uint>("Big Division 128/32", (x, y) => { BigInteger.DivRem(x, y, out _); }, bigValues, smallValues);
        }

        private void TimeOperations<T> (string name, Action<T, T> operation, IEnumerable<T> values) {
            Stopwatch t = Stopwatch.StartNew();
            foreach (T x in values) {
                foreach (T y in values) {
                    operation(x, y);
                }
            }
            t.Stop();
            Console.WriteLine($"{name}: {t.ElapsedMilliseconds}");
        }

        private void TimeOperations2<T1, T2> (string name, Action<T1, T2> operation, IEnumerable<T1> values1, IEnumerable<T2> values2) {
            Stopwatch t = Stopwatch.StartNew();
            foreach (T1 x in values1) {
                foreach (T2 y in values2) {
                    operation(x, y);
                }
            }
            t.Stop();
            Console.WriteLine($"{name}: {t.ElapsedMilliseconds}");
        }

        [TestMethod]
        public void CompareToBigInteger () {

            Random rng = new Random(1);

            //UInt64[] values64 = GetRandomUInt64(rng, 1000).ToArray();
            //UInt128[] values = new UInt128[values64.Length];
            //for (int i = 0; i < values.Length; i++) values[i] = values64[i];

            UInt128[] values = GetRandomUInt128(rng, 1000).ToArray();

            Stopwatch t1 = Stopwatch.StartNew();
            for (int i = 0; i < values.Length; i++) {
                for (int j = 0; j < values.Length; j++) {
                    //UInt128 sum = values[i] + values[j];
                    //UInt128 difference = values[i] - values[j];
                    //UInt128 product = values[i] * values[j];
                    UInt128 quotient = values[i] % values[j];
                }
            }
            t1.Stop();
            Console.WriteLine(t1.ElapsedMilliseconds);

            // BigInteger
            BigInteger[] bigs = new BigInteger[values.Length];
            for (int i = 0; i < bigs.Length; i++) {
                bigs[i] = BigInteger.Parse(values[i].ToString());
            }
            Stopwatch t2 = Stopwatch.StartNew();
            for (int i = 0; i < values.Length; i++) {
                for (int j = 0; j < values.Length; j++) {
                    //BigInteger sum = bigs[i] + bigs[j];
                    //BigInteger difference = bigs[i] - bigs[j];
                    //BigInteger product = bigs[i] * bigs[j];
                    BigInteger quotient = bigs[i] / bigs[j];
                }
            }
            t2.Stop();
            Console.WriteLine(t2.ElapsedMilliseconds);

            /*
            // UInt128B
            UInt128B[] bValues = new UInt128B[values.Length];
            for(int i = 0; i < bValues.Length; i++) {
                bValues[i] = UInt128B.Parse(values[i].ToString());
            }
            Stopwatch t3 = Stopwatch.StartNew();
            for (int i = 0; i < values.Length; i++) {
                for (int j = 0; j < values.Length; j++) {
                    //UInt128B sum = bValues[i] + bValues[j];
                    //UInt128B difference = bValues[i] - bValues[j];
                    //UInt128B product = bValues[i] * bValues[j];
                    UInt128B quotient = bValues[i] / bValues[j];
                }
            }
            t3.Stop();
            Console.WriteLine(t3.ElapsedMilliseconds);
            */

            Stopwatch t4 = Stopwatch.StartNew();
            for (int i = 0; i < values.Length; i++) {
                for (int j = 0; j < values.Length; j++) {
                    //UInt128 sum = values[i] + values[j];
                    //UInt128 difference = values[i] - values[j];
                    //UInt128 product = values[i] * values[j];
                    UInt128 quotient = values[i] / values[j];
                }
            }
            t4.Stop();
            Console.WriteLine(t4.ElapsedMilliseconds);


        }

        [TestMethod]
        public void DeMorgansLaws () {

            Random rng = new Random(6);
            foreach (UInt128 x in GetRandomUInt128(rng, 4)) {
                foreach (UInt128 y in GetRandomUInt128(rng, 4)) {
                    Assert.IsTrue(~(x & y) == (~x | ~y));
                    Assert.IsTrue(~(x | y) == (~x & ~y));
                    //Assert.IsTrue(~(x ^ y) == Convert.ToUInt64(x == y));

                    Assert.IsTrue(~(x + y) == ~x - y);
                    Assert.IsTrue(~(x - y) == ~x + y);

                    Assert.IsTrue((x ^ y) <= (x | y));
                }
                Assert.IsTrue(~(x + 1) == ~x - 1);
            }

        }

        [TestMethod]
        public void BitshiftAgreement () {

            UInt128 mask = ((UInt128) UInt64.MaxValue) << 64;
            Random rng = new Random(1);
            foreach (UInt64 u in GetRandomUInt64(8, rng)) {
                UInt128 v = (UInt128) u;
                for (int i = 0; i < 4; i++) {
                    int n = rng.Next(0, 63);
                    Assert.IsTrue((v >> n) == (u >> n));
                    Assert.IsTrue(((UInt64) (v << n)) == (u << n));
                    Assert.IsTrue(((v << 64) >> (64 + n)) == (u >> n)); 
                }

            }

        }


        // Binary long division
        // This tests division and bit operations

        [TestMethod]
        public void BinaryLongDivisionAgreement () {
            Random rng = new Random(12);
            foreach (UInt128 u in GetRandomUInt128(rng, 8)) {
                foreach (UInt128 v in GetRandomUInt128(rng, 8)) {

                    UInt128 q = BinaryLongDivision(u, v, out UInt128 r);

                    UInt128 q1 = UInt128.DivRem(u, v, out UInt128 r1);

                    Assert.IsTrue(q == q1);
                    Assert.IsTrue(r == r1);

                }
            }
        }

        private static UInt128 BinaryLongDivision (UInt128 u, UInt128 v, out UInt128 r) {
            Debug.Assert(v != UInt128.Zero);

            // Binary long division

            UInt128 q = UInt128.Zero;
            r = u;

            int nu = MostSignificantBit(u);
            int nv = MostSignificantBit(v);

            for (int i = nu - nv; i >= 0; i--) {
                UInt128 w = v << i;
                if (w <= r) {
                    q = SetBit(q, i);
                    r = r - w;
                }
            }

            Debug.Assert(r < v);
            Debug.Assert(q * v + r == u);

            return q;
        }

        private static int MostSignificantBit (UInt128 u) {
            int i = -1;
            while (u != UInt128.Zero) {
                u = u >> 1;
                i++;
            }
            return i;
        }

        private static UInt128 SetBit (UInt128 u, int i) {
            Debug.Assert(0 <= i && i < 128);
            UInt128 mask = UInt128.One << i;
            return u | mask;
        }

    }
}
