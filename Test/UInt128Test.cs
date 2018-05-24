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

        private static IEnumerable<UInt64> GetRandomUInt64 (Random rng, int count) {

            double scale = Math.Log(UInt64.MaxValue);

            for (int i = 0; i < count; i++) {

                double logValue = rng.NextDouble() * scale;
                UInt64 value = (UInt64) Math.Round(Math.Exp(logValue));
                yield return value;
            }

        }

        private static IEnumerable<UInt128> GetRandomUInt128 (Random rng, int count) {

            double scale = Math.Log(UInt64.MaxValue);

            for (int i = 0; i < count; i++) {
                double logValue1 = rng.NextDouble() * scale;
                double logValue0 = rng.NextDouble() * scale;
                UInt128 value1 = (UInt128) Math.Round(Math.Exp(logValue1));
                UInt128 value0 = (UInt128) Math.Round(Math.Exp(logValue0));
                UInt128 value = value1 * UInt64.MaxValue + value0;
                yield return value;
            }
        }

        [TestMethod]
        public void UInt128Equality () {

            UInt128B zero = UInt128B.Zero;
            UInt128B one = UInt128B.One;
            UInt128B m = UInt64.MaxValue;
            UInt128B n = m + 1;

            Assert.IsTrue(UInt128B.MinValue == UInt128B.Zero);
            Assert.IsTrue(UInt128.Equals(UInt128.MaxValue, UInt128.MaxValue));
            Assert.IsTrue(UInt128B.One.Equals((object) UInt128B.One));

            Assert.IsTrue(zero != one);
            Assert.IsTrue(!one.Equals(zero));

            Assert.IsTrue(m + one == n);
            Assert.IsTrue(n.Equals(m + one));
            Assert.IsTrue(!m.Equals(null));
            Assert.IsTrue((m + one).GetHashCode() == n.GetHashCode());
        }

        [TestMethod]
        public void UInt128Comparison () {

            UInt128B zero = UInt128B.Zero;
            UInt128B one = UInt128B.One;
            UInt128B n = ((UInt128B) UInt64.MaxValue) + 1;
            UInt128B p = UInt128B.MaxValue - 1;

            Assert.IsTrue(zero < one);
            Assert.IsTrue(zero < n);

            Assert.IsTrue(one > zero);
            Assert.IsTrue(one < n);

            Assert.IsTrue(n > one);
            Assert.IsTrue(n < p);

            Assert.IsTrue(p >= n);

            Assert.IsTrue(n.CompareTo(n) == 0);
            Assert.IsTrue(n.CompareTo(p) == -1);

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
        public void UInt128Overflow () {

            UInt128B two = 2;
            unchecked {
                Assert.IsTrue(UInt128B.MaxValue + UInt128B.One == UInt128B.Zero); // 99 + 1 = [1]00
                Assert.IsTrue(UInt128B.MaxValue + UInt128B.MaxValue == UInt128B.MaxValue - UInt128B.One); // 99 + 99 = [1]98
                Assert.IsTrue(UInt128B.Zero - UInt128B.One == UInt128B.MaxValue); // [1]00 - 1 = 99
                Assert.IsTrue(UInt128B.Zero - UInt128B.MaxValue == UInt128B.One); // [1]00 - 99 = 1
                Assert.IsTrue(UInt128B.MaxValue * two == UInt128B.MaxValue - UInt128B.One); // 99 * 2 = [1]98
                Assert.IsTrue(UInt128B.MaxValue * UInt128B.MaxValue == UInt128B.One); // 99 * 99 = [98]01
            }

        }

        [TestMethod]
        public void UInt128Parsing () {
            UInt128 q = UInt128.Parse("18446744073709551616");
            Assert.IsTrue((UInt128) UInt64.MaxValue + (UInt128) 1UL == q);
        }

        [TestMethod]
        public void UInt128Printing () {
            int u = 4032010;
            UInt128 v = (ulong) u;
            string text = v.ToString();
            Assert.IsTrue(text == u.ToString());

            string maxText = UInt128.MaxValue.ToString();
        }

        [TestMethod]
        public void UInt128UInt64Agreement () {
            Random rng = new Random(1);
            foreach (UInt64 a in GetRandomUInt64(rng, 4)) {
                UInt128 aa = a;
                Assert.IsTrue(aa.ToString() == a.ToString());
                foreach (UInt64 b in GetRandomUInt64(rng, 4)) {
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
        public void CompareToBigInteger () {

            Random rng = new Random(1);
            UInt128[] values = GetRandomUInt128(rng, 1000).ToArray();

            Stopwatch t1 = Stopwatch.StartNew();
            for (int i = 0; i < values.Length; i++) {
                for (int j = 0; j < values.Length; j++) {
                    UInt128 q = UInt128.DivRem(values[i], values[j], out _);
                    //UInt128 sum = values[i] - values[j];
                }
            }
            t1.Stop();
            Console.WriteLine(t1.ElapsedMilliseconds);

            BigInteger[] bigs = new BigInteger[values.Length];
            for (int i = 0; i < bigs.Length; i++) {
                bigs[i] = BigInteger.Parse(values[i].ToString());
            }

            Stopwatch t2 = Stopwatch.StartNew();
            for (int i = 0; i < values.Length; i++) {
                for (int j = 0; j < values.Length; j++) {
                    BigInteger r;
                    BigInteger q = BigInteger.DivRem(bigs[i], bigs[j], out r);
                    //BigInteger q = BigInteger.Add(bigs[i], bigs[j]);
                    //BigInteger sum = bigs[i] - bigs[j];
                }
            }
            t2.Stop();
            Console.WriteLine(t2.ElapsedMilliseconds);

            UInt128B[] bValues = new UInt128B[values.Length];
            for(int i = 0; i < bValues.Length; i++) {
                uint s0 = (uint) rng.Next();
                uint s1 = (uint) rng.Next();
                uint s2 = (uint) rng.Next();
                uint s3 = (uint) rng.Next();
                bValues[i] = new UInt128B(s3, s2, s1, s0);
            }

            Stopwatch t3 = Stopwatch.StartNew();
            for (int i = 0; i < values.Length; i++) {
                for (int j = 0; j < values.Length; j++) {
                    //UInt128B p = bValues[i] - bValues[j];
                    UInt128B q = bValues[i] / bValues[j];
                }
            }
            t3.Stop();
            Console.WriteLine(t3.ElapsedMilliseconds);

        }

        [TestMethod]
        public void BTest () {

            UInt128B m2 = UInt128B.Multiply(UInt128B.MaxValue, UInt128B.MaxValue);
            UInt128B n2 = UInt128B.Add(UInt128B.MaxValue, UInt128B.MaxValue);

            uint r;
            UInt128B.DivRem(UInt64.MaxValue, 12345, out r);

            UInt128B n = UInt128B.Parse("123456789012345678901234567891234567890");
            UInt128B d = UInt128B.Parse("9876543210987654321");

            UInt128B q = n / d;

        }

    }
}
