using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Extended;

namespace Test {

    [TestClass]
    public class Int128Test {

        private static IEnumerable<long> GetRandomInt64 (int count, Random rng = null) {
            if (rng is null) rng = new Random(1);
            for (int i = 0; i < count; i++) {
                long value = (long) Math.Round(Math.Exp(rng.NextDouble() * logMax64));
                if (rng.NextDouble() < 0.5) value = -value;
                yield return value;
            }
        }

        private static double logMax64 = Math.Log(Int64.MaxValue);

        private static IEnumerable<Int128> GetRandomInt128 (int count, Random rng = null) {
            if (rng is null) rng = new Random(1);
            for (int i = 0; i < count; i++) {
                Int128 u = (Int128) Math.Round(Math.Exp(rng.NextDouble() * logMax64));
                if (rng.NextDouble() > 0.5) {
                    // This isn't quite right. Should be UInt164
                    Int128 u0 = (Int128) Math.Round(rng.NextDouble() * Int64.MaxValue);
                    u = u * Int64.MaxValue + u0;
                }
                if (rng.NextDouble() < 0.5) u = -u;
                yield return u;
            }
        }

        [TestMethod]
        public void Int128Equality () {

            // Reflexivity a = a
            // Symmetry a = b => b = a
            // Transitivity a = b & b = c => a = c

            Int128[] values = new Int128[] { Int128.MinValue, Int64.MinValue, Int32.MinValue, -Int128.One, Int128.Zero, Int128.One, Int64.MaxValue, Int128.MaxValue };

            foreach (Int128 x in values) {
                Assert.IsTrue(x.Equals(x));
                Assert.IsFalse(x.Equals(null));
                foreach (Int128 y in values) {
                    if (x == y) {
                        Assert.IsTrue(y == x);
                        Assert.IsFalse(x != y);
                        Assert.IsTrue(y.Equals(x));
                        Assert.IsTrue(x.Equals((object) y));
                        Assert.IsTrue(x.GetHashCode() == y.GetHashCode());
                    } else {
                        Assert.IsTrue(x != y);
                        Assert.IsFalse(y.Equals(x));
                        Assert.IsFalse(x.Equals((object) y));
                    }
                }
            }

        }

        [TestMethod]
        public void Int128Comparison () {
            Int128[] orderedValues = new Int128[] { Int128.MinValue, Int64.MinValue, Int32.MinValue, -Int128.One, Int128.Zero, Int128.One, Int32.MaxValue, Int64.MaxValue, Int128.MaxValue };
            for (int i = 0; i < orderedValues.Length; i++) {
                for (int j = 0; j < orderedValues.Length; j++) {
                    if (i < j) {
                        Assert.IsTrue(Int128.Compare(orderedValues[i], orderedValues[j]) == -1);
                        Assert.IsTrue(orderedValues[i].CompareTo(orderedValues[j]) == -1);
                        Assert.IsTrue(orderedValues[i] < orderedValues[j]);
                        Assert.IsFalse(orderedValues[i] >= orderedValues[j]);
                    } else if (i == j) {
                        Assert.IsTrue(Int128.Compare(orderedValues[i], orderedValues[j]) == 0);
                        Assert.IsTrue(orderedValues[i].CompareTo(orderedValues[j]) == 0);
                        Assert.IsFalse(orderedValues[i] < orderedValues[j]);
                        Assert.IsFalse(orderedValues[i] > orderedValues[j]);
                        Assert.IsTrue(orderedValues[i] <= orderedValues[j]);
                    } else {
                        Assert.IsTrue(Int128.Compare(orderedValues[i], orderedValues[j]) == +1);
                        Assert.IsTrue(orderedValues[i].CompareTo(orderedValues[j]) == +1);
                        Assert.IsTrue(orderedValues[i] > orderedValues[j]);
                        Assert.IsFalse(orderedValues[i] <= orderedValues[j]);
                    }
                }
            }
        }

        [TestMethod]
        public void Int128Negation () {

            Random rng = new Random(127);
            foreach (Int128 x in GetRandomInt128(8, rng)) {
                foreach (Int128 y in GetRandomInt128(8, rng)) {

                    Assert.IsTrue(-(x + y) == -x - y);

                    Assert.IsTrue((x - y) == -(y - x));

                    Assert.IsTrue(x * (-y) == -(x * y));
                    Assert.IsTrue((-x) * y == -(x * y));
                    Assert.IsTrue((-x) * (-y) == x * y);

                    Assert.IsTrue(x / (-y) == -(x / y));
                    Assert.IsTrue((-x) / y == -(x / y));
                    Assert.IsTrue((-x) / (-y) == x / y);

                }
            }

        }

        [TestMethod]
        public void Int128Arithmetic () {

            Int128[] values = new Int128[] { Int128.MinValue, Int64.MinValue, Int128.Zero, Int128.One, Int64.MaxValue + Int128.One, Int128.MaxValue };

            foreach (Int128 u in values) {

                // Adding zero
                Assert.IsTrue(u + Int128.Zero == u);
                Assert.IsTrue(Int128.Zero + u == u);

                // Subtracting zero
                Assert.IsTrue(u - Int128.Zero == u);

                // Subtracting self
                Assert.IsTrue(u - u == Int128.Zero);

                // Additive inverse
                Assert.IsTrue(u + (-u) == Int128.Zero);
                Assert.IsTrue((-u) + u == Int128.Zero);

                // Multiplying by zero
                Assert.IsTrue(u * Int128.Zero == Int128.Zero);
                Assert.IsTrue(Int128.Zero * u == Int128.Zero);

                // Multiplying by one
                Assert.IsTrue(u * Int128.One == u);
                Assert.IsTrue(Int128.One * u == u);

                foreach (Int128 v in values) {

                    // Commutivity
                    Assert.IsTrue(u + v == v + u);
                    Assert.IsTrue(u * v == v * u);

                    foreach (Int128 w in values) {

                        // Associativity
                        Assert.IsTrue((u + v) + w == u + (v + w));
                        Assert.IsTrue((u - v) - w == u - (v + w));
                        Assert.IsTrue((u * v) * w == u * (v * w));

                        // Distributivity
                        Assert.IsTrue(u * (v + w) == u * v + u * w);
                    }

                }

            }

        }

        [TestMethod]
        public void Int128RandomArithmetic () {
            Random rng = new Random(314159);
            foreach (Int128 x in GetRandomInt128(8, rng)) {

                // Relation of addition to multiplication
                Assert.IsTrue(-x == -1 * x);
                Assert.IsTrue(0 == 0 * x);
                Assert.IsTrue(x == 1 * x);
                Assert.IsTrue(x + x == 2 * x);
                Assert.IsTrue(x + x + x == 3 * x);

                foreach (Int128 y in GetRandomInt128(8, rng)) {

                    // Subtraction and addition are inverses
                    Assert.IsTrue(y + (x - y) == x);

                    // Division and multiplication are inverses
                    Int128 q = Int128.DivRem(x, y, out Int128 r);
                    Assert.IsTrue(q == x / y);
                    Assert.IsTrue(q * y + r == x);

                }
            }
        }

        [TestMethod]
        public void DivisionByZero () {

            Assert.ThrowsException<DivideByZeroException>(() => Int128.MaxValue / Int128.Zero);
            Assert.ThrowsException<DivideByZeroException>(() => Int128.MinValue / Int128.Zero);
            Assert.ThrowsException<DivideByZeroException>(() => Int128.One / 0);
            Assert.ThrowsException<DivideByZeroException>(() => Int128.Zero / Int128.Zero);

        }


        [TestMethod]
        public void Int128Overflow () {

            Assert.IsTrue(Int128.MinValue - 1 == Int128.MaxValue);
            Assert.IsTrue(Int128.MaxValue + 1 == Int128.MinValue);

        }

        [TestMethod]
        public void Int128SepcialStrings () {

            Assert.IsTrue((-Int128.One).ToString() == "-1");
            Assert.IsTrue(new Int128("-1") == -Int128.One);

            Assert.IsTrue(Int128.Zero.ToString() == "0");
            Assert.IsTrue(new Int128("0") == Int128.Zero);

            Assert.IsTrue(Int128.One.ToString() == "1");
            Assert.IsTrue(new Int128("1") == Int128.One);

        }

        [TestMethod]
        public void Int128ParseRoundtrip () {

            // Random values are round-tripped.
            Random rng = new Random(128);
            foreach (Int128 x in GetRandomInt128(8, rng)) {
                Assert.IsTrue(Int128.Parse(x.ToString()) == x);
            }

        }

        [TestMethod]
        public void Int128ParseFailure () {

            Assert.IsFalse(Int128.TryParse(null, out _));
            Assert.ThrowsException<ArgumentNullException>(() => Int128.Parse(null));

            Assert.IsFalse(Int128.TryParse(String.Empty, out _));
            Assert.ThrowsException<FormatException>(() => Int128.Parse(String.Empty));

            Assert.IsFalse(Int128.TryParse("1x2", out _));
            Assert.ThrowsException<FormatException>(() => Int128.Parse("-1 2"));

            // Max is parsable, but one more isn't.
            Assert.IsTrue(Int128.TryParse(Int128.MaxValue.ToString(), out Int128 max));
            Assert.IsTrue(max == Int128.MaxValue);
            string oneOverMax = UInt128Test.MakeNumberOneBigger(Int128.MaxValue.ToString());
            Assert.IsFalse(Int128.TryParse(oneOverMax, out _));
            Assert.ThrowsException<OverflowException>(() => Int128.Parse(oneOverMax));

            // Min is parsable, but one less isn't.
            Assert.IsTrue(Int128.TryParse(Int128.MinValue.ToString(), out Int128 min));
            Assert.IsTrue(min == Int128.MinValue);
            string oneUnderMin = UInt128Test.MakeNumberOneBigger(Int128.MinValue.ToString());
            Assert.IsFalse(Int128.TryParse(oneUnderMin, out _));
            Assert.ThrowsException<OverflowException>(() => Int128.Parse(oneUnderMin));

        }

        [TestMethod]
        public void Int128Int64Agreement () {

            Random rng = new Random(12864);
            foreach (long u in GetRandomInt64(8, rng)) {
                foreach (long v in GetRandomInt64(8, rng)) {

                    Int128 sum = ((Int128) u) + ((Int128) v);
                    long sum1 = u + v;
                    Assert.IsTrue((long) sum == sum1);

                    Int128 difference = ((Int128) u) - ((Int128) v);
                    long difference1 = u - v;
                    Assert.IsTrue((long) difference == difference1);

                    Int128 product = ((Int128) u) * ((Int128) v);
                    long product1 = u * v;
                    Assert.IsTrue((long) product == product1);

                    Int128 quotient = ((Int128) u) / ((Int128) v);
                    long quotient1 = u / v;
                    Assert.IsTrue((long) quotient == quotient);

                }
            }

        }

        [TestMethod]
        public void Int128BigIntegerAgreement () {
            Random rng = new Random(3);
            foreach (Int128 u in GetRandomInt128(10, rng)) {
                BigInteger u1 = BigInteger.Parse(u.ToString());
                foreach (Int128 v in GetRandomInt128(10, rng)) {
                    BigInteger v1 = BigInteger.Parse(v.ToString());

                    // Comparison
                    Assert.IsTrue((u < v) == (u1 < v1));
                    Assert.IsTrue((u > v) == (u1 > v1));

                    // To discard possible higher order bits in BigInteger results,
                    // coerce the result into an Int128.

                    // Sum
                    Int128 s = u + v;
                    BigInteger s1 = u1 + v1;
                    Assert.IsTrue((Int128) s1 == s);

                    // Difference
                    Int128 d = u - v;
                    BigInteger d1 = u1 - v1;
                    Assert.IsTrue((Int128) d1 == d);

                    // Product
                    Int128 p = u * v;
                    BigInteger p1 = u1 * v1;
                    Assert.IsTrue((Int128) p1 == p);

                    // Quotient
                    // Since results are guaranteed to fit in Int128, we can safely compare to BigInteger directly.
                    Int128 q = Int128.DivRem(u, v, out Int128 r);
                    BigInteger q1 = BigInteger.DivRem(u1, v1, out BigInteger r1);
                    Assert.IsTrue((BigInteger) q == q1);
                    Assert.IsTrue((BigInteger) r == r1);

                }
            }

        }

        [TestMethod]
        public void Int128DoubleRoundtrip () {
            // If value is less than 2^52, we should be able to roundtrip to double.
            foreach (long x in GetRandomInt64(8)) {
                if (Math.Abs(x) < ((long) (1UL << 52))) {
                    Int128 y = (Int128) x;
                    Assert.IsTrue(((Int128) ((double) y)) == y);
                }
            }
        }

        [TestMethod]
        public void Int128Abs () {
            Assert.IsTrue(Int128.Abs(Int128.Zero) == Int128.Zero);
            Assert.ThrowsException<OverflowException>(() => Int128.Abs(Int128.MinValue));

            foreach (Int128 x in GetRandomInt128(8)) {
                Int128 y = Int128.Abs(x);
                if (x < Int128.Zero) {
                    Assert.IsTrue(y == -x);
                } else {
                    Assert.IsTrue(y == x);
                }
            }
        }

        [TestMethod]
        public void Int128TimingComparison () {

            Int128[] values = GetRandomInt128(1000).ToArray();
            TimeOperations<Int128>("Addition", (x, y) => { Int128 z = x + y; }, values);
            TimeOperations<Int128>("Subtraction", (x, y) => { Int128 z = x - y; }, values);
            TimeOperations<Int128>("Multiplication", (x, y) => { Int128 z = x * y; }, values);
            TimeOperations<Int128>("Division", (x, y) => { Int128 z = x / y; }, values);

            BigInteger[] bigValues = values.Select(x => BigInteger.Parse(x.ToString())).ToArray();
            TimeOperations<BigInteger>("Big Addition", (x, y) => { BigInteger z = x + y; }, bigValues);
            TimeOperations<BigInteger>("Big Subtraction", (x, y) => { BigInteger z = x - y; }, bigValues);
            TimeOperations<BigInteger>("Big Multiplication", (x, y) => { BigInteger z = x * y; }, bigValues);
            TimeOperations<BigInteger>("Big Division", (x, y) => { BigInteger z = x / y; }, bigValues);

            //Decimal[] dValues = values.Select(x => Decimal.Parse(x.ToString())).ToArray();
            //TimeOperations<Decimal>("Decimal Addition", (x, y) => { Decimal d = x + y; }, dValues);
            //TimeOperations<Decimal>("Decimal Multiplication", (x, y) => { Decimal d = x * y; }, dValues);
            //TimeOperations<Decimal>("Decimal Division", (x, y) => { Decimal d = x / y; }, dValues);

            long[] longValues = values.Select(x => (long) x).ToArray();
            TimeOperations<long>("Long Addition", (x, y) => { long z = x + y; }, longValues);
            TimeOperations<long>("Long Multiplication", (x, y) => { long z = x * y; }, longValues);
            TimeOperations<long>("Long Division", (x, y) => { long z = x / y; }, longValues);

        }

        private void TimeOperations<T>(string name, Action<T,T> operation, IEnumerable<T> values) {
            Stopwatch t = Stopwatch.StartNew();
            foreach (T x in values) {
                foreach (T y in values) {
                    operation(x, y);
                }
            }
            t.Stop();
            Console.WriteLine($"{name}: {t.ElapsedMilliseconds}");
        }

    }
}
