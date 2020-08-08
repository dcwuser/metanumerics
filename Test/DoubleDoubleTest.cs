using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Extended;

namespace Test {

    [TestClass]
    public class DoubleDoubleTest {

        public static bool IsNearlyEqual (DoubleDouble a, DoubleDouble b) {
            DoubleDouble magnitude = DoubleDouble.Abs(a) + DoubleDouble.Abs(b);
            DoubleDouble difference = a - b;
            return (DoubleDouble.Abs(difference) <= magnitude * 1.0E-30);
        }

        public static IEnumerable<DoubleDouble> GetRandomDoubleDoubles (double minimum, double maximum, int n, bool negative = false, Random rng = null) {
            if (rng == null) rng = new Random(1);
            double logMinimum = Math.Log(minimum);
            double logMaximum = Math.Log(maximum);
            for (int i = 0; i < n; i++) {
                DoubleDouble x = DoubleDouble.GetRandomValue(rng);
                DoubleDouble logValue = logMinimum * (1.0 - x) + logMaximum * x;
                DoubleDouble value = DoubleDouble.Exp(logValue);
                if (negative && rng.NextDouble() < 0.5) value = -value;
                yield return (value);
            }
        }

        [TestMethod]
        public void DoubleDoubleArithmeticAxioms () {

            Random rng = new Random(2);
            foreach (DoubleDouble a in GetRandomDoubleDoubles(1.0E-4, 1.0E4, 8, true, rng)) {

                Assert.IsTrue(a + DoubleDouble.Zero == a);
                Assert.IsTrue(DoubleDouble.Zero + a == a);

                Assert.IsTrue(a * DoubleDouble.One == a);
                Assert.IsTrue(DoubleDouble.One * a == a);

                DoubleDouble ma = -a;

                Assert.IsTrue(DoubleDouble.Zero - a == ma);
                Assert.IsTrue(-DoubleDouble.One * a == ma);

                Assert.IsTrue(a + ma == DoubleDouble.Zero);
                Assert.IsTrue(ma + a == DoubleDouble.Zero);

                Assert.IsTrue(a + a == 2.0 * a);

                foreach (DoubleDouble b in GetRandomDoubleDoubles(1.0E-4, 1.0E4, 4, true, rng)) {

                    Assert.IsTrue(a + b == b + a);
                    Assert.IsTrue(a * b == b * a);

                    Assert.IsTrue(a - b == -(b - a));

                    DoubleDouble a_over_b = a / b;
                    DoubleDouble b_over_a = b / a;
                    Assert.IsTrue(IsNearlyEqual(a_over_b * b, a));
                    Assert.IsTrue(IsNearlyEqual(b_over_a * a, b));
                    Assert.IsTrue(IsNearlyEqual(DoubleDouble.One / a_over_b, b_over_a));

                }
            }

        }

        [TestMethod]
        public void DoubleDoubleEquality () {
            DoubleDouble[] values = new DoubleDouble[] { DoubleDouble.NegativeInfinity, Double.MinValue, -DoubleDouble.One, DoubleDouble.Zero, DoubleDouble.E, DoubleDouble.Pi, Double.MaxValue, DoubleDouble.PositiveInfinity};
            for (int i = 0; i < values.Length; i++) {
                Assert.IsFalse(values[i].Equals(null));
                for (int j = 0; j < values.Length; j++) {
                    if (i == j) {
                        Assert.IsTrue(values[i] == values[j]);
                        Assert.IsFalse(values[i] != values[j]);
                        Assert.IsTrue(DoubleDouble.Equals(values[i], values[j]));
                        Assert.IsTrue(values[i].Equals(values[j]));
                        Assert.IsTrue(values[i].Equals(values[j] as object));
                        Assert.IsTrue(values[i].GetHashCode() == values[j].GetHashCode());
                    } else {
                        Assert.IsFalse(values[i] == values[j]);
                        Assert.IsTrue(values[i] != values[j]);
                        Assert.IsFalse(DoubleDouble.Equals(values[i], values[j]));
                        Assert.IsFalse(values[i].Equals(values[j]));
                        Assert.IsFalse(values[i].Equals(values[j] as object));
                        // Hash code inequality isn't strictly forbidden, but it should be so unlikely as to almost certainly not occur within our small set.
                        Assert.IsTrue(values[i].GetHashCode() != values[j].GetHashCode());
                    }
                }
            }
        }

        [TestMethod]
        public void DoubleDoubleComparison () {

            DoubleDouble[] orderedValues = new DoubleDouble[] { DoubleDouble.NegativeInfinity, -DoubleDouble.One, DoubleDouble.Zero, DoubleDouble.E, DoubleDouble.Pi, Double.MaxValue, DoubleDouble.PositiveInfinity };
            for (int i = 0; i < orderedValues.Length; i++) {
                for (int j = 0; j < orderedValues.Length; j++) {
                    if (i < j) {
                        Assert.IsTrue(DoubleDouble.Compare(orderedValues[i], orderedValues[j]) == -1);
                        Assert.IsTrue(orderedValues[i].CompareTo(orderedValues[j]) == -1);
                        Assert.IsTrue(orderedValues[i] < orderedValues[j]);
                        Assert.IsFalse(orderedValues[i] >= orderedValues[j]);
                    } else if (i == j) {
                        Assert.IsTrue(DoubleDouble.Compare(orderedValues[i], orderedValues[j]) == 0);
                        Assert.IsTrue(orderedValues[i].CompareTo(orderedValues[j]) == 0);
                        Assert.IsFalse(orderedValues[i] < orderedValues[j]);
                        Assert.IsFalse(orderedValues[i] > orderedValues[j]);
                        Assert.IsTrue(orderedValues[i] <= orderedValues[j]);
                    } else {
                        Assert.IsTrue(DoubleDouble.Compare(orderedValues[i], orderedValues[j]) == +1);
                        Assert.IsTrue(orderedValues[i].CompareTo(orderedValues[j]) == +1);
                        Assert.IsTrue(orderedValues[i] > orderedValues[j]);
                        Assert.IsFalse(orderedValues[i] <= orderedValues[j]);
                    }
                }
            }

        }

        [TestMethod]
        public void DoubleDoubleNaN () {

            // Equality
            Assert.IsFalse(DoubleDouble.NaN == DoubleDouble.NaN);
            Assert.IsTrue(DoubleDouble.NaN != DoubleDouble.NaN);
            Assert.IsTrue(DoubleDouble.IsNaN(DoubleDouble.NaN));

            // Comparison
            Assert.IsFalse(DoubleDouble.NaN < DoubleDouble.Zero);
            Assert.IsFalse(DoubleDouble.NaN <= DoubleDouble.Zero);
            Assert.IsFalse(DoubleDouble.NaN >= DoubleDouble.Zero);
            Assert.IsFalse(DoubleDouble.NaN > DoubleDouble.Zero);

        }

        [TestMethod]
        public void DoubleDoubleInfinity () {

            // Addition
            Assert.IsTrue(DoubleDouble.PositiveInfinity + DoubleDouble.One == DoubleDouble.PositiveInfinity);
            Assert.IsTrue(DoubleDouble.One + DoubleDouble.NegativeInfinity == DoubleDouble.NegativeInfinity);
            Assert.IsTrue(DoubleDouble.IsNaN(DoubleDouble.NegativeInfinity + DoubleDouble.PositiveInfinity));

            // Subtraction
            Assert.IsTrue(DoubleDouble.One - DoubleDouble.PositiveInfinity == DoubleDouble.NegativeInfinity);
            Assert.IsTrue(DoubleDouble.IsNaN(DoubleDouble.PositiveInfinity - DoubleDouble.PositiveInfinity));

            // Multiplication
            Assert.IsTrue(DoubleDouble.One * DoubleDouble.PositiveInfinity == DoubleDouble.PositiveInfinity);
            Assert.IsTrue(DoubleDouble.One * DoubleDouble.NegativeInfinity == DoubleDouble.NegativeInfinity);
            Assert.IsTrue(DoubleDouble.IsNaN(DoubleDouble.Zero * DoubleDouble.PositiveInfinity));

            // Division
            Assert.IsTrue(DoubleDouble.One / DoubleDouble.Zero == DoubleDouble.PositiveInfinity);
            Assert.IsTrue((-DoubleDouble.One) / DoubleDouble.Zero == DoubleDouble.NegativeInfinity);
            Assert.IsTrue(DoubleDouble.One / DoubleDouble.PositiveInfinity == DoubleDouble.Zero);
            Assert.IsTrue(DoubleDouble.One / DoubleDouble.NegativeInfinity == DoubleDouble.Zero);
            Assert.IsTrue(DoubleDouble.IsNaN(DoubleDouble.Zero / DoubleDouble.Zero));
            Assert.IsTrue(DoubleDouble.IsNaN(DoubleDouble.PositiveInfinity / DoubleDouble.PositiveInfinity));

        }

        [TestMethod]
        public void DoubleDoubleUnderflow () {
            DoubleDouble z = DoubleDouble.One;
            while (z != DoubleDouble.Zero) {
                z = z / 10.0;
            }
        }

        [TestMethod]
        public void DoubleDoubleSerializationRoundtrip () {
            Random rng = new Random(3);
            foreach (DoubleDouble x in GetRandomDoubleDoubles(1.0E-12, 1.0E12, 24, true, rng)) {
                string s = x.ToString();
                DoubleDouble y = DoubleDouble.Parse(s);
                Assert.IsTrue(IsNearlyEqual(x, y));
                bool r = DoubleDouble.TryParse(s, out DoubleDouble z);
                Assert.IsTrue(r);
                Assert.IsTrue(y == z);
            }

        }

        [TestMethod]
        public void DoubleDoubleToStringSpecial () {
            Assert.IsTrue(DoubleDouble.NaN.ToString() == Double.NaN.ToString());
        }

        [TestMethod]
        public void DoubleDoubleParseSpecialCases () {
            Assert.IsTrue(DoubleDouble.Parse("0.0") == DoubleDouble.Zero);
            Assert.IsTrue(DoubleDouble.Parse("1.0") == DoubleDouble.One);
            Assert.IsTrue(DoubleDouble.Parse("-1.0") == -DoubleDouble.One);
            Assert.IsTrue(DoubleDouble.IsNaN(DoubleDouble.Parse("NaN")));
        }

        [TestMethod]
        public void DoubleDoubleParseInvalid () {
            Assert.ThrowsException<ArgumentNullException>(() => DoubleDouble.Parse(null));
            Assert.ThrowsException<FormatException>(() => DoubleDouble.Parse(""));
            Assert.ThrowsException<FormatException>(() => DoubleDouble.Parse("-"));
            Assert.ThrowsException<FormatException>(() => DoubleDouble.Parse("."));
            Assert.ThrowsException<FormatException>(() => DoubleDouble.Parse("0.0.0"));
            Assert.ThrowsException<FormatException>(() => DoubleDouble.Parse("E0"));
            Assert.ThrowsException<FormatException>(() => DoubleDouble.Parse("1.0Eq"));

            Assert.IsFalse(DoubleDouble.TryParse(null, out _));
            Assert.IsFalse(DoubleDouble.TryParse("", out _));
        }

    }
}
