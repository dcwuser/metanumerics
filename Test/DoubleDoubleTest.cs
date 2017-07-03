using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Extended;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public class DoubleDoubleTest {

        public bool IsNearlyEqual (DoubleDouble a, DoubleDouble b) {
            DoubleDouble magnitude = DoubleDouble.Abs(a) + DoubleDouble.Abs(b);
            DoubleDouble difference = a - b;
            return (DoubleDouble.Abs(difference) <= magnitude * 1.0E-30);
        }

        IEnumerable<DoubleDouble> GetRandomPositiveDoubleDoubles (double minimum, double maximum, int n, Random rng = null) {
            if (rng == null) rng = new Random(1);
            double logMinimum = Math.Log(minimum);
            double logMaximum = Math.Log(maximum);
            for (int i = 0; i < n; i++) {
                DoubleDouble x = DoubleDouble.GetRandomValue(rng);
                DoubleDouble logValue = logMinimum * (1.0 - x) + logMaximum * x;
                DoubleDouble value = DoubleDouble.Exp(logValue);
                yield return (value);
            }
        }

        [TestMethod]
        public void DoubleDoubleArithmeticAxioms () {

            Random rng = new Random(2);
            foreach (DoubleDouble a in GetRandomPositiveDoubleDoubles(1.0E-4, 1.0E4, 8, rng)) {

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

                foreach (DoubleDouble b in GetRandomPositiveDoubleDoubles(1.0E-4, 1.0E4, 4, rng)) {

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

            DoubleDouble x = DoubleDouble.E;
            DoubleDouble y = DoubleDouble.Pi;

            Assert.IsTrue(x == x);
            Assert.IsTrue(DoubleDouble.Equals(x, x));
            Assert.IsTrue(x.Equals(x));

            Assert.IsTrue(x != y);
            Assert.IsFalse(DoubleDouble.Equals(x, y));
            Assert.IsFalse(x.Equals(y));
            Assert.IsFalse(x.Equals(null));

            Assert.IsTrue(x.GetHashCode() != y.GetHashCode());
        }

        [TestMethod]
        public void DoubleDoubleComparison () {

            DoubleDouble x = DoubleDouble.E;
            DoubleDouble y = DoubleDouble.Pi;

            Assert.IsTrue(x <= x);
            Assert.IsTrue(x >= x);
            Assert.IsFalse(x < x);
            Assert.IsFalse(x > x);
            Assert.IsTrue(x.CompareTo(x) == 0);

            Assert.IsTrue(x < y);
            Assert.IsTrue(x <= y);
            Assert.IsFalse(x > y);
            Assert.IsFalse(x >= y);
            Assert.IsTrue(x.CompareTo(y) == -1);
        }

        [TestMethod]
        public void DoubleDoubleSerialization () {

            foreach (DoubleDouble x in GetRandomPositiveDoubleDoubles(1.0E-12, 1.0E12, 24)) {

                string s = x.ToString();
                Console.WriteLine(s);

                DoubleDouble y = DoubleDouble.Parse(s);

                Assert.IsTrue(IsNearlyEqual(x, y));

            }

        }

        [TestMethod]
        public void DoubleDoubleParseSpecialCases () {
            Assert.IsTrue(DoubleDouble.Parse("0.0") == DoubleDouble.Zero);
            Assert.IsTrue(DoubleDouble.Parse("1.0") == DoubleDouble.One);
            Assert.IsTrue(DoubleDouble.Parse("-1.0") == -DoubleDouble.One);
        }

        [TestMethod]
        public void DoubleDoubleSqrtSpecialValues () {
            Assert.IsTrue(DoubleDouble.Sqrt(DoubleDouble.Zero) == DoubleDouble.Zero);
            Assert.IsTrue(DoubleDouble.Sqrt(DoubleDouble.One) == DoubleDouble.One);
        }

        [TestMethod]
        public void DoubleDoubleSqrt () {
            foreach (DoubleDouble x in GetRandomPositiveDoubleDoubles(1.0E-4, 1.0E4, 16)) {
                DoubleDouble y = DoubleDouble.Sqrt(x);
                Assert.IsTrue(IsNearlyEqual(y * y, x));
            }
        }

        [TestMethod]
        public void DoubleDoubleSqrtAgreement () {
            foreach (DoubleDouble x in GetRandomPositiveDoubleDoubles(1.0E-4, 1.0E4, 16)) {
                DoubleDouble xSqrt = DoubleDouble.Sqrt(x);
                double xSqrtAsDouble = (double) xSqrt;

                double xAsDouble = (double) x;
                double xAsDoubleSqrt = Math.Sqrt(xAsDouble);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(xSqrtAsDouble, xAsDoubleSqrt));
            }
        }

        [TestMethod]
        public void DoubleDoubleLogExpSpecialValues () {
            Assert.IsTrue(DoubleDouble.Log(DoubleDouble.One) == DoubleDouble.Zero);
            Assert.IsTrue(DoubleDouble.Exp(DoubleDouble.Zero) == DoubleDouble.One);

            Assert.IsTrue(IsNearlyEqual(DoubleDouble.Exp(DoubleDouble.One), DoubleDouble.E));

            Assert.IsTrue(IsNearlyEqual(DoubleDouble.Exp(-DoubleDouble.One), DoubleDouble.One / DoubleDouble.E));
        }

        [TestMethod]
        public void DoubleDoubleLogExp () {
            foreach (DoubleDouble x in GetRandomPositiveDoubleDoubles(1.0E-2, 1.0E6, 16)) {
                DoubleDouble xLog = DoubleDouble.Log(x);
                DoubleDouble xLogExp = DoubleDouble.Exp(xLog);
                Assert.IsTrue(IsNearlyEqual(xLogExp, x));
            }
        }

        [TestMethod]
        public void DoubleDoubleErrorFunctionSpecialCases () {
            Assert.IsTrue(AdvancedDoubleDoubleMath.Erf(DoubleDouble.Zero) == DoubleDouble.Zero);
            Assert.IsTrue(AdvancedDoubleDoubleMath.Erfc(DoubleDouble.Zero) == DoubleDouble.One);
        }

        [TestMethod]
        public void DoubleDoubleErrorFunctionReflection () {
            foreach (DoubleDouble x in GetRandomPositiveDoubleDoubles(1.0E-4, 1.0E4, 8)) {
                DoubleDouble xErf = AdvancedDoubleDoubleMath.Erf(x);
                Assert.IsTrue(-AdvancedDoubleDoubleMath.Erf(-x) == AdvancedDoubleDoubleMath.Erf(x));
            }
        }

        [TestMethod]
        public void DoubleDoubleErrorFunctionComplementiarity () {
            foreach (DoubleDouble x in GetRandomPositiveDoubleDoubles(1.0E-2, 1.0E2, 8)) {
                DoubleDouble xErf = AdvancedDoubleDoubleMath.Erf(x);
                DoubleDouble xErfc = AdvancedDoubleDoubleMath.Erfc(x);
                Assert.IsTrue(IsNearlyEqual(xErf + xErfc, DoubleDouble.One));
            }
        }

        [TestMethod]
        public void DoubleDoubleErrorFunctionAgreement () {
            foreach (DoubleDouble x in GetRandomPositiveDoubleDoubles(1.0E-2, 1.0E2, 8)) {

                DoubleDouble xErf = AdvancedDoubleDoubleMath.Erf(x);
                double xErfAsDouble = (double) xErf;

                double xAsDouble = (double) x;
                double xAsDoubleErf = AdvancedMath.Erf(xAsDouble);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(xErfAsDouble, xAsDoubleErf));

                DoubleDouble xErfc = AdvancedDoubleDoubleMath.Erfc(x);
                Assert.IsTrue(IsNearlyEqual(xErf + xErfc, DoubleDouble.One));

            }
        }

    }
}
