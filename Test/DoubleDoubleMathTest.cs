using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Extended;
using Meta.Numerics.Functions;
namespace Test {

    [TestClass]
    public class DoubleDoubleMathTest {

        [TestMethod]
        public void DoubleDoubleSqrtSpecialValues () {
            Assert.IsTrue(DoubleDouble.IsNaN(DoubleDouble.Sqrt(-DoubleDouble.One)));
            Assert.IsTrue(DoubleDouble.Sqrt(DoubleDouble.Zero) == DoubleDouble.Zero);
            Assert.IsTrue(DoubleDouble.Sqrt(DoubleDouble.One) == DoubleDouble.One);
            Assert.IsTrue(DoubleDouble.Sqrt(DoubleDouble.PositiveInfinity) == DoubleDouble.PositiveInfinity);
        }

        [TestMethod]
        public void DoubleDoubleSqrt () {
            foreach (DoubleDouble x in DoubleDoubleTest.GetRandomDoubleDoubles(1.0E-4, 1.0E4, 16)) {
                DoubleDouble y = DoubleDouble.Sqrt(x);
                Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(y * y, x));
            }
        }

        [TestMethod]
        public void DoubleDoubleSqrtAgreement () {
            foreach (DoubleDouble x in DoubleDoubleTest.GetRandomDoubleDoubles(1.0E-4, 1.0E4, 16)) {
                DoubleDouble xSqrt = DoubleDouble.Sqrt(x);
                double xSqrtAsDouble = (double) xSqrt;

                double xAsDouble = (double) x;
                double xAsDoubleSqrt = Math.Sqrt(xAsDouble);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(xSqrtAsDouble, xAsDoubleSqrt));
            }
        }

        [TestMethod]
        public void DoubleDoubleLogExpSpecialValues () {
            Assert.IsTrue(DoubleDouble.IsNaN(DoubleDouble.Log(DoubleDouble.NaN)));
            Assert.IsTrue(DoubleDouble.Log(DoubleDouble.Zero) == DoubleDouble.NegativeInfinity);
            Assert.IsTrue(DoubleDouble.Log(DoubleDouble.One) == DoubleDouble.Zero);
            Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(DoubleDouble.Log(DoubleDouble.E), DoubleDouble.One));
            Assert.IsTrue(DoubleDouble.Log(DoubleDouble.PositiveInfinity) == DoubleDouble.PositiveInfinity);

            Assert.IsTrue(DoubleDouble.Exp(DoubleDouble.Zero) == DoubleDouble.One);
            Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(DoubleDouble.Exp(DoubleDouble.One), DoubleDouble.E));
            Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(DoubleDouble.Exp(-DoubleDouble.One), DoubleDouble.One / DoubleDouble.E));
        }

        [TestMethod]
        public void DoubleDoubleLogExp () {
            foreach (DoubleDouble x in DoubleDoubleTest.GetRandomDoubleDoubles(1.0E-2, 1.0E6, 16)) {
                DoubleDouble xLog = DoubleDouble.Log(x);
                DoubleDouble xLogExp = DoubleDouble.Exp(xLog);
                Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(xLogExp, x));
            }
        }

        [TestMethod]
        public void DoubleDoubleTrigSpecialCases () {
            Assert.IsTrue(DoubleDouble.Sin(0.0) == DoubleDouble.Zero);
            Assert.IsTrue(DoubleDouble.Cos(0.0) == DoubleDouble.One);

            Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(DoubleDouble.Sin(DoubleDouble.Pi / 6), DoubleDouble.One / 2));
            Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(DoubleDouble.Cos(DoubleDouble.Pi / 6), DoubleDouble.Sqrt(3) / 2));

            Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(DoubleDouble.Cos(DoubleDouble.Pi / 5), (1 + DoubleDouble.Sqrt(5)) / 4));

            Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(DoubleDouble.Sin(DoubleDouble.Pi / 4), DoubleDouble.One / DoubleDouble.Sqrt(2)));
            Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(DoubleDouble.Cos(DoubleDouble.Pi / 4), DoubleDouble.One / DoubleDouble.Sqrt(2)));

            Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(DoubleDouble.Cos(DoubleDouble.Pi / 3), DoubleDouble.One / 2));
            Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(DoubleDouble.Sin(DoubleDouble.Pi / 3), DoubleDouble.Sqrt(3) / 2));

            Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(DoubleDouble.Sin(DoubleDouble.Pi / 2), DoubleDouble.One));

            Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(DoubleDouble.Cos(DoubleDouble.Pi), -DoubleDouble.One));

            // Not testing zeros because of need to distingish very small values
        }

        [TestMethod]
        public void DoubleDoubleTrigSymmetries () {
            foreach (DoubleDouble x in DoubleDoubleTest.GetRandomDoubleDoubles(1.0E-3, 1.0E3, 8)) {
                Assert.IsTrue(DoubleDouble.Sin(x) == -DoubleDouble.Sin(-x));
                Assert.IsTrue(DoubleDouble.Cos(x) == DoubleDouble.Cos(-x));
            }
        }

        [TestMethod]
        public void DoubleDoubleTrigIdentities () {
            Random rng = new Random(1);
            for (int i = 0; i < 16; i++) {
                DoubleDouble x = (2.0 * DoubleDouble.Pi) * (1.0 - 2.0 * DoubleDouble.GetRandomValue(rng));

                DoubleDouble s = DoubleDouble.Sin(x);
                DoubleDouble c = DoubleDouble.Cos(x);
                Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(s * s + c * c, DoubleDouble.One));

                DoubleDouble s2 = DoubleDouble.Sin(2.0 * x);
                Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(s2, 2.0 * s * c));

            }
        }

    }
}
