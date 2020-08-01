using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Test {

    /// <summary>
    ///This is a test class for AdvancedMathTest and is intended
    ///to contain all AdvancedMathTest Unit Tests
    ///</summary>
    [TestClass()]
    public class MoreMathTest {

        private TestContext testContextInstance;

        /// <summary>
        ///Gets or sets the test context which provides
        ///information about and functionality for the current test run.
        ///</summary>
        public TestContext TestContext {
            get {
                return testContextInstance;
            }
            set {
                testContextInstance = value;
            }
        }

        [TestMethod]
        public void IntegerPowerTest () {

            for (int k = -20; k <= 20; k++) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    MoreMath.Pow(Math.PI, k), Math.Pow(Math.PI, (double) k)
                ));
            }

        }

        [TestMethod]
        public void IntegerPowerSpecialCases () {

            double x = 2.7;

            Assert.IsTrue(MoreMath.Pow(x, 0) == 1.0);
            Assert.IsTrue(MoreMath.Pow(x, 1) == x);
            Assert.IsTrue(MoreMath.Pow(x, 2) == x * x);

            Assert.IsTrue(MoreMath.Pow(0.0, 2) == 0.0);
            Assert.IsTrue(MoreMath.Pow(0.0, 1) == 0.0);
            Assert.IsTrue(MoreMath.Pow(0.0, 0) == 1.0);

            // By convention, the power limit is taken before the argument limit, so 0^0 = 1
        }

        [TestMethod]
        public void HypotBig () {
            
            // construct 3-4-5 right triangle with side lengths whoose squares would overflow
            double x = Double.MaxValue / 8.0;

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                MoreMath.Hypot(3.0 * x, 4.0 * x), 5.0 * x,
            TestUtilities.RelativeTarget));

        }


        [TestMethod]
        public void HypotSmall () {

            // construct 3-4-5 right triangle with side lengths whoose squares would underflow
            double x = 8.0 / Double.MaxValue;

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                MoreMath.Hypot(3.0 * x, 4.0 * x), 5.0 * x,
            TestUtilities.RelativeTarget));
        }

        [TestMethod]
        public void HypotSpecialCases () {

            double x = 3.14;

            Assert.IsTrue(MoreMath.Hypot(0.0, 0.0) == 0.0);
            Assert.IsTrue(MoreMath.Hypot(0.0, x) == Math.Abs(x));
            Assert.IsTrue(MoreMath.Hypot(0.0, -x) == Math.Abs(x));
            Assert.IsTrue(MoreMath.Hypot(x, 0.0) == Math.Abs(x));
            Assert.IsTrue(MoreMath.Hypot(-x, 0.0) == Math.Abs(x));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(MoreMath.Hypot(x, -x), Math.Sqrt(2.0) * Math.Abs(x)));
        }

        [TestMethod]
        public void HypotExtremeValues () {
            Assert.IsTrue(Double.IsPositiveInfinity(MoreMath.Hypot(Double.MaxValue, Double.MaxValue)));
            Assert.IsTrue(Double.IsPositiveInfinity(MoreMath.Hypot(-Double.MaxValue, Double.MaxValue)));
            Assert.IsTrue(Double.IsPositiveInfinity(MoreMath.Hypot(Double.PositiveInfinity, Double.NegativeInfinity)));
            Assert.IsTrue(Double.IsPositiveInfinity(MoreMath.Hypot(Double.NegativeInfinity, Double.NegativeInfinity)));
            Assert.IsTrue(Double.IsNaN(MoreMath.Hypot(Double.NaN, 0.0)));
            Assert.IsTrue(Double.IsNaN(MoreMath.Hypot(Double.PositiveInfinity, Double.NaN)));
        }

        [TestMethod]
        public void BigTrig () {
            // An exactly representable double that is extremely close to a multiple of \pi/2
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                MoreMath.Cos(6381956970095103.0 * MoreMath.Pow(2.0, 797)), -4.6871659242546276E-19,
            TestUtilities.RelativeTarget));
        }

        // https://www.csee.umbc.edu/~phatak/645/supl/Ng-ArgReduction.pdf says
        //   6381956970095103 X 2^{797} = n (\pi / 2) + 4.687165924254624 X 10^{-19}
        // This is supposedly the closest any representable double gets to an exact multiple of \pi / 2.

        // http://perso.ens-lyon.fr/nathalie.revol/publis/BDKMR05.pdf says
        // 6411027962775774 X 2^{-48} = n (\pi / 4) + 3.094903 X 10^{-19}


        [TestMethod]
        public void SinPiLargeIntegers () {
            // For large integers, SinPi should be exactly zero
            foreach(int x in TestUtilities.GenerateIntegerValues(100, Int32.MaxValue, 8)) {
                Assert.IsTrue(MoreMath.SinPi((double) x) == 0.0);
                Assert.IsTrue(MoreMath.SinPi((double) -x) == 0.0);
            }
        }

        [TestMethod]
        public void SinPiSmallArguments () {
            // For small x, SinPi(x) and Sin(x * PI) should agree
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 4)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(Math.Sin(x * Math.PI), MoreMath.SinPi(x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(Math.Sin(-x * Math.PI), MoreMath.SinPi(-x)));
            }
        }

        [TestMethod]
        public void CosPiLargeIntegers () {
            // For large integers, CosPi should be exactly plus or minus one.
            foreach (int x in TestUtilities.GenerateIntegerValues(100, Int32.MaxValue, 8)) {
                Assert.IsTrue(MoreMath.CosPi((double) x) == (x % 2 == 0 ? 1 : -1));
                Assert.IsTrue(MoreMath.CosPi((double) -x) == (x % 2 == 0 ? 1 : -1));
                // In this range, the addition of 0.5 is exact, so CosPi at half integers should
                // be exactly zero.
                Assert.IsTrue(MoreMath.CosPi(((double) x) + 0.5) == 0.0);
            }
        }

        [TestMethod]
        public void CosPiSmallArguments () {
            // For small x, CosPi(x) and Cos(x * PI) should agree
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 4)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(Math.Cos(x * Math.PI), MoreMath.CosPi(x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(Math.Cos(-x * Math.PI), MoreMath.CosPi(-x)));
            }
        }

        [TestMethod]
        public void ExpMinusOneSpecialCases () {
            Assert.IsTrue(MoreMath.ExpMinusOne(0.0) == 0.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(MoreMath.ExpMinusOne(1.0), Math.E - 1.0));
        }

        [TestMethod]
        public void ExpMinusOneExtremeValues () {
            Assert.IsTrue(MoreMath.ExpMinusOne(Double.NegativeInfinity) == -1.0);
            Assert.IsTrue(MoreMath.ExpMinusOne(Double.MinValue) == -1.0);
            Assert.IsTrue(MoreMath.ExpMinusOne(-Double.Epsilon) == -Double.Epsilon);
            Assert.IsTrue(MoreMath.ExpMinusOne(Double.Epsilon) == Double.Epsilon);
            Assert.IsTrue(MoreMath.ExpMinusOne(Double.PositiveInfinity) == Double.PositiveInfinity);
            Assert.IsTrue(MoreMath.ExpMinusOne(Double.MaxValue) == Double.PositiveInfinity);
            Assert.IsTrue(Double.IsNaN(MoreMath.ExpMinusOne(Double.NaN)));
        }

        [TestMethod]
        public void ExpMinusOneExpAgreement () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E2, 8)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(MoreMath.ExpMinusOne(x) + 1.0, Math.Exp(x)));
            }
        }

        [TestMethod]
        public void LogOnePlusSpecialCases () {
            Assert.IsTrue(MoreMath.LogOnePlus(0.0) == 0.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(MoreMath.LogOnePlus(1.0), Math.Log(2.0)));

        }

        [TestMethod]
        public void LogOnePlusExtremeValues () {
            Assert.IsTrue(MoreMath.LogOnePlus(-1.0) == Double.NegativeInfinity);
            Assert.IsTrue(MoreMath.LogOnePlus(-Double.Epsilon) == -Double.Epsilon);
            Assert.IsTrue(MoreMath.LogOnePlus(Double.Epsilon) == Double.Epsilon);
            Assert.IsTrue(Double.IsNaN(MoreMath.LogOnePlus(Double.NaN)));
        }

        [TestMethod]
        public void LogOnePlusAgreement () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0, 1.0E2, 4)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(MoreMath.LogOnePlus(x), Math.Log(1.0 + x)));
            }
        }

    }

}