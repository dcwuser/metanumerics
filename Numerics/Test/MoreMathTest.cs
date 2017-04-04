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

    }

}