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
        public void BigHypotenuseTest () {
            
            // construct 3-4-5 right triangle with side lengths whoose squares would overflow
            double x = Double.MaxValue / 8.0;

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                MoreMath.Hypot(3.0 * x, 4.0 * x), 5.0 * x
            ));

        }


        [TestMethod]
        public void SmallHypotenuseTest () {

            // construct 3-4-5 right triangle with side lengths whoose squares would underflow
            double x = 8.0 / Double.MaxValue;

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                MoreMath.Hypot(3.0 * x, 4.0 * x), 5.0 * x
            ));
        }
    }

}