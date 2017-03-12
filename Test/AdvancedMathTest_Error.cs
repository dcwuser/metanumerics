using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_Error {

        [TestMethod]
        public void ErrorFunctionSpecialCases () {
            Assert.IsTrue(AdvancedMath.Erf(0.0) == 0.0);
            Assert.IsTrue(AdvancedMath.Erfc(0.0) == 1.0);
        }


        [TestMethod]
        public void ErrorFunctionExtremeValues () {

            Assert.IsTrue(AdvancedMath.Erf(Double.MaxValue) == 1.0);
            Assert.IsTrue(AdvancedMath.Erfc(Double.MaxValue) == 0.0);

            Assert.IsTrue(AdvancedMath.Erf(Double.PositiveInfinity) == 1.0);
            Assert.IsTrue(AdvancedMath.Erfc(Double.PositiveInfinity) == 0.0);

            Assert.IsTrue(AdvancedMath.Erf(-Double.MaxValue) == -1.0);

            Assert.IsTrue(AdvancedMath.Erf(Double.NegativeInfinity) == -1.0);

            Assert.IsTrue(Double.IsNaN(AdvancedMath.Erf(Double.NaN)));
            Assert.IsTrue(Double.IsNaN(AdvancedMath.Erfc(Double.NaN)));

        }

        [TestMethod]
        public void ErfReflection () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 4)) {
                Assert.IsTrue(AdvancedMath.Erf(-x) == -AdvancedMath.Erf(x));
            }
        }

        [TestMethod]
        public void ErrorFunctionComplementarity () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 8)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Erf(x) + AdvancedMath.Erfc(x), 1.0));
            }
        }

        [TestMethod]
        public void ErfcInequality () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 16)) {
                double erfc = AdvancedMath.Erfc(x);
                double factor = 2.0 / Math.Sqrt(Math.PI) * Math.Exp(-x * x);
                double lower = factor / (x + Math.Sqrt(x * x + 2.0));
                double upper = factor / (x + Math.Sqrt(x * x + 4.0 / Math.PI));
                Assert.IsTrue((lower <= erfc) && (erfc <= upper));
            }

        }

        [TestMethod]
        public void InverseErfSpecialValues () {

            Assert.IsTrue(AdvancedMath.InverseErf(-1.0) == Double.NegativeInfinity);
            Assert.IsTrue(AdvancedMath.InverseErf(0.0) == 0.0);
            Assert.IsTrue(AdvancedMath.InverseErf(1.0) == Double.PositiveInfinity);

            //Assert.IsTrue(Double.IsNaN(AdvancedMath.InverseErf(Double.NaN)));

        }

        [TestMethod]
        public void InverseErfcSpecialCases () {
            Assert.IsTrue(AdvancedMath.InverseErfc(0.0) == Double.PositiveInfinity);
            Assert.IsTrue(AdvancedMath.InverseErfc(1.0) == 0.0);
        }

        [TestMethod]
        public void InverseErfTest () {
            foreach (double P in TestUtilities.GenerateRealValues(1.0E-8, 1.0, 16)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Erf(AdvancedMath.InverseErf(P)), P));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Erfc(AdvancedMath.InverseErfc(P)), P));
            }
        }

    }
}
