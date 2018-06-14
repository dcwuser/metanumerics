using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Functions;

namespace Test {
    public partial class AdvancedMathTest_Gamma {

        [TestMethod]
        public void LogGammaSpecialValues () {
            Assert.IsTrue(Double.IsPositiveInfinity(AdvancedMath.LogGamma(0.0)));
            Assert.IsTrue(AdvancedMath.LogGamma(1.0) == 0.0);
            Assert.IsTrue(AdvancedMath.LogGamma(2.0) == 0.0);
        }

        [TestMethod]
        public void LogGammaDuplication () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 24)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.LogGamma(2.0 * x),
                    AdvancedMath.LogGamma(x) + AdvancedMath.LogGamma(x + 0.5) + (2.0 * x - 0.5) * Math.Log(2.0) - 0.5 * Math.Log(2.0 * Math.PI)
                ));
            }
        }

        [TestMethod]
        public void LogGammaTriplication () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 24)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.LogGamma(3.0 * x),
                    AdvancedMath.LogGamma(x) + AdvancedMath.LogGamma(x + 1.0 / 3.0) + AdvancedMath.LogGamma(x + 2.0 / 3.0) + (3.0 * x - 0.5) * Math.Log(3.0) - Math.Log(2.0 * Math.PI)
                ));
            }
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void LogGammaNegativeArgument () {
            AdvancedMath.LogGamma(-0.5);
        }

        [TestMethod]
        public void GammaRatioInequality () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 16)) {
                double LB = Math.Log(Math.Sqrt(x + 0.25));
                double UB = Math.Log((x + 0.5) / Math.Sqrt(x + 0.75));
                double R = AdvancedMath.LogGamma(x + 1.0) - AdvancedMath.LogGamma(x + 0.5);
                Assert.IsTrue(LB <= R);
                Assert.IsTrue(R <= UB);
            }
        }

    }
}
