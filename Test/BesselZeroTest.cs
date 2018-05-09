using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Analysis;

namespace Test {

    [TestClass]
    public class BesselZeroTest {

        [TestMethod]
        public void AiryZeros () {
            foreach (int k in new int[] { 1, 10, 100, 1000 }) {
                double a = AdvancedMath.AiryAiZero(k);
                double ya = AdvancedMath.AiryAi(a);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ya, 0.0, 2.0E-16 * k));

                double b = AdvancedMath.AiryBiZero(k);
                double yb = AdvancedMath.AiryBi(b);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(yb, 0.0, 2.0E-16 * k));
            }
        }
    }
}
