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

        [TestMethod]
        public void BesselJZeros () {
            foreach (double nu in new double[] { 0.0, 0.5, 1.0, 4.25, 16.0, 124.0 }) {
                foreach (int k in new int[] { 1, 2, 4, 64, 256 }) {

                    double j = AdvancedMath.BesselJZero(nu, k);
                    double J = AdvancedMath.BesselJ(nu, j);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(J, 0.0, 4.0E-16 * 3 * k));

                }
            }
        }

        [TestMethod]
        public void RayleighSums () {

            // Watson Chapter 15.51 gives expressions for the sum, due to Rayleigh
            //   \sigma^{(r)}{\nu} = \sum_{k=1}^{\infty} \frac{1}{j_{\nu, k}^{2r}}

            // Convergence is slower for for larger nu, so don't try with nu too big.

            foreach (double nu in TestUtilities.GenerateUniformRealValues(0.0, 16.0, 4)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(RayleighSum(nu, 1), 1.0 / (4.0 * (nu + 1.0)), 5.0E-4));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(RayleighSum(nu, 2), 1.0 / (16.0 * MoreMath.Sqr(nu + 1.0) * (nu + 2.0)), 1.0E-4));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(RayleighSum(nu, 3), 1.0 / (32.0 * MoreMath.Pow(nu + 1.0, 3) * (nu + 2.0) * (nu + 3.0)), 5.0E-5));
            }

        }

        private static double RayleighSum(double nu, int r) {
            double s = 0.0;
            for (int k = 1; k < 1024; k++) {
                double s_old = s;
                double ds = 1.0 / MoreMath.Pow(AdvancedMath.BesselJZero(nu, k), 2 * r);
                s += ds;
                if (ds < 1.0E-5 * s) {
                    return (s);
                }
            }
            throw new NonconvergenceException();
        }
    }
}
