using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_Carlson {

        [TestMethod]
        public void CarlsonNormalization () {

            // All Carlson integrals are normalized so that they are 1 when all arguments are 1.

            Assert.IsTrue(AdvancedMath.CarlsonF(1.0, 1.0, 1.0) == 1.0);
            Assert.IsTrue(AdvancedMath.CarlsonG(1.0, 1.0, 1.0) == 1.0);
            Assert.IsTrue(AdvancedMath.CarlsonD(1.0, 1.0, 1.0) == 1.0);

        }

        [TestMethod]
        public void CarlsonScaling () {

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-3, 1.0E1, 2)) {
                foreach (double y in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 2)) {
                    foreach (double z in TestUtilities.GenerateRealValues(1.0E-1, 1.0E3, 2)) {
                        foreach (double lambda in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 4)) {

                            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                                AdvancedMath.CarlsonF(lambda * x, lambda * y, lambda * z),
                                AdvancedMath.CarlsonF(x, y, z) / Math.Sqrt(lambda)    
                            ));

                            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                                AdvancedMath.CarlsonG(lambda * x, lambda * y, lambda * z),
                                AdvancedMath.CarlsonG(x, y, z) * Math.Sqrt(lambda)
                            ));

                            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                                AdvancedMath.CarlsonD(lambda * x, lambda * y, lambda * z),
                                AdvancedMath.CarlsonD(x, y, z) * Math.Pow(lambda, -1.5)
                            ));

                        }
                    }
                }
            }

        }

        [TestMethod]
        public void CarlsonFSpecialCases () {

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 8)) {

                // DLMF 19.20.1
                // R_F(x, x, x) = 1 / \sqrt{x}
                // R_F(0, x, x) = \frac{\pi}{2} \frac{1}{\sqrt{x}}
                // R_F(0, 0, x) = \infty


                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.CarlsonF(x, x, x),
                    1.0 / Math.Sqrt(x)
                ));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.CarlsonF(0, x, x),
                    Math.PI / 2.0 / Math.Sqrt(x)
                ));

            }

        }

        [TestMethod]
        public void CarlsonGSpecialCases () {

            // R_G(x, x, x) = \sqrt{x}
            // R_G(0, x, x) = \pi / 4 \sqrt{x}
            // R_G(0, 0, x) = 1/2 \sqrt{x}
            // R_G(0, 0, 0) = 0

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 8)) {

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.CarlsonG(x, x, x),
                    Math.Sqrt(x)
                ));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.CarlsonG(0.0, x, x),
                    Math.PI / 4.0 * Math.Sqrt(x)
                ));

                /*
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.CarlsonG(0.0, 0.0, x),
                    1.0 / 2.0 * Math.Sqrt(x)
                ));
                */

            }

        }

        [TestMethod]
        public void CarslonDSpecialCases () {

            // DLMF 19.20.18
            // R_D(x, x, x) = x^{-3/2}
            // R_D(0, x, x) = 3/4 \pi x^{-3/2}
            // R_D(0, 0, x) = \infty

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 8)) {

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.CarlsonD(x, x, x), Math.Pow(x, -3.0 / 2.0)
                ));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.CarlsonD(0, x, x), 3.0 * Math.PI / 4.0 * Math.Pow(x, -3.0 / 2.0)
                ));

            }

        }

        [TestMethod]
        public void CarlsonFInequality () {

            // DLMF 19.24.10
            // \frac{3}{\sqrt{x} + \sqrt{y} + \sqrt{z}} \le R_F(x, y, z) \le \frac{1}{(x y z)^{1/6}}

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-3, 1.0E3, 4)) {
                foreach (double y in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 4)) {
                    foreach (double z in TestUtilities.GenerateRealValues(1.0E-1, 1.0E5, 4)) {

                        double min = 3.0 / (Math.Sqrt(x) + Math.Sqrt(y) + Math.Sqrt(z));
                        double R_F = AdvancedMath.CarlsonF(x, y, z);
                        double max = 1.0 / Math.Pow(x * y * z, 1.0 / 6.0);

                        Assert.IsTrue(min <= R_F && R_F <= max);

                    }
                }
            }

        }

        [TestMethod]
        public void CarlsonGInequality () {

            // DLMF 19.24.12

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-3, 1.0E3, 4)) {
                foreach (double y in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 4)) {
                    foreach (double z in TestUtilities.GenerateRealValues(1.0E-1, 1.0E5, 4)) {

                        double min = (Math.Sqrt(x) + Math.Sqrt(y) + Math.Sqrt(z)) / 3.0;
                        double R_G = AdvancedMath.CarlsonG(x, y, z);
                        double max = Math.Min(
                            Math.Sqrt((x + y + z) / 3.0),
                            (x * x + y * y + z * z) / Math.Sqrt(x * y * z) / 3.0
                        );

                        Assert.IsTrue(min <= R_G && R_G <= max);

                    }
                }
            }

        }

    }

}
