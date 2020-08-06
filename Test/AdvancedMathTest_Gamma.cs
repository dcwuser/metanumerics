using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public partial class AdvancedMathTest_Gamma {

        [TestMethod]
        public void GammaRecurrsion () {
            // Limit x to avoid overflow.
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E2, 16)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(x + 1.0), x * AdvancedMath.Gamma(x)));
            }
        }

        [TestMethod]
        public void GammaAtNegativeIntegers () {
            Assert.IsTrue(Double.IsInfinity(AdvancedMath.Gamma(0.0)));
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                Assert.IsTrue(Double.IsInfinity(AdvancedMath.Gamma(-n)));
            }
        }

        [TestMethod]
        public void GammaSpecialCases () {
            // It would be nice to be able to make more of these comparisons exact.
            Assert.IsTrue(Double.IsPositiveInfinity(AdvancedMath.Gamma(0.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(0.5), Math.Sqrt(Math.PI)));
            Assert.IsTrue(AdvancedMath.Gamma(1.0) == 1.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(1.5), Math.Sqrt(Math.PI) / 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(2.0), 1.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(3.0), 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(4.0), 6.0));
            Assert.IsTrue(Double.IsPositiveInfinity(AdvancedMath.Gamma(Double.PositiveInfinity)));
        }

        [TestMethod]
        public void GammaExtremeValues () {
            Assert.IsTrue(Double.IsPositiveInfinity(AdvancedMath.Gamma(Double.MaxValue)));
            Assert.IsTrue(Double.IsNaN(AdvancedMath.Gamma(Double.NaN)));
        }

        [TestMethod]
        public void GammaReflection () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E2, 16)) {
                double GP = AdvancedMath.Gamma(x);
                double GN = AdvancedMath.Gamma(-x);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(-x * GN * GP, Math.PI / MoreMath.SinPi(x)));
            }
        }

        [TestMethod]
        public void GammaDuplication () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E2, 16)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.Gamma(2.0 * x),
                    AdvancedMath.Gamma(x) * AdvancedMath.Gamma(x + 0.5) * Math.Pow(2.0, 2.0 * x - 1.0) / Math.Sqrt(Math.PI)
                ));
            }
        }

        [TestMethod]
        public void GammaInequality () {
            foreach (double x in TestUtilities.GenerateRealValues(2.0, 1.0E2, 16)) {
                // for x >= 2
                double lower = Math.Pow(x / Math.E, x - 1.0);
                double upper = Math.Pow(x / 2.0, x - 1.0);
                double value = AdvancedMath.Gamma(x);
                Assert.IsTrue((lower <= value) && (value <= upper));
            }
        }

        [TestMethod]
        public void GammaIntegral () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0, 1.0E2, 4)) {
                Func<double, double> f = delegate (double t) {
                    return (Math.Pow(t, x - 1.0) * Math.Exp(-t));
                };
                Interval r = Interval.FromEndpoints(0.0, Double.PositiveInfinity);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(x), FunctionMath.Integrate(f, r)));
            }
        }

        [TestMethod]
        public void GammaTrottIdentity () {
            // http://mathworld.wolfram.com/GammaFunction.html
            double G1 = AdvancedMath.Gamma(1.0 / 24.0);
            double G5 = AdvancedMath.Gamma(5.0 / 24.0);
            double G7 = AdvancedMath.Gamma(7.0 / 24.0);
            double G11 = AdvancedMath.Gamma(11.0 / 24.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual((G1 * G11) / (G5 * G7), Math.Sqrt(3.0) * Math.Sqrt(2.0 + Math.Sqrt(3.0))));
        }

        private double GammaProduct (int r) {
            double p = 1.0;
            for (int i = 1; i < r; i++) {
                p *= AdvancedMath.Gamma(((double) i) / r);
            }
            return (p);
        }

        [TestMethod]
        public void GammaProducts () {
            // https://en.wikipedia.org/wiki/Particular_values_of_the_Gamma_function#Products
            // http://mathworld.wolfram.com/GammaFunction.html
            Assert.IsTrue(TestUtilities.IsNearlyEqual(GammaProduct(3), 2.0 * Math.PI / Math.Sqrt(3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(GammaProduct(4), Math.Sqrt(2.0 * Math.PI * Math.PI * Math.PI)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(GammaProduct(5), 4.0 * Math.PI * Math.PI / Math.Sqrt(5.0)));
        }

        [TestMethod]
        public void GammaMultiplication () {
            foreach (int k in new int[] { 2, 3, 4 }) {
                foreach (double z in TestUtilities.GenerateRealValues(1.0E-2, 10.0, 4)) {
                    double p = AdvancedMath.Gamma(z);
                    for (int i = 1; i < k; i++) {
                        p *= AdvancedMath.Gamma(z + ((double) i) / k);
                    }
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        p, Math.Pow(2.0 * Math.PI, (k - 1) / 2.0) * Math.Pow(k, 0.5 - k * z) * AdvancedMath.Gamma(k * z)
                    ));
                }
            }
        }


        [TestMethod]
        public void LogGammaSpecialValues () {
            Assert.IsTrue(Double.IsPositiveInfinity(AdvancedMath.LogGamma(0.0)));
            Assert.IsTrue(AdvancedMath.LogGamma(1.0) == 0.0);
            Assert.IsTrue(AdvancedMath.LogGamma(2.0) == 0.0);
            Assert.IsTrue(Double.IsPositiveInfinity(AdvancedMath.LogGamma(Double.PositiveInfinity)));
        }

        [TestMethod]
        public void LogGammaExtremeValues () {
            //Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.LogGamma(1.0 / Double.MaxValue), Math.Log(Double.MaxValue)));
            double z = AdvancedMath.LogGamma(Double.MaxValue);
            //double z = double.PositiveInfinity;
            //Console.WriteLine(z);
            Console.WriteLine(Double.IsPositiveInfinity(z));
            Console.WriteLine(AdvancedMath.LogGamma(Double.MaxValue));
            Console.WriteLine(Double.IsPositiveInfinity(AdvancedMath.LogGamma(Double.MaxValue)));
            Assert.IsTrue(Double.IsPositiveInfinity(AdvancedMath.LogGamma(Double.MaxValue)));
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
