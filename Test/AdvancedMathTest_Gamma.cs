using System;
using System.Linq;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using FluentAssertions;


using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public partial class AdvancedMathTest_Gamma {

        [TestMethod]
        public void GammaRecurrsion() {
            // DLMF 5.5.1
            // Limit x to avoid overflow.
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E2).Take(16)) {
                AdvancedMath.Gamma(x + 1.0).Should().BeNearly(x * AdvancedMath.Gamma(x));
            }
        }

        [TestMethod]
        public void GammaAtNegativeIntegers() {
            Assert.IsTrue(Double.IsInfinity(AdvancedMath.Gamma(0.0)));
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                Assert.IsTrue(Double.IsInfinity(AdvancedMath.Gamma(-n)));
            }
        }

        [TestMethod]
        public void GammaSpecialCases() {
            // It would be nice to be able to make more of these comparisons exact.
            Assert.IsTrue(Double.IsPositiveInfinity(AdvancedMath.Gamma(0.0)));
            AdvancedMath.Gamma(0.5).Should().BeNearly(Math.Sqrt(Math.PI));
            AdvancedMath.Gamma(1.0).Should().Be(1.0);
            AdvancedMath.Gamma(1.5).Should().BeNearly(Math.Sqrt(Math.PI) / 2.0);
            AdvancedMath.Gamma(2.0).Should().Be(1.0);
            AdvancedMath.Gamma(3.0).Should().Be(2.0);
            AdvancedMath.Gamma(4.0).Should().BeNearly(6.0);
            AdvancedMath.Gamma(Double.PositiveInfinity).Should().Be(Double.PositiveInfinity);
        }

        [TestMethod]
        public void GammaExtremeValues() {
            double eps = MoreMath.Pow(2.0, -1022); // Can't use 1.0 / Double.MaxValue because it's subnormal and looses accuracy
            AdvancedMath.Gamma(eps).Should().Be(1.0 / eps);
            AdvancedMath.Gamma(Double.MaxValue).Should().Be(Double.PositiveInfinity);
            AdvancedMath.Gamma(Double.NaN).Should().Be(Double.NaN);
        }

        [TestMethod]
        public void GammaReflection() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E2).Take(16)) {
                double GP = AdvancedMath.Gamma(x);
                double GN = AdvancedMath.Gamma(-x);
                (-x * GN * GP).Should().BeNearly(Math.PI / MoreMath.SinPi(x));
            }
        }

        [TestMethod]
        public void GammaDuplication() {
            // DLMF 5.5.5
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 85.0).Take(16)) {
                AdvancedMath.Gamma(2.0 * x).Should().BeNearly(
                    AdvancedMath.Gamma(x) * AdvancedMath.Gamma(x + 0.5) * Math.Pow(2.0, 2.0 * x - 1.0) / Math.Sqrt(Math.PI)
                );
            }
        }

        [TestMethod]
        public void GammaInequality() {
            foreach (double x in TestUtilities.GenerateRealValues(2.0, 1.0E2).Take(16)) {
                // for x >= 2
                double lower = Math.Pow(x / Math.E, x - 1.0);
                double upper = Math.Pow(x / 2.0, x - 1.0);
                double value = AdvancedMath.Gamma(x);
                Assert.IsTrue((lower <= value) && (value <= upper));
            }
        }

        [TestMethod]
        public void GammaAndRecriprocalInequality() {
            // DLMF 5.6.2
             foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(8)) {
                Assert.IsTrue(1.0 / AdvancedMath.Gamma(x) + 1.0 / AdvancedMath.Gamma(1.0 / x) <= 2.0);
            }
        }

        [TestMethod]
        public void GautschiInequality() {
            // DLMF 5.6.4
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(8)) {
                foreach (double s in TestUtilities.GenerateRealValues(1.0E-2, 1.0).Take(4)) {
                    double r = AdvancedMath.Gamma(x + 1.0) / AdvancedMath.Gamma(x + s);
                    r.Should().BeGreaterThan(Math.Pow(x, 1.0 - s));
                    r.Should().BeLessThan(Math.Pow(x + 1.0, 1.0 - s));
                }
            }
        }

        [TestMethod]
        public void GammaIntegral() {
            // DLMF 5.2.1
            foreach (double x in TestUtilities.GenerateRealValues(1.0, 1.0E2).Take(4)) {
                FunctionMath.Integrate(t => Math.Exp(-t) * Math.Pow(t, x - 1.0), 0.0, Double.PositiveInfinity).Value.Should().BeNearly(AdvancedMath.Gamma(x));
            }
        }

        [TestMethod]
        public void GammaTrottIdentity() {
            // http://mathworld.wolfram.com/GammaFunction.html
            double G1 = AdvancedMath.Gamma(1.0 / 24.0);
            double G5 = AdvancedMath.Gamma(5.0 / 24.0);
            double G7 = AdvancedMath.Gamma(7.0 / 24.0);
            double G11 = AdvancedMath.Gamma(11.0 / 24.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual((G1 * G11) / (G5 * G7), Math.Sqrt(3.0) * Math.Sqrt(2.0 + Math.Sqrt(3.0))));
        }

        private double GammaProduct(int r) {
            double p = 1.0;
            for (int i = 1; i < r; i++) {
                p *= AdvancedMath.Gamma(((double)i) / r);
            }
            return (p);
        }

        [TestMethod]
        public void GammaProducts() {
            // https://en.wikipedia.org/wiki/Particular_values_of_the_Gamma_function#Products
            // http://mathworld.wolfram.com/GammaFunction.html
            Assert.IsTrue(TestUtilities.IsNearlyEqual(GammaProduct(3), 2.0 * Math.PI / Math.Sqrt(3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(GammaProduct(4), Math.Sqrt(2.0 * Math.PI * Math.PI * Math.PI)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(GammaProduct(5), 4.0 * Math.PI * Math.PI / Math.Sqrt(5.0)));
        }

        [TestMethod]
        public void GammaMultiplication() {
            foreach (int k in new int[] { 2, 3, 4 }) {
                foreach (double z in TestUtilities.GenerateRealValues(1.0E-2, 10.0, 4)) {
                    double p = AdvancedMath.Gamma(z);
                    for (int i = 1; i < k; i++) {
                        p *= AdvancedMath.Gamma(z + ((double)i) / k);
                    }
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        p, Math.Pow(2.0 * Math.PI, (k - 1) / 2.0) * Math.Pow(k, 0.5 - k * z) * AdvancedMath.Gamma(k * z)
                    ));
                }
            }
        }


        [TestMethod]
        public void LogGammaSpecialValues() {
            AdvancedMath.LogGamma(0.0).Should().Be(Double.PositiveInfinity);
            AdvancedMath.LogGamma(0.5).Should().BeNearly(Math.Log(Math.Sqrt(Math.PI)));
            AdvancedMath.LogGamma(1.0).Should().Be(0.0);
            AdvancedMath.LogGamma(1.5).Should().BeNearly(Math.Log(Math.Sqrt(Math.PI) / 2.0));
            AdvancedMath.LogGamma(2.0).Should().Be(0.0);
            AdvancedMath.LogGamma(Double.PositiveInfinity).Should().Be(Double.PositiveInfinity);
        }

        [TestMethod]
        public void LogGammaExtremeValues() {
            double eps = MoreMath.Pow(2.0, -1022);
            AdvancedMath.LogGamma(eps).Should().Be(-Math.Log(eps)); // goes like -ln x for small x
            AdvancedMath.LogGamma(Double.MaxValue).Should().Be(Double.PositiveInfinity); // should overflow since goes like x log x
        }

        [TestMethod]
        public void LogGammaDuplication() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 24)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.LogGamma(2.0 * x),
                    AdvancedMath.LogGamma(x) + AdvancedMath.LogGamma(x + 0.5) + (2.0 * x - 0.5) * Math.Log(2.0) - 0.5 * Math.Log(2.0 * Math.PI)
                ));
            }
        }

        [TestMethod]
        public void LogGammaTriplication() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 24)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.LogGamma(3.0 * x),
                    AdvancedMath.LogGamma(x) + AdvancedMath.LogGamma(x + 1.0 / 3.0) + AdvancedMath.LogGamma(x + 2.0 / 3.0) + (3.0 * x - 0.5) * Math.Log(3.0) - Math.Log(2.0 * Math.PI)
                ));
            }
        }

        [TestMethod]
        public void LogGammaCosineIntegral () {
            // From Wolfram functions site
            FunctionMath.Integrate(t => AdvancedMath.LogGamma(t) * MoreMath.CosPi(2.0 * t), 0.0, 1.0).Value.Should().BeNearly(1.0 / 4.0);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void LogGammaNegativeArgument() {
            AdvancedMath.LogGamma(-0.5);
        }

        [TestMethod]
        public void GammaRatioInequality() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 16)) {
                double LB = Math.Log(Math.Sqrt(x + 0.25));
                double UB = Math.Log((x + 0.5) / Math.Sqrt(x + 0.75));
                double R = AdvancedMath.LogGamma(x + 1.0) - AdvancedMath.LogGamma(x + 0.5);
                Assert.IsTrue(LB <= R);
                Assert.IsTrue(R <= UB);
            }
        }


        [TestMethod]
        public void PsiSpecialCases() {
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(1.0), -AdvancedMath.EulerGamma));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(2.0), -AdvancedMath.EulerGamma + 1.0));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(1.0 / 2.0), -AdvancedMath.EulerGamma - 2.0 * Math.Log(2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(1.0 / 3.0), -AdvancedMath.EulerGamma - 3.0 * Math.Log(3.0) / 2.0 - Math.PI / 2.0 / Math.Sqrt(3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(1.0 / 4.0), -AdvancedMath.EulerGamma - 3.0 * Math.Log(2.0) - Math.PI / 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(1.0 / 6.0), -AdvancedMath.EulerGamma - 3.0 * Math.Log(3.0) / 2.0 - 2.0 * Math.Log(2.0) - Math.PI / 2.0 * Math.Sqrt(3.0)));

        }

        [TestMethod]
        public void PsiExtremeValues() {
            double eps = MoreMath.Pow(2.0, -1022);
            Assert.IsTrue(Double.IsInfinity(AdvancedMath.Psi(0.0)));
            AdvancedMath.Psi(eps).Should().Be(-1.0 / eps); // goes like -1/z near 0
            AdvancedMath.Psi(Double.MaxValue).Should().BeNearly(Math.Log(Double.MaxValue)); // goes like log(z) asymptotically
            Assert.IsTrue(AdvancedMath.Psi(Double.PositiveInfinity) == Double.PositiveInfinity);
            Assert.IsTrue(Double.IsNaN(AdvancedMath.Psi(Double.NaN)));
        }


        [TestMethod]
        public void PsiNegativeIntegers () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100).Take(4)) {
                Double.IsInfinity(AdvancedMath.Psi(-n)).Should().BeTrue();
            }
        }

        [TestMethod]
        public void PsiRecurrence() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(24)) {
                AdvancedMath.Psi(1.0 + x).Should().BeNearlySumOf(1.0 / x, AdvancedMath.Psi(x));
            }
        }

        [TestMethod]
        public void PsiReflection() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(16)) {
                AdvancedMath.Psi(1.0 - x).Should().BeNearlySumOf(Math.PI / MoreMath.TanPi(x), AdvancedMath.Psi(x));
            }
        }

        [TestMethod]
        public void PsiDuplication() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(24)) {
                AdvancedMath.Psi(2 * x).Should().BeNearly(AdvancedMath.Psi(x) / 2.0 + AdvancedMath.Psi(x + 0.5) / 2.0 + Math.Log(2.0));
            }
        }

        [TestMethod]
        public void PsiIntegration () {
            // Integral of Psi is LogGamma
            double[] points = TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(4).ToArray();
            foreach (double x in points) {
                foreach (double y in points) {
                    FunctionMath.Integrate(t => AdvancedMath.Psi(t), x, y).Value.Should().BeNearly(AdvancedMath.LogGamma(y) - AdvancedMath.LogGamma(x));
                }
            }
        }

        [TestMethod]
        public void PsiInequality() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(16)) {
                AdvancedMath.Psi(x).Should().BeInRange(Math.Log(x) - 1.0 / x, Math.Log(x) - 0.5 / x);
            }
        }

        [TestMethod]
        public void PsiSums () {
            // https://plouffe.fr/articles/The%20many%20faces%20of%20the%20polygamma%20function%202016.pdf
            (AdvancedMath.Psi(1.0 / 4.0) - AdvancedMath.Psi(3.0 / 4.0)).Should().BeNearly(-Math.PI);
            (AdvancedMath.Psi(1.0 / 5.0) - AdvancedMath.Psi(2.0 / 5.0) - AdvancedMath.Psi(3.0 / 5.0) + AdvancedMath.Psi(4.0 / 5.0)).Should().BeNearly(-2.0 * Math.Log(AdvancedMath.GoldenRatio) * Math.Sqrt(5.0));
        }

    }
}
