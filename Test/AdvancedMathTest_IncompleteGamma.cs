using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using FluentAssertions;
using System.Linq;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_IncompleteGamma {

        // Regularized incomplete Gamma

        [TestMethod]
        public void RegularizedIncompleteGammaLimits () {
            foreach (double a in TestUtilities.GenerateRealValues(1.0E-3, 1.0E3).Take(8)) {
                AdvancedMath.LeftRegularizedGamma(a, 0.0).Should().Be(0.0);
                AdvancedMath.LeftRegularizedGamma(a, Double.PositiveInfinity).Should().Be(1.0);
                AdvancedMath.RightRegularizedGamma(a, 0.0).Should().Be(1.0);
                AdvancedMath.RightRegularizedGamma(a, Double.PositiveInfinity).Should().Be(0.0);
            }

        }

        [TestMethod]
        public void RegularizedIncompleteGammaZeroShape () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-3, 1.0E3).Take(8)) {
                AdvancedMath.LeftRegularizedGamma(0.0, x).Should().Be(1.0);
                AdvancedMath.RightRegularizedGamma(0.0, x).Should().Be(0.0);
            }
        }

        [TestMethod]
        public void RegularizedIncompleteGammaExtremeValues () {
            foreach (double a in TestUtilities.GenerateRealValues(1.0E-2, 1.0E1).Take(4)) {
                AdvancedMath.LeftRegularizedGamma(a, TestUtilities.SmallestNormal).Should().BeNearly(Math.Pow(TestUtilities.SmallestNormal, a) / AdvancedMath.Gamma(a + 1.0));
            }
        }

        [TestMethod]
        public void RegularizedIncompleteGammaIntegerInequality() {
            // DLMF 8.10.13
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000).Take(8)) {
                AdvancedMath.RightRegularizedGamma(n, n).Should().BeLessThan(0.5);
                AdvancedMath.RightRegularizedGamma(n, n - 1).Should().BeGreaterThan(0.5);
            }
        }

        [TestMethod]
        public void RegularizedIncompleteGammaRecurrence() {
            // DLMF 8.8.5-8.8.6
            foreach (double a in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(5)) {
                    AdvancedMath.LeftRegularizedGamma(a, x).Should().BeNearly(AdvancedMath.LeftRegularizedGamma(a + 1.0, x) + Math.Pow(x, a) * Math.Exp(-x) / AdvancedMath.Gamma(a + 1.0));
                    AdvancedMath.RightRegularizedGamma(a + 1.0, x).Should().BeNearly(AdvancedMath.RightRegularizedGamma(a, x) + Math.Pow(x, a) * Math.Exp(-x) / AdvancedMath.Gamma(a + 1.0)); 
                }
            }
        }

        [TestMethod]
        public void RegularizedIncompleteGammaExponential() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.RightRegularizedGamma(1.0, x), Math.Exp(-x)
                ));
            }
        }

        [TestMethod]
        public void RegularizedIncompleteGammaUnitarity() {
            // P + Q = 1
            foreach (double a in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(6)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(6)) {
                    (AdvancedMath.LeftRegularizedGamma(a, x) + AdvancedMath.RightRegularizedGamma(a, x)).Should().BeNearly(1.0);
                }
            }
        }

        [TestMethod]
        public void RegularizedIncompleteGammaAgreement () {
            foreach (double a in TestUtilities.GenerateRealValues(1.0E-4, 1.0E2).Take(6)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4).Take(6)) {
                    AdvancedMath.LowerIncompleteGamma(a, x).Should().BeNearly(AdvancedMath.LeftRegularizedGamma(a, x) * AdvancedMath.Gamma(a));
                    AdvancedMath.UpperIncompleteGamma(a, x).Should().BeNearly(AdvancedMath.RightRegularizedGamma(a, x) * AdvancedMath.Gamma(a));
                }
            }
        }

        [TestMethod]
        public void RegularizedIncompleteGammaExpIntegral () {
            // DLMF 8.14.1
            foreach (double b in TestUtilities.GenerateRealValues(1.0E-1, 10E1).Take(3)) {
                IntegrationSettings s = new IntegrationSettings() { AbsolutePrecision = 0.0 };
                FunctionMath.Integrate(x => Math.Exp(-x) * AdvancedMath.LeftRegularizedGamma(b, x), 0.0, Double.PositiveInfinity, s).Value.Should().BeNearly(Math.Pow(2.0, -b));
                FunctionMath.Integrate(x => Math.Exp(-x) * AdvancedMath.RightRegularizedGamma(b, x), 0.0, Double.PositiveInfinity, s).Value.Should().BeNearly(1.0 - Math.Pow(2.0, -b));
            }
        }

        [TestMethod]
        public void LowerIncompleteGammaExtremeValues () {
            foreach (double a in TestUtilities.GenerateRealValues(1.0E-2, 1.0).Take(4)) {
                AdvancedMath.LowerIncompleteGamma(a, TestUtilities.SmallestNormal).Should().BeNearly(Math.Pow(TestUtilities.SmallestNormal, a) / a);
            }
        }

        [TestMethod]
        public void IncompleteGammaLimits () {
            foreach (double a in TestUtilities.GenerateRealValues(1.0E-3, 1.0E-2).Take(8)) {
                AdvancedMath.LowerIncompleteGamma(a, 0.0).Should().Be(0.0);
                AdvancedMath.LowerIncompleteGamma(a, Double.PositiveInfinity).Should().BeNearly(AdvancedMath.Gamma(a));
                AdvancedMath.UpperIncompleteGamma(a, 0.0).Should().BeNearly(AdvancedMath.Gamma(a));
                AdvancedMath.UpperIncompleteGamma(a, Double.PositiveInfinity).Should().Be(0.0);
            }
        }

        [TestMethod]
        public void IncompleteGammaInequality() {
            // DLMF 8.10.1
            foreach (double a in TestUtilities.GenerateRealValues(1.0E-4, 1).Take(5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E1).Take(6)) {
                    AdvancedMath.UpperIncompleteGamma(a, x).Should().BeLessOrEqualTo(Math.Exp(-x) * Math.Pow(x, a - 1.0));
                }
            }
        }

        [TestMethod]
        public void IncompleteGammaUnitarity () {
            // DLMF 8.2.3
            Random rng = new Random(3);
            foreach (double a in TestUtilities.GenerateRealValues(1.0E-2, 1.5E2, rng).Take(4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E3, rng).Take(5)) {
                    (AdvancedMath.LowerIncompleteGamma(a, x) + AdvancedMath.UpperIncompleteGamma(a, x)).Should().BeNearly(AdvancedMath.Gamma(a));
                }
            }

        }

        [TestMethod]
        public void IncompleteGammaRecurrence () {
            // DLMF 8.8.1 and 8.8.2
            Random rng = new Random(2);
            foreach (double a in TestUtilities.GenerateRealValues(1.0E-2, 1.5E2, rng).Take(4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, rng).Take(4)) {
                    (a * AdvancedMath.LowerIncompleteGamma(a, x)).Should().BeNearly(AdvancedMath.LowerIncompleteGamma(a + 1.0, x) + Math.Pow(x, a) * Math.Exp(-x));
                    AdvancedMath.UpperIncompleteGamma(a + 1.0, x).Should().BeNearly(a * AdvancedMath.UpperIncompleteGamma(a, x) + Math.Pow(x, a) * Math.Exp(-x));
                }
            }
        }

        [TestMethod]
        public void UpperIncompleteGammaIntegralE1Agreement() {
            // DLMF 8.4.4
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(16)) {
                AdvancedMath.UpperIncompleteGamma(0.0, x).Should().BeNearly(AdvancedMath.IntegralE(1, x));
            }
        }

        [TestMethod]
        public void UpperIncompleteGammaErfcAgreement() {
            // DLMF 8.4.6
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(16)) {
                AdvancedMath.UpperIncompleteGamma(0.5, x * x).Should().BeNearly(Math.Sqrt(Math.PI) * AdvancedMath.Erfc(x));
            }
        }

        [TestMethod]
        public void LowerIncompleteGammaErfAgreement() {
            // DLMF 8.4.6
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(16)) {
                AdvancedMath.LowerIncompleteGamma(0.5, x * x).Should().BeNearly(Math.Sqrt(Math.PI) * AdvancedMath.Erf(x));
            }
        }

        [TestMethod]
        public void UpperIncompleteGammaExpAgreement() {
            // DLMF 8.4.5
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(16)) {
                AdvancedMath.UpperIncompleteGamma(1.0, x).Should().BeNearly(Math.Exp(-x));
            }
        }

        [TestMethod]
        public void LowerIncompleteGammaExpMinusOneAgreement () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(16)) {
                AdvancedMath.LowerIncompleteGamma(1.0, x).Should().BeNearly(-MoreMath.ExpMinusOne(-x));
            }
        }

        [TestMethod]
        public void IncompleteGammaIntegralDefinition () {
            // DLMF 8.2.2
            Random rng = new Random(1);
            IntegrationSettings s = new IntegrationSettings() { EvaluationBudget = 10000, AbsolutePrecision = 0.0 };
            foreach (double a in TestUtilities.GenerateRealValues(0.5, 1.5E2, rng).Take(4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, rng).Take(4)) {
                    AdvancedMath.LowerIncompleteGamma(a, x).Should().BeNearly(
                        FunctionMath.Integrate((double t) => Math.Pow(t, a - 1.0) * Math.Exp(-t), 0.0, x, s)
                    );
                    AdvancedMath.UpperIncompleteGamma(a, x).Should().BeNearly(
                        FunctionMath.Integrate((double t) => Math.Pow(t, a - 1.0) * Math.Exp(-t), x, Double.PositiveInfinity, s)
                    );
                }
            }
        }

        [TestMethod]
        public void UpperIncompleteGammaExpIntegral () {
            // DLMF 8.14.2
            foreach (double b in TestUtilities.GenerateRealValues(1.0E-1, 1.0E1).Take(3)) {
                FunctionMath.Integrate((double x) => Math.Exp(-x) * AdvancedMath.UpperIncompleteGamma(b, x), 0.0, Double.PositiveInfinity).Value.Should().BeNearly(AdvancedMath.Gamma(b) * (1.0 - Math.Pow(2.0, -b)));
            }
        }

        [TestMethod]
        public void LowerIncompleteGammaPowerIntegral () {
            // DLMF 8.14.3
            // Since \gamma ~ x^a, need a > 2 to avoid accuracy-killing power-law singularity at origin
            foreach (double b in TestUtilities.GenerateRealValues(2.0, 1.5E2).Take(3)) {
                FunctionMath.Integrate((double x) => AdvancedMath.LowerIncompleteGamma(b, x) / (x * x), 0.0, Double.PositiveInfinity).Value.Should().BeNearly(AdvancedMath.Gamma(b - 1.0));
            }
        }

        [TestMethod]
        public void UpperIncompleteGammaIntegral () {
            // DLMF 8.14.4
            foreach (double b in TestUtilities.GenerateRealValues(1.0E-1, 1.0E1).Take(4)) {
                FunctionMath.Integrate((double x) => AdvancedMath.UpperIncompleteGamma(b, x), 0.0, Double.PositiveInfinity).Value.Should().BeNearly(AdvancedMath.Gamma(b + 1.0));
            }
        }


        [TestMethod]
        public void UpperIncompleteGammaSepcialValues () {
            //AdvancedMath.UpperIncompleteGamma(0.0, 0.0).Should().Be(Double.PositiveInfinity);
            AdvancedMath.UpperIncompleteGamma(0.0, Double.PositiveInfinity).Should().Be(0.0);
            AdvancedMath.UpperIncompleteGamma(0.5, 0.0).Should().BeNearly(Math.Sqrt(Math.PI));
            AdvancedMath.UpperIncompleteGamma(0.5, Double.PositiveInfinity).Should().Be(0.0);
            AdvancedMath.UpperIncompleteGamma(1.0, 0.0).Should().Be(1.0);
            AdvancedMath.UpperIncompleteGamma(1.0, 1.0).Should().BeNearly(1.0 / Math.E);
            AdvancedMath.UpperIncompleteGamma(2.0, 0.0).Should().Be(1.0);
            AdvancedMath.UpperIncompleteGamma(2.0, 1.0).Should().BeNearly(2.0 / Math.E);
        }

    }
}
