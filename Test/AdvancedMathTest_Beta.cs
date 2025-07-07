using System;
using System.Linq;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using FluentAssertions;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using System.Data;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_Beta {

        [TestMethod]
        public void BetaSpecialValues() {

            AdvancedMath.Beta(0.5, 0.5).Should().BeNearly(Math.PI);
            AdvancedMath.Beta(1.0, 1.0).Should().Be(1.0);
            AdvancedMath.Beta(1.0, 2.0).Should().BeNearly(1.0 / 2.0);
            AdvancedMath.Beta(2.0, 2.0).Should().BeNearly(1.0 / 6.0);

        }

        [TestMethod]
        public void BetaArgumentZero () {
            foreach (double y in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(4)) {
                Double.IsInfinity(AdvancedMath.Beta(0.0, y)).Should().BeTrue();
            }
        }

        [TestMethod]
        public void BetaArgumentOne () {
            foreach (double y in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(4)) {
                AdvancedMath.Beta(1.0, y).Should().BeNearly(1.0 / y);
            }
        }

        [TestMethod]
        public void BetaSymmetry () {
            Random rng = new Random(2);
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, rng).Take(4)) {
                foreach (double y in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, rng).Take(4)) {
                    AdvancedMath.Beta(x, y).Should().Be(AdvancedMath.Beta(y, x));
                }
            }
        }

        [TestMethod]
        public void BetaIntegral() {

            // B(a, b) = \int_0^1 \! dt \, t^{a-1} (1 - t)^{b-1}

            // a and b values less than one correspond to singular integrals,
            // which are numerically difficult.

            // The range of these integrals varies enoromously, and there is no cancelation,
            // so specify a pure relative accuracy target.
            IntegrationSettings settings = new IntegrationSettings() {
                AbsolutePrecision = 0.0
            };

            Random rng = new Random(3);
            foreach (double a in TestUtilities.GenerateRealValues(1.0, 100.0, rng).Take(4)) {
                foreach (double b in TestUtilities.GenerateRealValues(1.0, 100.0, rng).Take(4)) {
                    FunctionMath.Integrate(
                        t => Math.Pow(t, a - 1.0) * Math.Pow(1.0 - t, b - 1.0), 0.0, 1.0, settings
                    ).Value.Should().BeNearly(AdvancedMath.Beta(a, b));
                }
            }
        }


        [TestMethod]
        public void BetaSumRecurrence() {
            // B(a, b) = B(a + 1, b) + B(a, b + 1)
            Random rng = new Random(4);
            foreach (double a in TestUtilities.GenerateRealValues(1.0E-2, 1.0E3, rng).Take(4)) {
                foreach (double b in TestUtilities.GenerateRealValues(1.0E-3, 1.0E2, rng).Take(4)) {
                    AdvancedMath.Beta(a, b).Should().BeNearly(AdvancedMath.Beta(a + 1.0, b) + AdvancedMath.Beta(a, b + 1.0));
                }
            }
        }

        [TestMethod]
        public void BetaProductRecurrence() {
            // B(a, b + 1) = b / (a + b) B(a, b)
            Random rng = new Random(5);
            foreach (double a in TestUtilities.GenerateRealValues(1.0E-3, 1.0E2, rng).Take(4)) {
                foreach (double b in TestUtilities.GenerateRealValues(1.0E-2, 1.0E3, rng).Take(4)) {
                    AdvancedMath.Beta(a, b + 1.0).Should().BeNearly(b / (a + b) * AdvancedMath.Beta(a, b));
                }
            }
        }

        [TestMethod]
        public void BetaReflection() {
            foreach (double a in TestUtilities.GenerateRealValues(0.1, 10.0, 3)) {
                foreach (double b in TestUtilities.GenerateUniformRealValues(0.0, 1.0, 3)) {
                    double B1 = AdvancedMath.Beta(a, b);
                    double B2 = AdvancedMath.Beta(a + b, 1.0 - b);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(B1 * B2, Math.PI / a / Math.Sin(Math.PI * b)));
                }
            }
        }

 
        [TestMethod]
        public void BetaTripleRelation () {
            // B(p, q) B(p + q, r) = B(q, r) B(q + r, p)
            Random rng = new Random(6);
            foreach (double p in TestUtilities.GenerateRealValues(1.0E-3, 1.0E1, rng).Take(4)) {
                foreach (double q in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, rng).Take(4)) {
                    foreach (double r in TestUtilities.GenerateRealValues(1.0E-1, 1.0E3, rng).Take(4)) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(
                            AdvancedMath.Beta(p, q) * AdvancedMath.Beta(p + q, r),
                            AdvancedMath.Beta(q, r) * AdvancedMath.Beta(q + r, p)
                        ));
                        //(AdvancedMath.Beta(p, q) * AdvancedMath.Beta(p + q, r)).Should().BeNearly(
                        //    AdvancedMath.Beta(q, r) * AdvancedMath.Beta(q + r, p), new EvaluationSettings() { RelativePrecision = 2.0E-14, AbsolutePrecision = 0.0 }
                        //);
                    }
                }
            }
        }

        [TestMethod]
        public void BetaMappings () {
            for (double x = -1.0; x <= 0.0; x += 0.5) {
                for (double y = -1.0; y <= 1.0; y += 0.5) {
                    Console.WriteLine($"({x},{y}) => ({1 - x},{x + y}) and ({1 - x - y},{y})");
                }
            }
        }


        [TestMethod]
        public void BetaInequality() {
            Random rng = new Random(7);
            foreach (double a in TestUtilities.GenerateRealValues(10.0E-2, 10.0E3, rng).Take(8)) {
                foreach (double b in TestUtilities.GenerateRealValues(10.0E-2, 10.0E3, rng).Take(8)) {
                    AdvancedMath.Beta(a, b).Should().BeGreaterOrEqualTo(Math.Sqrt(AdvancedMath.Beta(a, a) * AdvancedMath.Beta(b, b)));
                }
            }
        }

        [TestMethod]
        public void BetaInequalityOnUnitSquare() {

            // The most famous inequality on the unit square is Dragomir's B(a, b) \le \frac{1}{ab}.
            // It can be improved and bounded on both sides via Ivady's
            //   a + b - ab <= a b B(a, b) <= \frac{a + b}{1 + ab}
            // see http://nntdm.net/papers/nntdm-21/NNTDM-21-2-01-07.pdf

            foreach (double a in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 4)) {
                foreach (double b in TestUtilities.GenerateRealValues(1.0E-3, 1.0, 4)) {

                    double abB = a * b * AdvancedMath.Beta(a, b);

                    Assert.IsTrue((a + b - a * b) <= abB);
                    Assert.IsTrue(abB <= (a + b) / (1.0 + a * b));

                }
            }

        }

        [TestMethod]
        public void BetaGammaAgreement () {
            Random rng = new Random(4);
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, rng).Take(4)) {
                foreach (double y in TestUtilities.GenerateRealValues(1.0E-2, 1.0E1, rng).Take(4)) {
                    AdvancedMath.Beta(x, y).Should().BeNearly(AdvancedMath.Gamma(x) / AdvancedMath.Gamma(x + y) * AdvancedMath.Gamma(y));
                }
            }
        }

        [TestMethod]
        public void LogBetaSpecialCases() {
            AdvancedMath.LogBeta(1.0, 1.0).Should().Be(0.0);
        }

        [TestMethod]
        public void LogBetaExtremeValues() {
            Assert.IsTrue(AdvancedMath.LogBeta(0.0, 0.0) == Double.PositiveInfinity);
            Assert.IsTrue(AdvancedMath.LogBeta(1.0E-2, 0.0) == Double.PositiveInfinity);
            Assert.IsTrue(AdvancedMath.LogBeta(0.0, 1.0E2) == Double.PositiveInfinity);
            Assert.IsTrue(Double.IsNaN(AdvancedMath.Beta(1.0E2, Double.NaN)));
            Assert.IsTrue(Double.IsNaN(AdvancedMath.Beta(Double.NaN, 1.0E-2)));
        }

        [TestMethod]
        public void LogBetaAgreement() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 4)) {
                foreach (double y in TestUtilities.GenerateRealValues(1.0E-3, 1.0E1, 4)) {
                    AdvancedMath.LogBeta(x, y).Should().BeNearly(Math.Log(AdvancedMath.Beta(x, y)));
                }
            }
        }

    }
}
