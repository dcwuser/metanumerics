using System;
using System.Linq;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using FluentAssertions;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_Pochhammer {

        [TestMethod]
        public void LogPochhamerSepcialValues () {
            AdvancedMath.LogPochhammer(0.0, 0.0).Should().Be(0.0); // because (0)_0 = 1
            AdvancedMath.LogPochhammer(0.0, 1.0).Should().Be(Double.NegativeInfinity); // (0)_1 = 0
            AdvancedMath.LogPochhammer(1.0, 0.0).Should().Be(0.0); // (1)_0 = 1
            AdvancedMath.LogPochhammer(1.0, 1.0).Should().Be(0.0); // (1)_1 = 1
        }

        [TestMethod]
        public void LogPochhammerExtremeValues () {
            double eps = Math.Pow(1.0, -1022);
            AdvancedMath.LogPochhammer(0.0, Double.PositiveInfinity).Should().Be(Double.NegativeInfinity);
            AdvancedMath.LogPochhammer(1.0, Double.PositiveInfinity).Should().Be(Double.PositiveInfinity);
            AdvancedMath.LogPochhammer(0.0, Double.MaxValue).Should().Be(Double.NegativeInfinity);
            AdvancedMath.LogPochhammer(1.0, Double.MaxValue).Should().Be(Double.PositiveInfinity); // Stirling says \ln\Gamma(z) > z \ln z for large z
        }

        [TestMethod]
        public void LogPochhammerZeroExponent () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(8)) {
                AdvancedMath.LogPochhammer(x, 0.0).Should().Be(0.0);
            }
        }

        [TestMethod]
        public void LogPochhammerOneExponent() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(8)) {
                AdvancedMath.LogPochhammer(x, 1.0).Should().Be(Math.Log(x));
            }
        }

        [TestMethod]
        public void LogPochhammerZeroBase () {
            foreach (double y in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(8)) {
                AdvancedMath.LogPochhammer(0.0, y).Should().Be(Double.NegativeInfinity);
            }
        }

        [TestMethod]
        public void LogPochhammerLogGammaAgreement () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(8)) {
                foreach (double y in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4).Take(8)) {
                    AdvancedMath.LogGamma(x + y).Should().BeNearlySumOf(AdvancedMath.LogGamma(x) + AdvancedMath.LogPochhammer(x, y));
                }
            }
        }

        [TestMethod]
        public void LogPochhammerPochhammerAgreement () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(4)) {
                foreach (double y in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(4)) {
                    Math.Log(AdvancedMath.Pochhammer(x, y)).Should().BeNearly(AdvancedMath.LogPochhammer(x, y));
                }
            }
        }

        [TestMethod]
        public void PochhammerExtremeValues () {
            double eps = Math.Pow(2.0, -1022);
            AdvancedMath.Pochhammer(eps, eps).Should().Be(0.5);
            AdvancedMath.Pochhammer(eps, Double.MaxValue).Should().Be(Double.PositiveInfinity);
            AdvancedMath.Pochhammer(eps, Double.PositiveInfinity).Should().Be(Double.PositiveInfinity);
            AdvancedMath.Pochhammer(Double.MaxValue, eps).Should().Be(1.0);
            AdvancedMath.Pochhammer(Double.MaxValue, Double.MaxValue).Should().Be(Double.PositiveInfinity);
            AdvancedMath.Pochhammer(Double.MaxValue, Double.PositiveInfinity).Should().Be(Double.PositiveInfinity); // NaN
        }

        // Note Pochhammer(0,0) is ill-defined.
        // Should be zero because Pochhammer(0,y) = 0.
        // Should be one because Pochhammer(x,0) = 1.
        // Mathematica defines it as one, presumably by analogy with Pow(x, y) behavior.

        [TestMethod]
        public void PochhammerSpecialIntegerArguments () {
            // Arguments for which numerator or denominator Gamma vanishes.
            // This tests whether we are taking ratio correctly.
            AdvancedMath.Pochhammer(0, 0).Should().Be(1.0);
            AdvancedMath.Pochhammer(0, 1).Should().Be(0.0);
            AdvancedMath.Pochhammer(0, 2).Should().Be(0.0);
            AdvancedMath.Pochhammer(1, 0).Should().Be(1.0);
            AdvancedMath.Pochhammer(1, 1).Should().Be(1.0);
            AdvancedMath.Pochhammer(1, 2).Should().Be(2.0);
            AdvancedMath.Pochhammer(2, 0).Should().Be(1.0);
            AdvancedMath.Pochhammer(2, 1).Should().Be(2.0);
        }

        [TestMethod]
        public void PochhhammerSpecialNegativeIntegerArguments () {
            AdvancedMath.Pochhammer(-1, 2).Should().Be(0.0);
            AdvancedMath.Pochhammer(-1, 1).Should().BeNearly(-1.0);
            AdvancedMath.Pochhammer(-1, 0).Should().BeNearly(1.0);
            AdvancedMath.Pochhammer(-1, -1).Should().BeNearly(-1.0 / 2.0);
            AdvancedMath.Pochhammer(-1, -2).Should().BeNearly(1.0 / 6.0);
            AdvancedMath.Pochhammer(-2, 2).Should().BeNearly(2.0);
            AdvancedMath.Pochhammer(-2, 1).Should().BeNearly(-2.0);
            AdvancedMath.Pochhammer(-2, 0).Should().BeNearly(1.0);
            AdvancedMath.Pochhammer(-2, -1).Should().BeNearly(-1.0 / 3.0);
        }

        [TestMethod]
        public void PochhammerZeroExponent() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(4)) {
                AdvancedMath.Pochhammer(x, 0.0).Should().Be(1.0);
                AdvancedMath.Pochhammer(-x, 0.0).Should().Be(1.0);
            }
        }

        [TestMethod]
        public void PochhammerOneExponent() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(4)) {
                AdvancedMath.Pochhammer(x, 1.0).Should().Be(x);
                AdvancedMath.Pochhammer(-x, 1.0).Should().BeNearly(-x);
            }
        }

        [TestMethod]
        public void PochhammerZeroBase() {
            // This is not true for negative integer y, but integers are measure zero.
            foreach(double y in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(4)) {
                AdvancedMath.Pochhammer(0.0, y).Should().Be(0.0);
                AdvancedMath.Pochhammer(0.0, -y).Should().Be(0.0);
            }
        }

        [TestMethod]
        public void PochhammerZeroBaseIntegerExponent() {
            foreach (int y in TestUtilities.GenerateIntegerValues(1, 100).Take(4)) {
                AdvancedMath.Pochhammer(0.0, y).Should().Be(0.0);
                AdvancedMath.Pochhammer(0.0, -y).Should().BeNearly((y % 2 == 0 ? 1.0 : -1.0) / AdvancedIntegerMath.Factorial(y));
            }
        }

        [TestMethod]
        public void PochhammerFactorialAgreement() {
            // One base -> Factorial
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100).Take(4)) {
                AdvancedMath.Pochhammer(1, n).Should().BeNearly(AdvancedIntegerMath.Factorial(n));
            }
        }

        [TestMethod]
        public void PochhammerDoubleFactorialAgreement () {
            // One-half base -> Double Factorial
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 200).Take(4)) {
                AdvancedMath.Pochhammer(0.5, n).Should().BeNearly(AdvancedIntegerMath.DoubleFactorial(2 * n - 1) / MoreMath.Pow(2.0, n));
            }
        }

        [TestMethod]
        public void PochhammerGammaAgreement() {
            // DLMF 5.2.5 (x)_y = \Gamma(x+y) / \Gamma(x)
            foreach (double x in TestUtilities.GenerateRealValues(0.05, 50.0).Take(4)) {
                foreach (double y in TestUtilities.GenerateRealValues(0.1, 100.0).Take(4)) {
                    AdvancedMath.Pochhammer(x, y).Should().BeNearly(AdvancedMath.Gamma(x + y) / AdvancedMath.Gamma(x));
                }
            }
        }

        [TestMethod]
        public void PochhammerDoubling() {
            // DLMF 5.2.8
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(4)) {
                foreach (int n in TestUtilities.GenerateIntegerValues(5, 50).Take(4)) {
                    AdvancedMath.Pochhammer(2.0 * x, 2 * n).Should().BeNearly(MoreMath.Pow(4.0, n) * AdvancedMath.Pochhammer(x, n) * AdvancedMath.Pochhammer(x + 0.5, n));
                    AdvancedMath.Pochhammer(x, 2 * n).Should().BeNearly(MoreMath.Pow(2.0, 2 * n) * AdvancedMath.Pochhammer(x / 2.0, n) * AdvancedMath.Pochhammer((x + 1.0) / 2.0, n));
                    AdvancedMath.Pochhammer(x, 2 * n + 1).Should().BeNearly(MoreMath.Pow(2.0, 2 * n + 1) * AdvancedMath.Pochhammer(x / 2.0, n + 1) * AdvancedMath.Pochhammer((x + 1.0) / 2.0, n));
                }
            }
        }

        [TestMethod]
        public void PochhammerLowPowers() {
            // DLMF 5.2.4
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4).Take(6)) {
                AdvancedMath.Pochhammer(x, -1).Should().BeNearly(1.0 / (x - 1.0));
                AdvancedMath.Pochhammer(-x, -1).Should().BeNearly(1.0 / (-x - 1.0));
                AdvancedMath.Pochhammer(x, 0).Should().Be(1.0);
                AdvancedMath.Pochhammer(-x, 0).Should().Be(1.0);
                AdvancedMath.Pochhammer(x, 1).Should().Be(x);
                AdvancedMath.Pochhammer(-x, 1).Should().Be(-x);
                AdvancedMath.Pochhammer(x, 2).Should().BeNearly(x * (x + 1.0));
                AdvancedMath.Pochhammer(-x, 2).Should().BeNearly(-x * (-x + 1.0));
            }
        }

        [TestMethod]
        public void PochhammerAddition() {
            // (a)_{m + n} = (a)_m (a + m)_{n}
            Random rng = new Random(1);
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, rng).Take(4)) {
                foreach (int n in TestUtilities.GenerateIntegerValues(1, 50, rng).Take(4)) {
                    foreach (int m in TestUtilities.GenerateIntegerValues(1, 50, rng).Take(4)) {
                            AdvancedMath.Pochhammer(x, m + n).Should().BeNearly(
                                AdvancedMath.Pochhammer(x, m) * AdvancedMath.Pochhammer(x + m, n)
                            );
                    }
                }
            }
        }

        [TestMethod]
        public void PochammerInequality() {
            // Wendel inequality
            // x > 0, 0 < y < 1
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(4)) {
                foreach (double y in TestUtilities.GenerateRealValues(1.0E-3, 1.0).Take(4)) {
                    (AdvancedMath.Pochhammer(x, y) / Math.Pow(x, y)).Should().BeInRange(
                        Math.Pow(x / (x + y), 1.0 - y), 1.0
                    );
                }
            }
        }

        [TestMethod]
        public void PochhammerIntegerReflection() {
            // DLMF 5.2.6
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2).Take(4)) {
                foreach (int n in TestUtilities.GenerateIntegerValues(1, 100).Take(4)) {
                    AdvancedMath.Pochhammer(-x, n).Should().BeNearly((n % 2 == 0 ? 1 : -1) * AdvancedMath.Pochhammer(x - (n - 1), n));
                }
            }
        }

        [TestMethod]
        public void PochhammerIntegerArguments () {
            foreach (int x in TestUtilities.GenerateUniformIntegerValues(-50, +50, 4)) {
                foreach (int y in TestUtilities.GenerateUniformIntegerValues(-50, +50, 4)) {
                    double p = AdvancedMath.Pochhammer(x, y);
                    if (x <= 0) {
                        if (x+y > 0) {
                            p.Should().Be(0.0);
                        } else {
                            p.Should().BeNearly((y % 2 == 0 ? 1 : -1) / AdvancedMath.Pochhammer(1 - x, -y));
                        }
                    } else {
                        if (x + y <= 0) {
                            Assert.IsTrue(Double.IsInfinity(p));
                        } else if (y >= 0) {
                            p.Should().BeNearly(RisingFactorial(x, y));
                        }
                        // Add test for x>0, -x<y<0, and add more rising and falling factorial
                    }
                }
            }


        }

        private double RisingFactorial (int n, int m) {
            double r = 1.0;
            for (int i = 0; i < m; i++) r *= (n + i);
            return r;
        }

        [TestMethod]
        public void PochammerWendelInequality() {
            // Wendel inequality, cited in Qi & Luo, "Bounds for the radio of two Gamma functions",
            // Banch Journal of Mathematical Analysis, 6 (2012) 132-158
            // ( x / x + y)^(1-y) <= (x)_y / x^y <= 1
            // for x > 0, 0 < y < 1
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(4)) {
                foreach (double s in TestUtilities.GenerateRealValues(1.0E-3, 1.0).Take(4)) {
                    (AdvancedMath.Pochhammer(x, s) / Math.Pow(x, s)).Should().BeInRange(
                        Math.Pow(x / (x + s), 1.0 - s), 1.0
                    );
                }
            }
        }

        public void TestNewPochhammer () {
            for (double x = -4.0; x <= 4.0; x += 0.5) {
                for (double y = -4.0; y <= 4.0; y+= 0.5) {
                    double z = AdvancedMath.Pochhammer(x, y);
                    Console.WriteLine($"{x} {y} {z}");
                }
            }
        }

    }
}
