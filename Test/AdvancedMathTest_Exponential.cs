using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using FluentAssertions;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using System.Xml.Serialization;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_Exponential {

        [TestMethod]
        public void IntegralEiSpecialValues() {
            Assert.IsTrue(AdvancedMath.IntegralEi(0.0) == Double.NegativeInfinity);
            Assert.IsTrue(AdvancedMath.IntegralEi(Double.PositiveInfinity) == Double.PositiveInfinity);
            Assert.IsTrue(Double.IsNaN(AdvancedMath.IntegralEi(Double.NaN)));
        }

        [TestMethod]
        public void IntegralEiZero() {
            // Value of zero documented at https://dlmf.nist.gov/6.13
            double x0 = FunctionMath.FindZero(AdvancedMath.IntegralEi, 1.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(x0, 0.37250741078136663446));
        }

        [TestMethod]
        public void IntegralEAtZero() {
            // A & S 5.1.23
            foreach (int n in TestUtilities.GenerateIntegerValues(2, 100, 8)) {
                Assert.IsTrue(AdvancedMath.IntegralE(n, 0.0) == 1.0 / (n - 1.0));
            }
        }

        [TestMethod]
        public void IntegralE0() {
            // A & S 5.1.24
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E3, 8)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.IntegralE(0, x), Math.Exp(-x) / x));
            }
        }

        [TestMethod]
        public void IntegralEDefinition() {
            // A & S 5.1.4
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 10.0, 8)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.IntegralE(n, x),
                        FunctionMath.Integrate(t => Math.Exp(-x * t) / MoreMath.Pow(t, n), 1.0, Double.PositiveInfinity)
                    ));
                }
            }
        }

        [TestMethod]
        public void IntegralESpecialCases() {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.IntegralE(n, 0.0), 1.0 / (n - 1)));
                Assert.IsTrue(AdvancedMath.IntegralE(n, Double.MaxValue) == 0.0);
                Assert.IsTrue(AdvancedMath.IntegralE(n, Double.PositiveInfinity) == 0.0);
                Assert.IsTrue(Double.IsNaN(AdvancedMath.IntegralE(n, Double.NaN)));
            }
        }

        [TestMethod]
        public void IntegralERecurrence() {
            // A & S 5.1.14
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 24)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        n * AdvancedMath.IntegralE(n + 1, x) + x * AdvancedMath.IntegralE(n, x),
                        Math.Exp(-x)
                    ));
                }
            }
        }

        [TestMethod]
        public void IntegralE1Inequality() {
            // A & S 5.1.20
            // Choose maximum x to keep e^x from overflowing.
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, Math.Log(Double.MaxValue / 10.0), 8)) {
                double lower = Math.Log(1.0 + 2.0 / x) / 2.0;
                double upper = Math.Log(1.0 + 1.0 / x);
                double value = Math.Exp(x) * AdvancedMath.IntegralE(1, x);
                Assert.IsTrue((lower < value) && (value < upper));
            }
        }

        [TestMethod]
        public void IntegralEInequality() {
            // A & S 5.1.19
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 8)) {
                // Choose maximum x to keep e^x from overflowing.
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, Math.Log(Double.MaxValue / 10.0), 8)) {
                    double lower = 1.0 / (x + n);
                    double upper = 1.0 / (x + (n - 1));
                    double value = Math.Exp(x) * AdvancedMath.IntegralE(n, x);
                    Assert.IsTrue((lower < value) && (value <= upper));
                }
            }
        }

        [TestMethod]
        public void IntegralE1Integral() {
            // A & S 5.1.33
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate(t => MoreMath.Sqr(AdvancedMath.IntegralE(1, t)), 0.0, Double.PositiveInfinity),
                2.0 * Math.Log(2.0)
            ));
        }

        [TestMethod]
        public void EiE1Integrals() {
            // Here are a couple unusual integrals involving both Ei and E1
            // https://nvlpubs.nist.gov/nistpubs/jres/73B/jresv73Bn3p191_A1b.pdf
            FunctionMath.Integrate(x => Math.Exp(-x) * AdvancedMath.IntegralE(1, x) * AdvancedMath.IntegralEi(x), 0.0, Double.PositiveInfinity).Value.Should().BeNearly(-Math.PI * Math.PI / 12.0);
        }

        [TestMethod]
        public void SiSpecialCases() {
            Assert.IsTrue(AdvancedMath.IntegralSi(Double.NegativeInfinity) == -Math.PI / 2.0);
            Assert.IsTrue(AdvancedMath.IntegralSi(0.0) == 0.0);
            Assert.IsTrue(AdvancedMath.IntegralSi(Double.PositiveInfinity) == Math.PI / 2.0);
            Assert.IsTrue(Double.IsNaN(AdvancedMath.IntegralSi(Double.NaN)));
        }

        [TestMethod]
        public void SiExtremeValues() {
            AdvancedMath.IntegralSi(1.0 / Double.MaxValue).Should().BeNearly(1.0 / Double.MaxValue);
            AdvancedMath.IntegralSi(Double.MaxValue / 2.0).Should().BeNearly(Math.PI / 2.0);
        }

        [TestMethod]
        public void SiReflection() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4).Take(8)) {
                Assert.IsTrue(AdvancedMath.IntegralSi(-x) == -AdvancedMath.IntegralSi(x));
            }
        }

        [TestMethod]
        public void SiDefinition() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(8)) {
                FunctionMath.Integrate(t => Math.Sin(t) / t, 0.0, x).Value.Should().BeNearly(AdvancedMath.IntegralSi(x));
            }
        }

        [TestMethod]
        public void CiSpecialCases() {
            AdvancedMath.IntegralCi(0.0).Should().Be(Double.NegativeInfinity);
            AdvancedMath.IntegralCi(Double.PositiveInfinity).Should().Be(0.0);
            AdvancedMath.IntegralCi(Double.NaN).Should().Be(Double.NaN);
        }

        [TestMethod]
        public void CiExtremeValues() {
            AdvancedMath.IntegralCi(1.0 / Double.MaxValue).Should().BeNearly(Math.Log(1.0 / Double.MaxValue) + AdvancedMath.EulerGamma);
            Math.Abs(AdvancedMath.IntegralCi(Double.MaxValue / 2.0)).Should().BeLessOrEqualTo(2.0 / Double.MaxValue);
        }

        // Ci defining integral does not converge numerically.

        [TestMethod]
        public void CiExpIntegral () {
            // DLMF 6.14.2
            foreach (double x in TestUtilities.GenerateRealValues(0.1, 10.0).Take(4)) {
                FunctionMath.Integrate(t => Math.Exp(-x * t) * AdvancedMath.IntegralCi(t), 0.0, Double.PositiveInfinity).Value.Should().BeNearly(-1.0 / 2.0 * MoreMath.LogOnePlus(x * x) / x);
            }
        }

        [TestMethod]
        public void LittleSiSpecialCases() {
            AdvancedMath.IntegralLittleSi(0.0).Should().Be(-Math.PI / 2.0);
            AdvancedMath.IntegralLittleSi(Double.PositiveInfinity).Should().Be(0.0);
        }

        [TestMethod]
        public void LittleSiBigSiRelation () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(6)) {
                AdvancedMath.IntegralSi(x).Should().BeNearly(AdvancedMath.IntegralLittleSi(x) + Math.PI / 2.0);
            }
        }

        [TestMethod]
        public void LittleSiFromIntegral () {
            // DLMF 6.7.9
            foreach (double x in TestUtilities.GenerateRealValues(0.1, 10.0).Take(4)) {
                FunctionMath.Integrate(t => Math.Exp(-x * MoreMath.Cos(t)) * MoreMath.Cos(x * MoreMath.Sin(t)), 0.0, Math.PI / 2.0).Value.Should().BeNearly(-AdvancedMath.IntegralLittleSi(x));
            }
        }

        // LittleSi squared integral DLMF 6.14.6 does not converge numerically

        [TestMethod]
        public void LittleSiExpIntegral () {
            // DLMF 6.14.3
            foreach (double x in TestUtilities.GenerateRealValues(0.1, 10.0).Take(4)) {
                FunctionMath.Integrate(t => Math.Exp(-x * t) * AdvancedMath.IntegralLittleSi(t), 0.0, Double.PositiveInfinity).Value.Should().BeNearly(-Math.Atan(x) / x);
            }
        }

        // si * ci integral DLMF 6.14.7 does not converge numberically 

        [TestMethod]
        public void CinSpecialCases() {
            AdvancedMath.IntegralCin(0.0).Should().Be(0.0);
            AdvancedMath.IntegralCin(Double.PositiveInfinity).Should().Be(Double.PositiveInfinity);
            AdvancedMath.IntegralCin(Double.NegativeInfinity).Should().Be(Double.PositiveInfinity);
        }

        [TestMethod]
        public void CinReflection() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4).Take(8)) {
                Assert.IsTrue(AdvancedMath.IntegralCin(-x) == AdvancedMath.IntegralCin(x));
            }
        }

        [TestMethod]
        public void CinDefinition() {
            // To avoid cancellation error, near cos(x) ~ 1, use half-angle formula
            //   1 - \cos(x) = 2 \sin^2(x/2)
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(8)) {
                FunctionMath.Integrate(t => 2.0 * MoreMath.Sqr(MoreMath.Sin(t / 2.0)) / t, 0.0, x).Value.Should().BeNearly(AdvancedMath.IntegralCin(x));
            }
        }

        [TestMethod]
        public void IntegralCiCinRelation() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0, 4)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.EulerGamma + Math.Log(x) - AdvancedMath.IntegralCin(x),
                    AdvancedMath.IntegralCi(x)
                ));
            }
            foreach (double x in TestUtilities.GenerateRealValues(1.0, 1.0E4, 4)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.EulerGamma + Math.Log(x) - AdvancedMath.IntegralCi(x),
                    AdvancedMath.IntegralCin(x)
                ));
            }
        }

        [TestMethod]
        public void IntegralCiSiExpIntegrals() {

            // these integrals are from Oldham et al, An Atlas of Functions
            // and Wolfram functions site

            // A & S 5.2.28
            FunctionMath.Integrate(t => AdvancedMath.IntegralCi(t) * Math.Exp(-t), 0.0, Double.PositiveInfinity).Value.Should().BeNearly(-Math.Log(2.0) / 2.0);

            FunctionMath.Integrate(t => AdvancedMath.IntegralSi(t) * Math.Exp(-t), 0.0, Double.PositiveInfinity).Value.Should().BeNearly(Math.PI / 4.0);

            /*
            // the integral of [Ci(t)]^2 does not converge numerically 
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate(t => MoreMath.Sqr(AdvancedMath.IntegralCi(t)), Interval.FromEndpoints(0.0, Double.PositiveInfinity)),
                Math.PI / 2.0
            ));
            */
        }

        [TestMethod]
        public void ShiSpecialCases() {
            AdvancedMath.IntegralShi(0.0).Should().Be(0.0);
            AdvancedMath.IntegralShi(Double.PositiveInfinity).Should().Be(Double.PositiveInfinity);
            AdvancedMath.IntegralShi(Double.NegativeInfinity).Should().Be(Double.NegativeInfinity);
            AdvancedMath.IntegralShi(Double.NaN).Should().Be(Double.NaN);
        }

        [TestMethod]
        public void ShiExtremeValues() {
            AdvancedMath.IntegralShi(1.0 / Double.MaxValue).Should().Be(1.0 / Double.MaxValue);
            AdvancedMath.IntegralShi(Double.MaxValue).Should().Be(Double.PositiveInfinity);
        }

        [TestMethod]
        public void ShiReflection() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(4)) {
                Assert.IsTrue(AdvancedMath.IntegralShi(-x) == -AdvancedMath.IntegralShi(x));
            }
        }

        [TestMethod]
        public void ShiDefintion() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E1).Take(4)) {
                FunctionMath.Integrate(t => Math.Sinh(t) / t, 0.0, x).Value.Should().BeNearly(AdvancedMath.IntegralShi(x));
            }
        }

        [TestMethod]
        public void ShiExpIntegral() {
            // \int_0^{\infty} \! dt \, shi(x) e^{-a t} = arctanh(a) / a according to Mathematica
            FunctionMath.Integrate(t => AdvancedMath.IntegralShi(t) * Math.Exp(-2.0 * t), 0.0, Double.PositiveInfinity).Value.Should().BeNearly(Math.Log(3.0) / 4.0);
        }

        [TestMethod]
        public void ShiEiE1Relation () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(5)) {
                AdvancedMath.IntegralShi(x).Should().BeNearly((AdvancedMath.IntegralEi(x) + AdvancedMath.IntegralE(1, x)) / 2.0);
            }
        }

        [TestMethod]
        public void ChiSpecialValues() {
            AdvancedMath.IntegralChi(0.0).Should().Be(Double.NegativeInfinity);
            AdvancedMath.IntegralChi(Double.PositiveInfinity).Should().Be(Double.PositiveInfinity);
            AdvancedMath.IntegralChi(Double.NaN).Should().Be(Double.NaN);
        }

        [TestMethod]
        public void ChiExtremeValues() {
            AdvancedMath.IntegralChi(1.0 / Double.MaxValue).Should().Be(Math.Log(1.0 / Double.MaxValue) + AdvancedMath.EulerGamma);
            AdvancedMath.IntegralShi(Double.MaxValue).Should().Be(Double.PositiveInfinity);
        }

        [TestMethod]
        public void ChiExpIntegral() {
            // \int_0^{\infty} \! dt \, chi(x) e^{-a t} = -\frac{\ln(a^2 - 1)}{2a}
            FunctionMath.Integrate(t => AdvancedMath.IntegralChi(t) * Math.Exp(-2.0 * t), 0.0, Double.PositiveInfinity).Value.Should().BeNearly(-Math.Log(3.0) / 4.0);
        }

        [TestMethod]
        public void ChiRoot() {
            FunctionMath.FindZero(t => AdvancedMath.IntegralChi(t), Interval.FromEndpoints(0.1, 10.0)).Should().BeNearly(0.523822571389864);
        }

        [TestMethod]
        public void ChiEiE1Relation() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(5)) {
                AdvancedMath.IntegralChi(x).Should().BeNearly((AdvancedMath.IntegralEi(x) - AdvancedMath.IntegralE(1, x)) / 2.0);
            }
        }

    }
}
