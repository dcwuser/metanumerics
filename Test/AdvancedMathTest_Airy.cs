using System;
using System.Linq;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using FluentAssertions;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Analysis;
using System.Diagnostics;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_Airy {


        [TestMethod]
        public void AiryAiExtremeValues () {
            AdvancedMath.AiryAi(Double.NegativeInfinity).Should().Be(Double.NaN);
            AdvancedMath.AiryAi(Double.MaxValue).Should().Be(0.0);
            AdvancedMath.AiryAi(Double.PositiveInfinity).Should().Be(0.0);
            AdvancedMath.AiryAi(Double.NaN).Should().Be(Double.NaN);
        }

        [TestMethod]
        public void AiryBiExtremeValues () {
            Assert.IsTrue(Double.IsNaN(AdvancedMath.AiryBi(Double.NegativeInfinity)));
            Assert.IsTrue(AdvancedMath.AiryBi(Double.MaxValue) == Double.PositiveInfinity);
            Assert.IsTrue(AdvancedMath.AiryBi(Double.PositiveInfinity) == Double.PositiveInfinity);
            Assert.IsTrue(Double.IsNaN(AdvancedMath.AiryBi(Double.NaN)));
        }

        [TestMethod]
        public void AiryAiIntegral () {
            // DLMF 9.10.11
            FunctionMath.Integrate(t => AdvancedMath.AiryAi(t), 0.0, Double.PositiveInfinity).Value.Should().BeNearly(1.0 / 3.0);
        }

        [TestMethod]
        public void AiryZeroArgument () {
            // DLMF 9.2.3 - 9.2.6
            AdvancedMath.AiryAi(0.0).Should().BeNearly(1.0 / Math.Pow(3.0, 2.0 / 3.0) / AdvancedMath.Gamma(2.0 / 3.0));
            AdvancedMath.AiryBi(0.0).Should().BeNearly(1.0 / Math.Pow(3.0, 1.0 / 6.0) / AdvancedMath.Gamma(2.0 / 3.0));

            SolutionPair s = AdvancedMath.Airy(0.0);
            s.FirstSolutionValue.Should().BeNearly(1.0 / Math.Pow(3.0, 2.0 / 3.0) / AdvancedMath.Gamma(2.0 / 3.0));
            s.FirstSolutionDerivative.Should().BeNearly(-1.0 / Math.Pow(3.0, 1.0 / 3.0) / AdvancedMath.Gamma(1.0 / 3.0));
            s.SecondSolutionValue.Should().BeNearly(1.0 / Math.Pow(3.0, 1.0 / 6.0) / AdvancedMath.Gamma(2.0 / 3.0));
            s.SecondSolutionDerivative.Should().BeNearly(Math.Pow(3.0, 1.0 / 6.0) / AdvancedMath.Gamma(1.0 / 3.0));
        }

        [TestMethod]
        public void AiryPowerIntegrals () {

            // Bernard Laurenzi, "Moment Integrals of Powers of Airy Functions", ZAMP 44 (1993) 891

            double I2 = FunctionMath.Integrate(t => MoreMath.Pow(AdvancedMath.AiryAi(t), 2), 0.0, Double.PositiveInfinity);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(I2, MoreMath.Sqr(AdvancedMath.Airy(0.0).FirstSolutionDerivative)));

            // I3 involves difference of 2F1 functions

            double I4 = FunctionMath.Integrate(t => MoreMath.Pow(AdvancedMath.AiryAi(t), 4), 0.0, Double.PositiveInfinity);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(I4, Math.Log(3.0) / (24.0 * Math.PI * Math.PI)));

        }

        [TestMethod]
        public void AiryModulusMomentIntegrals () {

            // https://math.stackexchange.com/questions/507425/an-integral-involving-airy-functions-int-0-infty-fracxp-operatornameai

            double I0 = FunctionMath.Integrate(t => {
                SolutionPair s = AdvancedMath.Airy(t);
                return (1.0 / (MoreMath.Sqr(s.FirstSolutionValue) + MoreMath.Sqr(s.SecondSolutionValue)));
            }, 0.0, Double.PositiveInfinity);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(I0, Math.PI * Math.PI / 6.0));

            double I3 = FunctionMath.Integrate(t => {
                SolutionPair s = AdvancedMath.Airy(t);
                return (MoreMath.Pow(t, 3) / (MoreMath.Sqr(s.FirstSolutionValue) + MoreMath.Sqr(s.SecondSolutionValue)));
            }, 0.0, Double.PositiveInfinity);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(I3, 5.0 * Math.PI * Math.PI / 32.0));

        }

        [TestMethod]
        public void AiryBairyIntegral () {

            // for small x: the integral will not converse for x <~ 5 because of the 1/sqrt(t) divergence near 0
            // for large x: Bi explodes exponentially
            // but in between, it is a decent test of Ai and Bi together

            foreach (double x in TestUtilities.GenerateRealValues(5.0, 50.0).Take(8)) {

                double I = FunctionMath.Integrate(
                    t => { return (Math.Exp(x * t - t * t * t / 12.0) / Math.Sqrt(t)); },
                    Interval.FromEndpoints(0.0, Double.PositiveInfinity)
                );

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    MoreMath.Pow(AdvancedMath.AiryAi(x), 2) + MoreMath.Pow(AdvancedMath.AiryBi(x), 2), I / Math.Pow(Math.PI, 3.0 / 2.0)
                ));

            }

        }

        [TestMethod]
        public void AiryWronskian () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(8)) {
                SolutionPair p = AdvancedMath.Airy(x);
                (1.0 / Math.PI).Should().BeNearlySumOf(p.FirstSolutionValue * p.SecondSolutionDerivative, -p.SecondSolutionValue * p.FirstSolutionDerivative);
                SolutionPair q = AdvancedMath.Airy(-x);
                (1.0 / Math.PI).Should().BeNearlySumOf(q.FirstSolutionValue * q.SecondSolutionDerivative, -q.SecondSolutionValue * q.FirstSolutionDerivative);
            }
        }

        [TestMethod]
        public void AiryAgreement () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(4)) {
                SolutionPair sp = AdvancedMath.Airy(x);
                sp.FirstSolutionValue.Should().BeNearly(AdvancedMath.AiryAi(x));
                sp.SecondSolutionValue.Should().BeNearly(AdvancedMath.AiryBi(x));
                SolutionPair sn = AdvancedMath.Airy(-x);
                sn.FirstSolutionValue.Should().BeNearly(AdvancedMath.AiryAi(-x));
                sn.SecondSolutionValue.Should().BeNearly(AdvancedMath.AiryBi(-x));
            }
        }

        [TestMethod]
        public void ComplexAiryAgreement () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2).Take(6)) {
                Complex ap = AdvancedComplexMath.AiryAi(x);
                ap.Re.Should().BeNearly(AdvancedMath.AiryAi(x));
                ap.Im.Should().Be(0.0);
                Complex an = AdvancedComplexMath.AiryAi(-x);
                an.Re.Should().BeNearly(AdvancedMath.AiryAi(-x));
                an.Im.Should().Be(0.0);
                Complex bp = AdvancedComplexMath.AiryBi(x);
                bp.Re.Should().BeNearly(AdvancedMath.AiryBi(x));
                bp.Im.Should().Be(0.0);
                Complex bn = AdvancedComplexMath.AiryBi(-x);
                bn.Re.Should().BeNearly(AdvancedMath.AiryBi(-x));
                bn.Im.Should().Be(0.0);
            }
        }

        [TestMethod]
        public void ComplexAiryIdentity () {
            // DLMF 9.2.12
            Complex omega = ComplexMath.Exp(2.0 / 3.0 * Math.PI * Complex.I);
            foreach (Complex z in TestUtilities.GenerateComplexValues(min: 1.0E-2, max: 1.0E2).Take(12)) {
                TestUtilities.IsSumNearlyEqual(new Complex[] {
                    AdvancedComplexMath.AiryAi(z) + omega * AdvancedComplexMath.AiryAi(omega * z) + AdvancedComplexMath.AiryAi(z / omega) / omega
                }, Complex.Zero);
            }
        }

        [TestMethod]
        public void ComplexBairyIdentity() {
            // DLMF 9.2.10
            Complex omega3 = ComplexMath.Exp(2.0 * Math.PI * Complex.I / 3.0);
            Complex omega12 = ComplexMath.Exp(2.0 * Math.PI * Complex.I / 12.0);
            foreach (Complex z in TestUtilities.GenerateComplexValues(min: 1.0E-2, max: 1.0E2).Take(10)) {
                TestUtilities.IsSumNearlyEqual(new Complex[] {
                    omega12 * AdvancedComplexMath.AiryAi(omega3 * z) + omega12.Conjugate * AdvancedComplexMath.AiryAi(omega3.Conjugate * z)
                }, AdvancedComplexMath.AiryBi(z));
            }
        }

        [TestMethod]
        public void ComplexAiryNegation() {
            // DLMF 9.2.14 & 9.2.15
            Complex omega3 = ComplexMath.Exp(2.0 * Math.PI * Complex.I / 3.0);
            Complex omega6 = ComplexMath.Exp(2.0 * Math.PI * Complex.I / 6.0);
            Complex omega12 = ComplexMath.Exp(2.0 * Math.PI * Complex.I / 12.0);
            foreach (Complex z in TestUtilities.GenerateComplexValues(min: 1.0E-2, max: 1.0E2).Take(8)) {
                TestUtilities.IsSumNearlyEqual(new Complex[] {
                    omega6 * AdvancedComplexMath.AiryAi(omega6 * z) + omega6.Conjugate * AdvancedComplexMath.AiryAi(omega6.Conjugate * z)
                }, AdvancedComplexMath.AiryAi(-z));
                TestUtilities.IsSumNearlyEqual(new Complex[] {
                    omega12.Conjugate * AdvancedComplexMath.AiryAi(omega6 * z) + omega12 * AdvancedComplexMath.AiryAi(omega6.Conjugate * z)
                }, AdvancedComplexMath.AiryBi(-z));
            }
        }

        [TestMethod]
        public void ComplexAiryConjugationSymmetry () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-1, 1.0E2).Take(6)) {
                AdvancedComplexMath.AiryAi(z.Conjugate).Should().Be(AdvancedComplexMath.AiryAi(z).Conjugate);
                AdvancedComplexMath.AiryBi(z.Conjugate).Should().Be(AdvancedComplexMath.AiryBi(z).Conjugate);
            }
        }

    }
}
