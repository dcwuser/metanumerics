using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Analysis;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_Airy {


        [TestMethod]
        public void AiryAiExtremeValues () {
            Assert.IsTrue(Double.IsNaN(AdvancedMath.AiryAi(Double.NegativeInfinity)));
            Assert.IsTrue(AdvancedMath.AiryAi(Double.MaxValue) == 0.0);
            Assert.IsTrue(AdvancedMath.AiryAi(Double.PositiveInfinity) == 0.0);
            Assert.IsTrue(Double.IsNaN(AdvancedMath.AiryAi(Double.NaN)));
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
            Func<double, double> f = delegate (double t) {
                return (AdvancedMath.AiryAi(t));
            };
            Interval r = Interval.FromEndpoints(0.0, Double.PositiveInfinity);
            double I = FunctionMath.Integrate(f, r);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(I, 1.0 / 3.0));

        }

        [TestMethod]
        public void AiryZeroArgument () {

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.AiryAi(0.0), 1.0 / Math.Pow(3.0, 2.0 / 3.0) / AdvancedMath.Gamma(2.0 / 3.0)
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.AiryBi(0.0), 1.0 / Math.Pow(3.0, 1.0 / 6.0) / AdvancedMath.Gamma(2.0 / 3.0)
            ));
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

            foreach (double x in TestUtilities.GenerateRealValues(5.0, 50.0, 8)) {

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

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 8)) {

                SolutionPair p = AdvancedMath.Airy(x);
                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(p.FirstSolutionValue * p.SecondSolutionDerivative, -p.SecondSolutionValue * p.FirstSolutionDerivative, 1.0 / Math.PI));

                SolutionPair q = AdvancedMath.Airy(-x);
                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(q.FirstSolutionValue * q.SecondSolutionDerivative, -q.SecondSolutionValue * q.FirstSolutionDerivative, 1.0 / Math.PI));

            }

        }

        [TestMethod]
        public void AiryAgreement () {

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 4)) {
                SolutionPair s = AdvancedMath.Airy(-x);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s.FirstSolutionValue, AdvancedMath.AiryAi(-x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s.SecondSolutionValue, AdvancedMath.AiryBi(-x)));
            }

        }


    }
}
