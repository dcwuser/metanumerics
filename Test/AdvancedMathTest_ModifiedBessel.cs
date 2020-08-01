using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_ModifiedBessel {

        [TestMethod]
        public void ModifiedBesselAtZero () {
            Assert.IsTrue(AdvancedMath.ModifiedBesselI(0, 0.0) == 1.0);
            Assert.IsTrue(AdvancedMath.ModifiedBesselI(1, 0.0) == 0.0);
            Assert.IsTrue(AdvancedMath.ModifiedBesselK(0, 0.0) == Double.PositiveInfinity);
            Assert.IsTrue(AdvancedMath.ModifiedBesselK(1, 0.0) == Double.PositiveInfinity);
        }

        [TestMethod]
        public void ModifiedBesselInequalities () {

            // I and K are positive and increase monotonically with x
            // I decreases monotonically with nu and K increases monotonically with nu
            foreach (double nu in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 4)) {
                    double I = AdvancedMath.ModifiedBesselI(nu, x);
                    Assert.IsTrue(I > 0.0);
                    double K = AdvancedMath.ModifiedBesselK(nu, x);
                    Assert.IsTrue(K > 0.0);
                    foreach (double y in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 2, 2)) {
                        double I1 = AdvancedMath.ModifiedBesselI(nu, y);
                        Assert.IsTrue((x < y) ^ !(I < I1));
                        double K1 = AdvancedMath.ModifiedBesselK(nu, y);
                        Assert.IsTrue((x < y) ^ !(K > K1));
                    }
                }
            }

        }

        [TestMethod]
        public void ScaledModifiedBesselWronskian () {
            // A&S 9.6.15
            // The same relation holds for scaled values, because each term gets an e^{x} and e^{-x}
            // We therefore use this to test the scaled functions at arguments for which the unscaled functions
            // would over- and under-flow.
            foreach (double nu in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E6, 8)) {

                    double I = AdvancedMath.ScaledModifiedBesselI(nu, x);
                    double Ip1 = AdvancedMath.ScaledModifiedBesselI(nu + 1.0, x);
                    double K = AdvancedMath.ScaledModifiedBesselK(nu, x);
                    double Kp1 = AdvancedMath.ScaledModifiedBesselK(nu + 1.0, x);
                    // No need to use IsSumNearlyEqual; both terms are positive so no cancelation is possible.
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(I * Kp1 + Ip1 * K, 1.0 / x));

                }
            }
        }

        [TestMethod]
        public void ModifiedBesselHalfIntegerOrder () {

            // A&S 10.2.13, 10.2.14

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 8)) {

                double F = Math.Sqrt(Math.PI / 2.0 / x);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(F * AdvancedMath.ModifiedBesselI(0.5, x), Math.Sinh(x) / x));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(F * AdvancedMath.ModifiedBesselI(-0.5, x), Math.Cosh(x) / x));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.ModifiedBesselK(0.5, x), F * Math.Exp(-x)));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.ModifiedBesselK(1.5, x), F * Math.Exp(-x) * (1.0 + 1.0 / x)));

            }

        }

        [TestMethod]
        public void ModifiedBesselTowers () {

            // e^z = I_0 + 2 I_1 + 2 I_2 + 2 I_3 + ...
            // cosh z = I_0 2 + 2 I_2 + 2 I_4 + 2 I_6 + ...
            // sinh z = 2 I_1 + 2 I_3 + 2 I_5

            // it takes about n ~ x terms to converge, so don't let x get too big

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 16)) {

                double I = AdvancedMath.ModifiedBesselI(0, x);

                double exp = I;
                double cosh = I;
                double sinh = 0.0;

                for (int n = 1; n < 100; n++) {

                    double exp_old = exp;
                    I = 2.0 * AdvancedMath.ModifiedBesselI(n, x);
                    exp += I;
                    if (n % 2 == 0) {
                        cosh += I;
                    } else {
                        sinh += I;
                    }
                    if (exp == exp_old) {
                        break;
                    }
                }

                Assert.IsTrue(TestUtilities.IsNearlyEqual(exp, Math.Exp(x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(cosh, Math.Cosh(x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(sinh, Math.Sinh(x)));

            }

        }

        [TestMethod]
        public void ModifiedBesselGeneratingFunction () {

            // DLMF 10.35.1	
            // e^{(z/2)(t + 1/t)} = \sum_{k=-\infty}^{+\infty} t^k I_k(z)

            foreach (double t in TestUtilities.GenerateRealValues(1.0E-1, 1.0E-1, 4)) {
                foreach (double z in TestUtilities.GenerateRealValues(1.0E-2, 1.0E1, 4)) {

                    double tn = 1.0;
                    double S = AdvancedMath.ModifiedBesselI(0, z);

                    for (int n = 1; n < 100; n++) {
                        double S_old = S;
                        tn *= t;
                        S += tn * AdvancedMath.ModifiedBesselI(n, z) + AdvancedMath.ModifiedBesselI(-n, z) / tn;
                        if (S == S_old) goto Test;
                    }
                    throw new NonconvergenceException();

                    Test:
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        S, Math.Exp(z / 2.0 * (t + 1.0 / t))
                    ));
                }
            }

        }

        [TestMethod]
        public void ModifiedBesselIntegral () {

            // DLMF 10.32.1
            // I_0(z) = \frac{1}{\pi} \int_{0}^{\pi} \! d\theta \, e^{z \cos(\theta)}  

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 5)) {

                double I = FunctionMath.Integrate(
                    (double t) => Math.Exp(x * Math.Cos(t)),
                    Interval.FromEndpoints(0.0, Math.PI)
                );

                Assert.IsTrue(TestUtilities.IsNearlyEqual(I / Math.PI, AdvancedMath.ModifiedBesselI(0, x)));

            }

        }

        [TestMethod]
        public void ModifiedBesselAgreement () {

            foreach (double nu in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 16)) {

                    SolutionPair s = AdvancedMath.ModifiedBessel(nu, x);
                    double I = AdvancedMath.ModifiedBesselI(nu, x);
                    double K = AdvancedMath.ModifiedBesselK(nu, x);

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(s.FirstSolutionValue, AdvancedMath.ModifiedBesselI(nu, x)));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(s.SecondSolutionValue, AdvancedMath.ModifiedBesselK(nu, x)));

                }
            }

        }

        [TestMethod]
        public void ModifiedBesselDerivativeWronskian () {

            foreach (double nu in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 16)) {

                    SolutionPair s = AdvancedMath.ModifiedBessel(nu, x);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        s.SecondSolutionValue * s.FirstSolutionDerivative - s.FirstSolutionValue * s.SecondSolutionDerivative,
                        1.0 / x
                    ));
                    // there is no cancelation, because I, I', K > 0 and K' < 0

                }
            }

        }

        [TestMethod]
        public void BesselModifiedBesselRelationship () {

            // This is a very interesting test because it relates (normal Bessel) J and (modified Bessel) I
            //   I_{\nu}(x) = \sum_{k=0}^{\infty} J_{\nu + k}(x) \frac{x^k}{k!}
            // We want x not too big, so that \frac{x^k}{k!} converges, and x <~ \nu, so that J_{\nu+k} decreases rapidly

            foreach (double nu in TestUtilities.GenerateRealValues(0.1, 10.0, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(0.1, nu, 4)) {

                    double I = 0.0;
                    for (int k = 0; k < 100; k++) {

                        double I_old = I;
                        I += AdvancedMath.BesselJ(nu + k, x) * MoreMath.Pow(x, k) / AdvancedIntegerMath.Factorial(k);
                        if (I == I_old) goto Test;

                    }
                    throw new NonconvergenceException();

                    Test:
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(I, AdvancedMath.ModifiedBesselI(nu, x)));

                }
            }

        }

    }
}
