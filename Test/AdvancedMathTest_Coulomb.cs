using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_Coulomb {

        [TestMethod]
        public void BesselInflectionValue () {

            // Abromowitz derived a series for J_{\nu}(\nu) and Y_{\nu}(\nu),
            // the value at the inflection point.

            // These appear in A&S 9.3.31 - 9.3.34 with numerical values for the coefficients.
            // Symbolic values for the coefficients appear in Jentschura & Loetstedt 2011
            // (https://arxiv.org/abs/1112.0072)

            // This is a fairly onerous test to code, but it is a good one, because
            // the Bessel functions are difficult to calculate in this region.

            double A = Math.Pow(2.0, 1.0 / 3.0) / Math.Pow(3.0, 2.0 / 3.0) / AdvancedMath.Gamma(2.0 / 3.0);

            double B = Math.Pow(2.0, 2.0 / 3.0) / Math.Pow(3.0, 1.0 / 3.0) / AdvancedMath.Gamma(1.0 / 3.0);

            double[] a = new double[] {
                1.0,
                -1.0 / 225.0,
                151439.0 / 218295000.0,
                -887278009.0 / 2504935125000.0,
                1374085664813273149.0 / 3633280647121125000000.0
            };

            double[] b = new double[] {
                1.0 / 70.0,
                -1213.0 / 1023750.0,
                16542537833.0 / 37743205500000.0,
                -9597171184603.0 / 25476663712500000.0,
                53299328587804322691259.0 / 91182706744837207500000000.0
            };

            foreach (double nu in TestUtilities.GenerateRealValues(10.0, 100.0, 8)) {

                double Af = A / Math.Pow(nu, 1.0 / 3.0);
                double Bf = B / Math.Pow(nu, 5.0 / 3.0);

                double J = 0.0;
                double dJ = 0.0;

                for (int k = 0; k < Math.Min(a.Length, b.Length); k++) {

                    double nuk = MoreMath.Pow(nu, 2 * k);

                    dJ = Af * a[k] / nuk;
                    J += dJ;

                    dJ = Bf * b[k] / nuk;
                    J -= dJ;

                }

                EvaluationSettings e = new EvaluationSettings() {
                    RelativePrecision = TestUtilities.TargetPrecision,
                    AbsolutePrecision = Math.Abs(dJ)
                };

                TestUtilities.IsNearlyEqual(AdvancedMath.BesselJ(nu, nu), J, e);

            }

        }

        [TestMethod]
        public void CoulombRhoZero () {

            // L = 0
            foreach (double eta in TestUtilities.GenerateRealValues(1.0E-2, 1.0E3, 4)) {
                foreach (int sign in new int[] { 1, -1 }) {
                    // F is zero
                    Assert.IsTrue(AdvancedMath.CoulombF(0, sign * eta, 0.0) == 0.0);
                    // G is 1/C_{\ell}(\eta), F' = C_{\ell}(eta), so product is 1
                    SolutionPair s = AdvancedMath.Coulomb(0, sign * eta, 0.0);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(s.FirstSolutionDerivative * s.SecondSolutionValue, 1.0));
                    Assert.IsTrue(Double.IsNegativeInfinity(s.SecondSolutionDerivative));
                }
            }


            // L > 0
            foreach (int L in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (int sign in new int[] {+1, -1}) {
                    foreach (double eta in TestUtilities.GenerateRealValues(1.0E-2, 1.0E3, 8)) {
                        // F is zero, G is infinite
                        Assert.IsTrue(AdvancedMath.CoulombF(L, sign * eta, 0.0) == 0.0);
                        Assert.IsTrue(Double.IsInfinity(AdvancedMath.CoulombG(L, sign * eta, 0.0)));
                    }
                }
            }

        }


        [TestMethod]
        public void CoulombEtaZero () {

            // For \eta = 0, Coulomb wave functions reduce to spherical Bessel functions

            foreach (int L in TestUtilities.GenerateIntegerValues(1, 100, 8)) {
                foreach (double rho in TestUtilities.GenerateRealValues(1.0E-3, 1.0E5, 16)) {

                    double F = AdvancedMath.CoulombF(L, 0, rho);
                    double j = AdvancedMath.SphericalBesselJ(L, rho);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(F, rho * j));

                    double G = AdvancedMath.CoulombG(L, 0, rho);
                    double y = AdvancedMath.SphericalBesselY(L, rho);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(G, -rho * y, 4.0 * TestUtilities.TargetPrecision));

                }
            }

        }

        [TestMethod]
        public void CoulombAgreement () {

            // The individual CoulombF and CoulombG functions should agree with the Coulomb
            // function that computes both.

            foreach (int L in TestUtilities.GenerateIntegerValues(1, 100, 8)) {
                foreach (int sign in new int[] { +1, -1 }) {
                    foreach (double eta in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 16)) {
                        foreach (double rho in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 32)) {

                            double F = AdvancedMath.CoulombF(L, sign * eta, rho);
                            double G = AdvancedMath.CoulombG(L, sign * eta, rho);
                            SolutionPair s = AdvancedMath.Coulomb(L, sign * eta, rho);

                            // Deep in transition region, there are cases where we integrate for full solution pair in order to get G,
                            // then use Wronskian from that G to get F. But if we are just computing F, we use series to get a different
                            // (and more accurate) answer. So we relax agreement criteria a bit. Eventually, would be great to do
                            // series for G everywhere we can do series for F.

                            if (!TestUtilities.IsNearlyEqual(F, s.FirstSolutionValue, 4.0 * TestUtilities.TargetPrecision)) {
                                Console.WriteLine($"L={L} eta={eta} rho={rho} F={F} s.F={s.FirstSolutionValue} b={Double.IsNaN(s.SecondSolutionDerivative)}");
                            }

                            Assert.IsTrue(TestUtilities.IsNearlyEqual(F, s.FirstSolutionValue, 4.0 * TestUtilities.TargetPrecision));

                            Assert.IsTrue(TestUtilities.IsNearlyEqual(G, s.SecondSolutionValue));
                        }
                    }
                }
            }
        }

        [TestMethod]
        public void CoulombRecursion () {

            foreach (int L in TestUtilities.GenerateIntegerValues(1, 100, 8)) {
                foreach (double eta in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 16)) {
                    foreach (double rho in TestUtilities.GenerateRealValues(1.0E-3, 1.0E3, 32)) {

                        CoulombRecursionHelper(L, eta, rho,
                            AdvancedMath.CoulombF(L - 1, eta, rho), AdvancedMath.CoulombF(L, eta, rho), AdvancedMath.CoulombF(L + 1, eta, rho)
                        );

                        CoulombRecursionHelper(L, -eta, rho,
                            AdvancedMath.CoulombF(L - 1, -eta, rho), AdvancedMath.CoulombF(L, -eta, rho), AdvancedMath.CoulombF(L + 1, -eta, rho)
                        );

                        CoulombRecursionHelper(L, eta, rho,
                            AdvancedMath.CoulombG(L - 1, eta, rho), AdvancedMath.CoulombG(L, eta, rho), AdvancedMath.CoulombG(L + 1, eta, rho)
                        );

                        CoulombRecursionHelper(L, -eta, rho,
                            AdvancedMath.CoulombG(L - 1, -eta, rho), AdvancedMath.CoulombG(L, -eta, rho), AdvancedMath.CoulombG(L + 1, -eta, rho)
                        );

                    }
                }
            }


        }

        private static void CoulombRecursionHelper (double L, double eta, double rho, double UM, double U, double UP) {

            // Had issues with comparison failing with subnormal numbers that can't possibly deliver target accuracy. Skip 'em.
            const double subnormal = 1.0 / Double.MaxValue;
            if ((Math.Abs(UM) < subnormal) && (Math.Abs(U) < subnormal) && (Math.Abs(UP) < subnormal)) return;

            double am = (L + 1) * Math.Sqrt(L * L + eta * eta);
            double a = (2 * L + 1) * (eta + L * (L + 1) / rho);
            double ap = L * Math.Sqrt((L + 1) * (L + 1) + eta * eta);

            Assert.IsTrue(TestUtilities.IsSumNearlyEqual(am * UM, ap * UP, a * U, 4.0 * TestUtilities.TargetPrecision));

        }

        [TestMethod]
        public void CoulombRecursionWithDerivatives () {

            // Abromowitz & Stegun 14.2.1-14.2.3 recursion relations

            foreach (int L in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (int sign in new int[] {+1 , -1}) {
                    foreach (double eta in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 8)) {
                        foreach (double rho in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 16)) {
                            CoulombRecursionWithDerivativesHelper(L, sign * eta, rho);
                        }
                    }
                }
            }

        }

        private void CoulombRecursionWithDerivativesHelper (int L, double eta, double rho) {

            SolutionPair sm = AdvancedMath.Coulomb(L - 1, eta, rho);
            SolutionPair s0 = AdvancedMath.Coulomb(L, eta, rho);
            SolutionPair sp = AdvancedMath.Coulomb(L + 1, eta, rho);

            double eps = 8.0 * TestUtilities.TargetPrecision;

            // Relating u_{L}' to u_{L-1} and u_{L}

            Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                MoreMath.Hypot(L, eta) * sm.FirstSolutionValue,
                -(L * L / rho + eta) * s0.FirstSolutionValue,
                L * s0.FirstSolutionDerivative,
                eps
            ));

            Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                MoreMath.Hypot(L, eta) * sm.SecondSolutionValue,
                -(L * L / rho + eta) * s0.SecondSolutionValue,
                L * s0.SecondSolutionDerivative,
                eps
            ));

            // Relating u_{L}' to u_{L} and u_{L+1}

            // Relating u_{L+1}, u_{L}, u_{L-1}

            Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                (2 * L + 1) * (eta + L * (L + 1) / rho) * s0.FirstSolutionValue,
                -(L + 1) * MoreMath.Hypot(L, eta) * sm.FirstSolutionValue,
                L * MoreMath.Hypot(L + 1, eta) * sp.FirstSolutionValue,
                eps
            ));

            Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                (2 * L + 1) * (eta + L * (L + 1) / rho) * s0.SecondSolutionValue,
                -(L + 1) * MoreMath.Hypot(L, eta) * sm.SecondSolutionValue,
                L * MoreMath.Hypot(L + 1, eta) * sp.SecondSolutionValue,
                eps
            ));

        }

        [TestMethod]
        public void CoulombWronskianWithDerivatives () {

            foreach (int L in TestUtilities.GenerateIntegerValues(1, 100, 8)) {
                foreach (int sign in new int[] {+1, -1}) {
                    foreach (double eta in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 16)) {
                        foreach (double rho in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 32)) {

                            SolutionPair s = AdvancedMath.Coulomb(L, eta, rho);

                            // If F and G underflow and overflow, skip check
                            if (s.FirstSolutionValue == 0.0) continue;

                            Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                                s.FirstSolutionDerivative * s.SecondSolutionValue, -s.FirstSolutionValue * s.SecondSolutionDerivative, 1.0,
                                TestUtilities.TargetPrecision * 10.0
                            ));

                        }
                    }
                }
            }

        }

        [TestMethod]
        public void CoulombInequality () {

        }

    }

}
