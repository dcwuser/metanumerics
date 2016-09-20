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
        public void CoulombSpecialValue () {

            // F is zero at origin for all L and eta
            Assert.IsTrue(AdvancedMath.CoulombF(0, -1.2, 0.0) == 0.0);
            Assert.IsTrue(AdvancedMath.CoulombF(5, 4.3, 0.0) == 0.0);

        }


        [TestMethod]
        public void CoulombEtaZero () {

            // For \eta = 0, Coulomb wave functions reduce to spherical Bessel functions

            foreach (int L in TestUtilities.GenerateIntegerValues(1, 100, 8)) {
                foreach (double rho in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 16)) {

                    double F = AdvancedMath.CoulombF(L, 0, rho);
                    double j = AdvancedMath.SphericalBesselJ(L, rho);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(F, rho * j, 16.0 * TestUtilities.TargetPrecision));

                    double G = AdvancedMath.CoulombG(L, 0, rho);
                    double y = AdvancedMath.SphericalBesselY(L, rho);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(G, -rho * y, 32.0 * TestUtilities.TargetPrecision));

                }
            }

        }

        [TestMethod]
        public void CoulombAgreement () {

            // The individual CoulombF and CoulombG functions should agree with the Coulomb
            // function that computes both.

            foreach (int L in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (double eta in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 8)) {
                    foreach (double rho in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 16)) {

                        double F = AdvancedMath.CoulombF(L, eta, rho);
                        double G = AdvancedMath.CoulombG(L, eta, rho);
                        SolutionPair s = AdvancedMath.Coulomb(L, eta, rho);

                        Assert.IsTrue(TestUtilities.IsNearlyEqual(F, s.FirstSolutionValue));

                        // We are producing NaNs for divergent G rathar than infinities
                        if (F == 0.0) continue;

                        Assert.IsTrue(TestUtilities.IsNearlyEqual(G, s.SecondSolutionValue));

                    }
                }
            }
        }

        [TestMethod]
        public void CoulombRecursion () {

            // Abromowitz & Stegun 14.2.1-14.2.3 recursion relations

            foreach (int L in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (double eta in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 4)) {
                    foreach (double rho in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 8)) {
                        CoulombRecursionTest(L, eta, rho);
                        CoulombRecursionTest(L, -eta, rho);
                    }
                }
            }

         }

        private void CoulombRecursionTest (int L, double eta, double rho) {

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
        public void CoulombWronskian () {

            foreach (int L in TestUtilities.GenerateIntegerValues(1, 100, 8)) {
                foreach (double eta in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 16)) {
                    foreach (double rho in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 32)) {

                        SolutionPair rp = AdvancedMath.Coulomb(L, eta, rho);
                        Debug.WriteLine("{0} {1} {2} : {3} {4} {5} {6} : {7}",
                            L, eta, rho,
                            rp.FirstSolutionValue, rp.FirstSolutionDerivative,
                            rp.SecondSolutionValue, rp.SecondSolutionDerivative,
                            rp.FirstSolutionDerivative * rp.SecondSolutionValue - rp.FirstSolutionValue * rp.SecondSolutionDerivative
                        );

                        // if F and G underflow and overflow, skip check
                        if (rp.FirstSolutionValue == 0.0) continue;

                        Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                            rp.FirstSolutionDerivative * rp.SecondSolutionValue, -rp.FirstSolutionValue * rp.SecondSolutionDerivative, 1.0,
                            TestUtilities.TargetPrecision * 10.0
                        ));


                        SolutionPair rm = AdvancedMath.Coulomb(L, -eta, rho);
                        Debug.WriteLine("{0} {1} {2} : {3} {4} {5} {6} : {7}",
                            L, -eta, rho,
                            rm.FirstSolutionValue, rm.FirstSolutionDerivative,
                            rm.SecondSolutionValue, rm.SecondSolutionDerivative,
                            rm.FirstSolutionDerivative * rm.SecondSolutionValue - rm.FirstSolutionValue * rm.SecondSolutionDerivative
                        );

                        Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                            rm.FirstSolutionDerivative * rm.SecondSolutionValue, -rm.FirstSolutionValue * rm.SecondSolutionDerivative, 1.0,
                            TestUtilities.TargetPrecision * 10.0
                        ));

                    }
                }
            }

        }
    }
}
