using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;

using Meta.Numerics.Extended;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_Bessel {

        [TestMethod]
        public void IntegerBesselSpecialCase () {
            Assert.IsTrue(AdvancedMath.BesselJ(0, 0.0) == 1.0);
            Assert.IsTrue(AdvancedMath.BesselJ(1, 0.0) == 0.0);
            Assert.IsTrue(AdvancedMath.BesselY(0, 0.0) == Double.NegativeInfinity);
            Assert.IsTrue(AdvancedMath.BesselY(1, 0.0) == Double.NegativeInfinity);
        }

        [TestMethod]
        public void IntegerBesselNegativeOrder () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E6, 8)) {
                    int s = (n % 2 == 0) ? +1 : -1;
                    Assert.IsTrue(AdvancedMath.BesselJ(-n, x) == s * AdvancedMath.BesselJ(n, x));
                    Assert.IsTrue(AdvancedMath.BesselY(-n, x) == s * AdvancedMath.BesselY(n, x));
                }
            }
        }

        [TestMethod]
        public void IntegerBesselNegativeArgument () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E6, 8)) {
                    int s = (n % 2 == 0) ? +1 : -1;
                    Assert.IsTrue(AdvancedMath.BesselJ(n, -x) == s * AdvancedMath.BesselJ(n, x));
                }
            }
        }

        [TestMethod]
        public void IntegerBesselRecurrence () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E6, 16)) {
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        AdvancedMath.BesselJ(n - 1, x), AdvancedMath.BesselJ(n + 1, x), 2 * n / x * AdvancedMath.BesselJ(n, x)
                    ));
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        AdvancedMath.BesselY(n - 1, x), AdvancedMath.BesselY(n + 1, x), 2 * n / x * AdvancedMath.BesselY(n, x)
                    ));
                }
            }
        }

        [TestMethod]
        public void IntegerBesselCrossProduct () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E6, 8)) {
                    double Jn = AdvancedMath.BesselJ(n, x);
                    double Jp = AdvancedMath.BesselJ(n + 1, x);
                    double Yn = AdvancedMath.BesselY(n, x);
                    double Yp = AdvancedMath.BesselY(n + 1, x);
                    if (Double.IsInfinity(Yn)) continue;
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        Jp * Yn, -Jn * Yp, 2.0 / Math.PI / x
                    ));
                }
            }
        }

        [TestMethod]
        public void BesselJ0Integral () {
            // Abromowitz & Stegun 9.1.18
            // J_0(x) = \int_{0}^{\pi} \cos( x \sin(t) ) \, dt
            // don't let x get too big, or the integral becomes too oscillatory to do accurately
            Interval r = Interval.FromEndpoints(0.0, Math.PI);
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 8)) {
                Func<double, double> f = delegate (double t) {
                    return (Math.Cos(x * Math.Sin(t)));
                };
                double J0 = FunctionMath.Integrate(f, r).Estimate.Value / Math.PI;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.BesselJ(0, x), J0));
            }
        }

        [TestMethod]
        public void BesselY0Integral () {
            // Abromowitz & Stegen 9.1.19
            Interval r = Interval.FromEndpoints(0.0, Math.PI / 2.0);
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 8)) {
                Func<double, double> f = delegate (double t) {
                    double s = Math.Sin(t);
                    return (Math.Cos(x * Math.Cos(t)) * (AdvancedMath.EulerGamma + Math.Log(2.0 * x * s * s)));
                };
                double Y0 = 4.0 / (Math.PI * Math.PI) * FunctionMath.Integrate(f, r).Estimate.Value;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.BesselY(0, x), Y0));
            }
        }

        [TestMethod]
        public void IntegerBesselJIntegral () {
            // Abromowitz & Stegun 9.1.21
            IntegrationSettings settings = new IntegrationSettings() { AbsolutePrecision = 2.5E-15, RelativePrecision = 0.0 };
            foreach (double x in TestUtilities.GenerateRealValues(0.1, 10.0, 8)) {
                foreach (int n in TestUtilities.GenerateIntegerValues(1, 10, 4)) {
                    double J = FunctionMath.Integrate(
                        t => MoreMath.Cos(x * MoreMath.Sin(t) - n * t),
                        0.0, Math.PI, settings
                    ).Estimate.Value / Math.PI;
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.BesselJ(n, x), J, settings
                    ), $"n={n}, x={x}, J={J}");
                    // The integral can produce significant cancelation, so use an absolute criterion.
                }
            }
        }

        [TestMethod]
        public void BesselKapteynIntegral () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 8)) {
                Func<double, double> f = delegate (double t) {
                    return (Math.Cos(x - t) * AdvancedMath.BesselJ(0, t));
                };
                Interval r = Interval.FromEndpoints(0.0, x);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(x * AdvancedMath.BesselJ(0, x), FunctionMath.Integrate(f, r).Estimate.Value));
            }
        }

        [TestMethod]
        public void BesselLipshitzIntegral () {
            // \int_{0}^{\infty} e^{-x} J_0(x) \, dx = \frac{1}{\sqrt{2}}
            Func<double, double> f = delegate (double t) {
                return (Math.Exp(-t) * AdvancedMath.BesselJ(0, t));
            };
            Interval r = Interval.FromEndpoints(0.0, Double.PositiveInfinity);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(FunctionMath.Integrate(f, r).Estimate.Value, 1.0 / Math.Sqrt(2.0)));
        }

        [TestMethod]
        public void BesselWeberIntegral () {
            Func<double, double> f = delegate (double t) {
                return (Math.Exp(-t * t) * AdvancedMath.BesselJ(0, t) * t);
            };
            Interval r = Interval.FromEndpoints(0.0, Double.PositiveInfinity);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(FunctionMath.Integrate(f, r).Estimate.Value, Math.Exp(-1.0 / 4.0) / 2.0));
        }

        


        [TestMethod]
        public void BesselAtZero () {

            // Normal Bessel \nu = 0

            Assert.IsTrue(AdvancedMath.BesselJ(0, 0.0) == 1.0);
            Assert.IsTrue(AdvancedMath.BesselJ(0.0, 0.0) == 1.0);

            Assert.IsTrue(Double.IsNegativeInfinity(AdvancedMath.BesselY(0, 0.0)));
            Assert.IsTrue(Double.IsNegativeInfinity(AdvancedMath.BesselY(0.0, 0.0)));

            SolutionPair jy0 = AdvancedMath.Bessel(0.0, 0.0);
            Assert.IsTrue(jy0.FirstSolutionValue == 1.0);
            Assert.IsTrue(jy0.FirstSolutionDerivative == 0.0);
            Assert.IsTrue(Double.IsNegativeInfinity(jy0.SecondSolutionValue));
            Assert.IsTrue(Double.IsPositiveInfinity(jy0.SecondSolutionDerivative));

            // Normal Bessel 0 < \nu < 1

            Assert.IsTrue(AdvancedMath.BesselJ(0.1, 0.0) == 0.0);

            Assert.IsTrue(Double.IsNegativeInfinity(AdvancedMath.BesselY(0.1, 0.0)));

            SolutionPair jyf = AdvancedMath.Bessel(0.9, 0.0);
            Assert.IsTrue(jyf.FirstSolutionValue == 0.0);
            Assert.IsTrue(jyf.FirstSolutionDerivative == Double.PositiveInfinity);
            Assert.IsTrue(Double.IsNegativeInfinity(jyf.SecondSolutionValue));
            Assert.IsTrue(Double.IsPositiveInfinity(jyf.SecondSolutionDerivative));

            // Normal Bessel \nu = 1

            Assert.IsTrue(AdvancedMath.BesselJ(1, 0.0) == 0.0);
            Assert.IsTrue(AdvancedMath.BesselJ(1.0, 0.0) == 0.0);

            Assert.IsTrue(Double.IsNegativeInfinity(AdvancedMath.BesselY(1, 0.0)));
            Assert.IsTrue(Double.IsNegativeInfinity(AdvancedMath.BesselY(1.0, 0.0)));

            SolutionPair jy1 = AdvancedMath.Bessel(1.0, 0.0);
            Assert.IsTrue(jy1.FirstSolutionValue == 0.0);
            Assert.IsTrue(jy1.FirstSolutionDerivative == 0.5);
            Assert.IsTrue(Double.IsNegativeInfinity(jy1.SecondSolutionValue));
            Assert.IsTrue(Double.IsPositiveInfinity(jy1.SecondSolutionDerivative));

            // Normal Bessel \nu > 1

            Assert.IsTrue(AdvancedMath.BesselJ(2, 0.0) == 0.0);
            Assert.IsTrue(AdvancedMath.BesselJ(1.2, 0.0) == 0.0);

            Assert.IsTrue(Double.IsNegativeInfinity(AdvancedMath.BesselY(2, 0.0)));
            Assert.IsTrue(Double.IsNegativeInfinity(AdvancedMath.BesselY(1.2, 0.0)));

            SolutionPair jy2 = AdvancedMath.Bessel(1.7, 0.0);
            Assert.IsTrue(jy2.FirstSolutionValue == 0.0);
            Assert.IsTrue(jy2.FirstSolutionDerivative == 0.0);
            Assert.IsTrue(Double.IsNegativeInfinity(jy2.SecondSolutionValue));
            Assert.IsTrue(Double.IsPositiveInfinity(jy2.SecondSolutionDerivative));

        }

        [TestMethod]
        public void ModifiedBesselAtZero () {

            // Modified Bessel \nu = 0

            Assert.IsTrue(AdvancedMath.ModifiedBesselI(0.0, 0.0) == 1.0);

            Assert.IsTrue(AdvancedMath.ModifiedBesselK(0.0, 0.0) == Double.PositiveInfinity); 

            SolutionPair ik0 = AdvancedMath.ModifiedBessel(0.0, 0.0);
            Assert.IsTrue(ik0.FirstSolutionValue == 1.0);
            Assert.IsTrue(ik0.FirstSolutionDerivative == 0.0);
            Assert.IsTrue(Double.IsNegativeInfinity(ik0.SecondSolutionValue));
            Assert.IsTrue(Double.IsPositiveInfinity(ik0.SecondSolutionDerivative));

            // Modified Bessel 0 < \nu < 1

            Assert.IsTrue(AdvancedMath.ModifiedBesselI(0.1, 0.0) == 0.0);

            Assert.IsTrue(AdvancedMath.ModifiedBesselK(0.1, 0.0) == Double.PositiveInfinity);

            SolutionPair ikf = AdvancedMath.ModifiedBessel(0.9, 0.0);
            Assert.IsTrue(ikf.FirstSolutionValue == 0.0);
            Assert.IsTrue(ikf.FirstSolutionDerivative == Double.PositiveInfinity);
            Assert.IsTrue(Double.IsNegativeInfinity(ikf.SecondSolutionValue));
            Assert.IsTrue(Double.IsPositiveInfinity(ikf.SecondSolutionDerivative));

            // Modified Bessel \nu = 1

            Assert.IsTrue(AdvancedMath.ModifiedBesselI(1.0, 0.0) == 0.0);

            Assert.IsTrue(AdvancedMath.ModifiedBesselK(1.0, 0.0) == Double.PositiveInfinity);

            SolutionPair ik1 = AdvancedMath.ModifiedBessel(1.0, 0.0);
            Assert.IsTrue(ik1.FirstSolutionValue == 0.0);
            Assert.IsTrue(ik1.FirstSolutionDerivative == 0.5);
            Assert.IsTrue(Double.IsNegativeInfinity(ik1.SecondSolutionValue));
            Assert.IsTrue(Double.IsPositiveInfinity(ik1.SecondSolutionDerivative));

            // Modified Bessel \nu > 1

            Assert.IsTrue(AdvancedMath.ModifiedBesselI(1.2, 0.0) == 0.0);

            Assert.IsTrue(AdvancedMath.ModifiedBesselK(1.2, 0.0) == Double.PositiveInfinity);

            SolutionPair ik2 = AdvancedMath.ModifiedBessel(1.7, 0.0);
            Assert.IsTrue(ik2.FirstSolutionValue == 0.0);
            Assert.IsTrue(ik2.FirstSolutionDerivative == 0.0);
            Assert.IsTrue(Double.IsNegativeInfinity(ik2.SecondSolutionValue));
            Assert.IsTrue(Double.IsPositiveInfinity(ik2.SecondSolutionDerivative));

        }

    }

}