using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_SphericalBessel {

        [TestMethod]
        public void SphericalBesselSpecialCase () {
            Assert.IsTrue(AdvancedMath.SphericalBesselJ(0, 0.0) == 1.0);
            Assert.IsTrue(AdvancedMath.SphericalBesselJ(1, 0.0) == 0.0);
            Assert.IsTrue(AdvancedMath.SphericalBesselY(0, 0.0) == Double.NegativeInfinity);
            Assert.IsTrue(AdvancedMath.SphericalBesselY(1, 0.0) == Double.NegativeInfinity);
        }

        [TestMethod]
        public void SphericalBesselRecurrence () {
            // A&S 10.1.19
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E6, 16)) {

                    double jm1 = AdvancedMath.SphericalBesselJ(n - 1, x);
                    double jp1 = AdvancedMath.SphericalBesselJ(n + 1, x);
                    double j = AdvancedMath.SphericalBesselJ(n, x);
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(jm1, jp1, (2 * n + 1) * j / x));

                    double ym1 = AdvancedMath.SphericalBesselY(n - 1, x);
                    double yp1 = AdvancedMath.SphericalBesselY(n + 1, x);
                    double y = AdvancedMath.SphericalBesselY(n, x);
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(ym1, yp1, (2 * n + 1) * y / x));

                }
            }
        }

        [TestMethod]
        public void SphericalBesselNegativeOrder () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E6, 8)) {
                    int s = n % 2 == 0 ? -1 : +1;
                    Assert.IsTrue(AdvancedMath.SphericalBesselY(n, x) == s * AdvancedMath.SphericalBesselJ(-n - 1, x));
                    Assert.IsTrue(AdvancedMath.SphericalBesselY(-n, x) == s * AdvancedMath.SphericalBesselJ(n - 1, x));
                }
            }
        }

        [TestMethod]
        public void SphericalBesselModulus () {

            // A&S 10.1.28-10.1.30

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E6, 16)) {

                double x2 = x * x;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    MoreMath.Sqr(AdvancedMath.SphericalBesselJ(0, x)) + MoreMath.Sqr(AdvancedMath.SphericalBesselY(0, x)), 1.0 / x2
                ));

                double x4 = x2 * x2;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    MoreMath.Sqr(AdvancedMath.SphericalBesselJ(1, x)) + MoreMath.Sqr(AdvancedMath.SphericalBesselY(1, x)), 1.0 / x2 + 1.0 / x4
                ));

                double x6 = x4 * x2;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    MoreMath.Sqr(AdvancedMath.SphericalBesselJ(2, x)) + MoreMath.Sqr(AdvancedMath.SphericalBesselY(2, x)), 1.0 / x2 + 3.0 / x4 + 9.0 / x6
                ));

            }

        }

        [TestMethod]
        public void SphericalBesselCrossProduct () {
            // A&S 10.1.31
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E6, 8)) {
                    double jn = AdvancedMath.SphericalBesselJ(n, x);
                    double ym = AdvancedMath.SphericalBesselY(n - 1, x);
                    double jm = AdvancedMath.SphericalBesselJ(n - 1, x);
                    double yn = AdvancedMath.SphericalBesselY(n, x);
                    if (Double.IsInfinity(ym)) continue;
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        jn * ym, -jm * yn, 1.0 / (x * x)
                    ));
                }
            }
        }

        [TestMethod]
        public void SphericalBesselRealBesselAgreement () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 16)) {

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                       AdvancedMath.SphericalBesselJ(n, x),
                       Math.Sqrt(Math.PI / 2.0 / x) * AdvancedMath.BesselJ(n + 0.5, x)
                   ));

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.SphericalBesselY(n, x),
                        Math.Sqrt(Math.PI / 2.0 / x) * AdvancedMath.BesselY(n + 0.5, x)
                    ));
                }
            }
        }

        [TestMethod]
        public void SphericalBesselMomentIntegral () {
            // \int_0^{\infty} dx x^n j_{n+1}(x) = 1/(2n+1)!!
            // which follows from derivative formula
            // (https://math.stackexchange.com/questions/805379/integral-involving-the-spherical-bessel-function-of-the-first-kind-int-0)
            foreach (int n in TestUtilities.GenerateIntegerValues(4, 100, 4)) {
                IntegrationResult M = FunctionMath.Integrate(x => MoreMath.Pow(x, -n) * AdvancedMath.SphericalBesselJ(n + 1, x), 0.0, Double.PositiveInfinity);
                Assert.IsTrue(M.Estimate.ConfidenceInterval(0.95).ClosedContains(1.0 / AdvancedIntegerMath.DoubleFactorial(2 * n + 1)));
            }

        }

        [TestMethod]
        public void SphericalBesselTower () {

            // A&S 10.1.50-10.1.52
            // x shouldn't be too big or we will need too many terms before the series converges

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 4)) {

                double s1 = 0.0;
                double s2 = 0.0;
                double s3 = 0.0;

                int n = 0;
                while (true) {

                    double s1_old = s1;
                    double s2_old = s2;
                    double s3_old = s3;

                    double j = AdvancedMath.SphericalBesselJ(n, x);

                    double ds1 = j * j;
                    double ds2 = (2 * n + 1) * ds1;
                    double ds3 = (n % 2 == 0) ? ds2 : -ds2;

                    s1 += ds1;
                    s2 += ds2;
                    s3 += ds3;

                    if ((s1_old == s1) && (s2_old == s2) && (s3_old == s3)) break;

                    n++;
                    if (n > 1000) throw new NonconvergenceException();

                }

                Assert.IsTrue(TestUtilities.IsNearlyEqual(s1, AdvancedMath.IntegralSi(2.0 * x) / (2.0 * x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s2, 1.0));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s3, MoreMath.Sinc(2.0 * x)));

            }

        }
    }
}
