using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_Jacobi {

        [TestMethod]
        public void JacobiSpecial () {

            // Verify values at u = 0, K/2, K given in DMLF 22.5i and A&S

            double k = 1.0 / 3.1;
            double k1 = Math.Sqrt(1.0 - k * k);
            double K = AdvancedMath.EllipticK(k);

            // u = 0
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.JacobiSn(0.0, k), 0.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.JacobiCn(0.0, k), 1.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.JacobiDn(0.0, k), 1.0));

            // u = K / 2 
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.JacobiSn(K / 2.0, k), 1.0 / Math.Sqrt(1.0 + k1)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.JacobiCn(K / 2.0, k), Math.Sqrt(k1 / (1.0 + k1))));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.JacobiDn(K / 2.0, k), Math.Sqrt(k1)));

            // u = K
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.JacobiSn(K, k), 1.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.JacobiCn(K, k), 0.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.JacobiDn(K, k), k1));

        }

        [TestMethod]
        public void JacobiReflection () {

            double k = 0.31415;
            foreach (double u in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 3)) {

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.JacobiSn(u, k),
                    -AdvancedMath.JacobiSn(-u, k)
                ));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.JacobiCn(u, k),
                    AdvancedMath.JacobiCn(-u, k)
                ));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.JacobiDn(u, k),
                    AdvancedMath.JacobiDn(-u, k)
                ));

            }

        }

        [TestMethod]
        public void JacobiPeriodic () {


            double u = 3.1;
            double k = 0.2;

            double K = AdvancedMath.EllipticK(k);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.JacobiSn(u, k),
                -AdvancedMath.JacobiSn(u + 2.0 * K, k) 
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.JacobiSn(u, k),
                AdvancedMath.JacobiSn(u + 4.0 * K, k)
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.JacobiCn(u, k),
                -AdvancedMath.JacobiCn(u + 2.0 * K, k)
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.JacobiCn(u, k),
                AdvancedMath.JacobiCn(u + 4.0 * K, k)
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.JacobiDn(u, k),
                AdvancedMath.JacobiDn(u + 2.0 * K, k)
            ));

        }

        [TestMethod]
        public void JacobiIdentities () {

            // Verify domain and simple identities

            foreach (double u in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 8)) {
                foreach (double k in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 4)) {

                    double sn = AdvancedMath.JacobiSn(u, k);
                    double cn = AdvancedMath.JacobiCn(u, k);
                    double dn = AdvancedMath.JacobiDn(u, k);

                    Assert.IsTrue(Math.Abs(sn) <= 1.0);
                    Assert.IsTrue(Math.Abs(cn) <= 1.0);
                    Assert.IsTrue((0.0 <= dn) && (dn <= 1.0));

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(sn * sn + cn * cn, 1.0));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(MoreMath.Sqr(k * sn) + dn * dn, 1.0));

                }
            }

        }


        [TestMethod]
        public void JacobiLimits () {

            // Verify limits at m=0, 1 given in DLMF 22.5ii and A&S

            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.JacobiSn(1.1, 0.0), Math.Sin(1.1)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.JacobiCn(2.2, 0.0), Math.Cos(2.2)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.JacobiDn(3.3, 0.0), 1.0));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.JacobiSn(3.3, 1.0), Math.Tanh(3.3)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.JacobiCn(4.4, 1.0), 1.0 / Math.Cosh(4.4)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.JacobiDn(5.5, 1.0), 1.0 / Math.Cosh(5.5)));

        }


        [TestMethod]
        public void JacobiDoubleAndHalf () {

            foreach (double k in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 3)) {
                foreach (double u in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 6)) {

                    double sn = AdvancedMath.JacobiSn(u, k);
                    double cn = AdvancedMath.JacobiCn(u, k);
                    double dn = AdvancedMath.JacobiDn(u, k);

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.JacobiSn(2.0 * u, k),
                        2.0 * sn * cn * dn / (1.0 - k * k * MoreMath.Pow(sn, 4))
                    ));

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.JacobiCn(2.0 * u, k),
                        (cn * cn - sn * sn * dn * dn) / (1.0 - k * k * MoreMath.Pow(sn, 4))
                    ));

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        MoreMath.Sqr(AdvancedMath.JacobiSn(u / 2.0, k)),
                        (1.0 - cn) / (1.0 + dn)
                    ));

                }
            }

        }

        [TestMethod]
        public void JacobiIntegrals () {

            foreach (double k in TestUtilities.GenerateRealValues(1.0E-1, 1.0, 4)) {

                // DLMF 22.14.18 says
                //   \int_{0}^{K(k)} \ln(\dn(t, k)) \, dt = \frac{1}{2} K(k) \ln k'

                double K = AdvancedMath.EllipticK(k);

                double k1 = Math.Sqrt(1.0 - k * k);

                double I1 = FunctionMath.Integrate(
                    x => Math.Log(AdvancedMath.JacobiDn(x, k)),
                    Interval.FromEndpoints(0.0, K)
                );

                Assert.IsTrue(TestUtilities.IsNearlyEqual(I1, K / 2.0 * Math.Log(k1)));

                // If k is small, log(k1) is subject to cancellation error, so don't pick k too small.

                // It also gives values for the analogous integrals with \sn and \cn,
                // but their values involve K'(k), for which we don't have a method.

                // Mathworld (http://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html) says
                //   \int_{0}^{K(k)} \dn^2(t, k) = E(k)

                double I2 = FunctionMath.Integrate(
                    u => MoreMath.Sqr(AdvancedMath.JacobiDn(u, k)),
                    Interval.FromEndpoints(0.0, K)
                );

                Assert.IsTrue(TestUtilities.IsNearlyEqual(I2, AdvancedMath.EllipticE(k)));

            }

        }

         
        [TestMethod]
        public void JacobiLandenTransformation () {

            // DMLF 22.7 covers Landen transformations, which change both u and k.

            double u = 3.2;
            double k = 0.54;

            double sn = AdvancedMath.JacobiSn(u, k);
            double cn = AdvancedMath.JacobiCn(u, k);
            double dn = AdvancedMath.JacobiDn(u, k);

            // Descending (ks < k)

            double k1 = Math.Sqrt(1.0 - k * k);
            double ks = (1.0 - k1) / (1.0 + k1);
            double us = u / (1.0 + ks);

            double sns = AdvancedMath.JacobiSn(us, ks);
            double cns = AdvancedMath.JacobiCn(us, ks);
            double dns = AdvancedMath.JacobiDn(us, ks);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                sn,
                (1.0 + ks) * sns / (1.0 + ks * sns * sns)
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                cn,
                cns * dns / (1.0 + ks * sns * sns)
            ));

            // Ascending (kt > k)

            double kt1 = (1.0 - k) / (1.0 + k);
            double kt = Math.Sqrt(1.0 - kt1 * kt1);
            double ut = u / (1.0 + kt1);

            double snt = AdvancedMath.JacobiSn(ut, kt);
            double cnt = AdvancedMath.JacobiCn(ut, kt);
            double dnt = AdvancedMath.JacobiDn(ut, kt);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                sn,
                (1.0 + kt1) * snt * cnt / dnt
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                cn,
                (1.0 + kt1) * (dnt * dnt - kt1) / (kt * kt * dnt)
            ));

        }

        [TestMethod]
        public void JacobiAsInverse () {

            foreach (double k in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 3)) {
                foreach (double u in TestUtilities.GenerateRealValues(1.0E-1, Math.PI, 3)) {

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.EllipticF(Math.Asin(AdvancedMath.JacobiSn(u, k)), k), u
                    ));

                }
            }

        }

    }
}
