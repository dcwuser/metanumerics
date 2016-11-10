﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_Elliptic {

        [TestMethod]
        public void EllipticKSpecialCases () {

            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticK(0.0), Math.PI / 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticK(1.0 / Math.Sqrt(2.0)), Math.Pow(AdvancedMath.Gamma(1.0 / 4.0), 2.0) / Math.Sqrt(Math.PI) / 4.0));
            Assert.IsTrue(Double.IsPositiveInfinity(AdvancedMath.EllipticK(1.0)));

        }

        [TestMethod]
        public void EllipticKInequality () {
            // DLMF 19.9.1
            foreach (double k in TestUtilities.GenerateUniformRealValues(0.0, 1.0, 4)) {
                double S = AdvancedMath.EllipticK(k) + Math.Log(Math.Sqrt(1.0 - k * k));
                Assert.IsTrue((Math.Log(4.0) <= S) && (S <= Math.PI / 2.0));
            }
        }

        [TestMethod]
        public void EllipticKIntegration () {

            // The defining integral for K

            Interval i = Interval.FromEndpoints(0.0, Math.PI / 2.0);

            foreach (double k in TestUtilities.GenerateRealValues(0.01, 1.0, 8)) {

                Func<double, double> f = delegate (double t) {
                    double z = k * Math.Sin(t);
                    return (1.0 / Math.Sqrt(1.0 - z * z));
                };

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    FunctionMath.Integrate(f, i), AdvancedMath.EllipticK(k)
                ));

            }

        }

        [TestMethod]
        public void EllipticKCatalanIntegral () {

            System.Diagnostics.Stopwatch timer = System.Diagnostics.Stopwatch.StartNew();


            Interval i = Interval.FromEndpoints(0.0, 1.0);

            // http://mathworld.wolfram.com/CatalansConstant.html equation 37

            Func<double, double> f1 = delegate (double k) {
                return (AdvancedMath.EllipticK(k));
            };
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate(f1, i), 2.0 * AdvancedMath.Catalan
            ));

            Func<double, double> f2 = delegate (double k) {
                return (AdvancedMath.EllipticK(k) * k);
            };
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate(f2, i), 1.0
            ));

            timer.Stop();
            Console.WriteLine(timer.ElapsedTicks);

        }

        [TestMethod]
        public void EllipticFSpecialCases () {
            foreach (double phi in TestUtilities.GenerateUniformRealValues(0.0, Math.PI / 2.0, 4)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticF(phi, 0.0), phi));
            }
        }

        [TestMethod]
        public void EllipticFIntegration () {

            foreach (double k in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 8)) {

                Func<double, double> f = delegate (double t) {
                    double z = k * Math.Sin(t);
                    return (1.0 / Math.Sqrt(1.0 - z * z));
                };

                foreach (double phi in TestUtilities.GenerateUniformRealValues(0.0, Math.PI / 2.0, 8)) {

                    Interval i = Interval.FromEndpoints(0.0, phi);

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        FunctionMath.Integrate(f, i), AdvancedMath.EllipticF(phi, k)
                    ));

                }

            }

        }

        [TestMethod]
        public void EllipticFBetaRelationship () {

            foreach (double phi in TestUtilities.GenerateUniformRealValues(0.0, Math.PI / 2.0, 8)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.EllipticF(phi, 1.0 / Math.Sqrt(2.0)),
                    AdvancedMath.Beta(1.0 / 2.0, 1.0 / 4.0, 1.0 - MoreMath.Pow(Math.Cos(phi), 4)) / Math.Sqrt(8.0)
                ));
            }

        }

        [TestMethod]
        public void EllipticFIntegral () {

            foreach (double k in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 8)) {
                foreach (double phi in TestUtilities.GenerateUniformRealValues(0.0, Math.PI / 2.0, 4)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        FunctionMath.Integrate(t => AdvancedMath.EllipticF(t, k) / Math.Sqrt(1.0 - MoreMath.Pow(k * Math.Sin(t), 2)), Interval.FromEndpoints(0.0, phi)),
                        MoreMath.Pow(AdvancedMath.EllipticF(phi, k), 2) / 2.0
                    ));
                }
            }

        }

        [TestMethod]
        public void EllipticESPecialCases () {
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticE(0.0), Math.PI / 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticE(1.0), 1.0));

            double g2 = Math.Pow(AdvancedMath.Gamma(1.0 / 4.0), 2.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.EllipticE(1.0 / Math.Sqrt(2.0)),
                Math.Pow(Math.PI, 3.0 / 2.0) / g2 + g2 / 8.0 / Math.Sqrt(Math.PI)
            ));

        }

        [TestMethod]
        public void EllipticLegendreRelation () {
            foreach (double k in TestUtilities.GenerateRealValues(0.01, 1.0, 8)) {
                double kp = Math.Sqrt(1.0 - k * k);
                double E = AdvancedMath.EllipticE(k);
                double EP = AdvancedMath.EllipticE(kp);
                double K = AdvancedMath.EllipticK(k);
                double KP = AdvancedMath.EllipticK(kp);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(K * EP + KP * E, Math.PI / 2.0 + K * KP));
            }
        }

        [TestMethod]
        public void EllipticEIntegration () {

            foreach (double k in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 8)) {

                // define the integrand
                Func<double, double> f = delegate (double t) {
                    double z = k * Math.Sin(t);
                    return (Math.Sqrt(1.0 - z * z));
                };

                // test complete integral
                double C = FunctionMath.Integrate(f, Interval.FromEndpoints(0.0, Math.PI / 2.0));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticE(k), C));

                // test incomplete integral
                foreach (double phi in TestUtilities.GenerateUniformRealValues(0.0, Math.PI / 2.0, 8)) {
                    double I = FunctionMath.Integrate(f, Interval.FromEndpoints(0.0, phi));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticE(phi, k), I));
                }

            }

        }

        [TestMethod]
        public void EllipticLandenTransform () {

            foreach (double k in TestUtilities.GenerateUniformRealValues(0.0, 1.0, 8)) {

                double kp = Math.Sqrt(1.0 - k * k);
                double k1 = (1.0 - kp) / (1.0 + kp);

                // For complete 1st kind
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.EllipticK(k),
                    (1.0 + k1) * AdvancedMath.EllipticK(k1)
                ));

                // For complete 2nd kind
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.EllipticE(k) + kp * AdvancedMath.EllipticK(k),
                    (1.0 + kp) * AdvancedMath.EllipticE(k1)
                ));

                /*
                foreach (double t in TestUtilities.GenerateUniformRealValues(0.0, Math.PI / 2.0, 4)) {

                    double s = Math.Sin(t);
                    double t1 = Math.Asin((1.0 + kp) * s * Math.Sqrt(1.0 - s * s) / Math.Sqrt(1.0 - k * k * s * s));

                    // Since t1 can be > pi/2, we need to handle a larger range of angles before we can test the incomplete cases
 

                }

                */
            }

        }

    }
}
