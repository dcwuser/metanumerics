using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;

using Meta.Numerics.Extended;
using Meta.Numerics.Statistics;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_Elliptic {

        [TestMethod]
        public void AccuracyTest () {

            SummaryStatistics s1 = new SummaryStatistics();
            SummaryStatistics s2 = new SummaryStatistics();

            foreach (double omx in TestUtilities.GenerateRealValues(1.0E-24, 1.0, 1000000)) {
                double x = 1.0 - omx;

                double f1 = Math.Sqrt((1.0 - x) * (1.0 + x));
                double f2 = Math.Sqrt(1.0 - x * x);

                DoubleDouble xe = (DoubleDouble) x;
                DoubleDouble fc = DoubleDouble.Sqrt(DoubleDouble.One - xe * xe);

                double e1 = Math.Abs((double) (f1 - fc));
                double e2 = Math.Abs((double) (f2 - fc));

                s1.Add(e1);
                s2.Add(e2);

            }

        }

        [TestMethod]
        public void EllipticKSpecialCases () {

            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticK(0.0), Math.PI / 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticK(1.0 / Math.Sqrt(2.0)), Math.Pow(AdvancedMath.Gamma(1.0 / 4.0), 2.0) / Math.Sqrt(Math.PI) / 4.0));
            Assert.IsTrue(Double.IsPositiveInfinity(AdvancedMath.EllipticK(1.0)));

        }

        [TestMethod]
        public void EllipticKExtremeValues () {
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticK(Double.Epsilon), Math.PI / 2.0));
            Assert.IsTrue(AdvancedMath.EllipticK(1.0 - 2.0E-16) < Double.PositiveInfinity);
            Assert.IsTrue(Double.IsNaN(AdvancedMath.EllipticK(Double.NaN)));
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
        public void EllipticKIntegrals () {

            Interval i = Interval.FromEndpoints(0.0, 1.0);

            // http://mathworld.wolfram.com/CatalansConstant.html equation 37
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate(
                    k => AdvancedMath.EllipticK(k),
                    Interval.FromEndpoints(0.0, 1.0)),
                2.0 * AdvancedMath.Catalan
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate(
                    k => AdvancedMath.EllipticK(k) * k,
                    Interval.FromEndpoints(0.0, 1.0)),
                1.0
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate(
                    k => AdvancedMath.EllipticK(k) / (1.0 + k),
                    Interval.FromEndpoints(0.0, 1.0)),
                Math.PI * Math.PI / 8.0
            ));

        }

        [TestMethod]
        public void CompleteEllipticHypergeometricAgreement () {

            foreach (double k in TestUtilities.GenerateRealValues(1.0E-4, 1.0, 4)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.EllipticK(k),
                    Math.PI / 2.0 * AdvancedMath.Hypergeometric2F1(0.5, 0.5, 1.0, k * k)
                ));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.EllipticE(k),
                    Math.PI / 2.0 * AdvancedMath.Hypergeometric2F1(0.5, -0.5, 1.0, k * k)
                ));
            }

        }

        [TestMethod]
        public void EllipticFSpecialCases () {
            foreach (double phi in TestUtilities.GenerateUniformRealValues(-10.0, 10.0, 8)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticF(phi, 0.0), phi));
            }
        }

        [TestMethod]
        public void EllipticFCompleteAgreement () {
            foreach (double k in TestUtilities.GenerateRealValues(0.01, 1.0, 8)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.EllipticF(Math.PI / 2.0, k),
                    AdvancedMath.EllipticK(k)
                ));
            }
        }

        [TestMethod]
        public void EllipticFDisplacement () {
            foreach (int m in TestUtilities.GenerateUniformIntegerValues(-100, 100, 4)) {
                foreach (double phi in TestUtilities.GenerateUniformRealValues(-Math.PI / 2.0, +Math.PI / 2.0, 4)) {
                    foreach (double k in TestUtilities.GenerateRealValues(1.0E-4, 1.0, 4)) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(
                            AdvancedMath.EllipticF(m * Math.PI + phi, k),
                            2 * m * AdvancedMath.EllipticK(k) + AdvancedMath.EllipticF(phi, k)
                        ));
                    }
                }
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
        public void EllipticESpecialCases () {
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticE(0.0), Math.PI / 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticE(1.0), 1.0));

            double g2 = Math.Pow(AdvancedMath.Gamma(1.0 / 4.0), 2.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.EllipticE(1.0 / Math.Sqrt(2.0)),
                Math.Pow(Math.PI, 3.0 / 2.0) / g2 + g2 / 8.0 / Math.Sqrt(Math.PI)
            ));
        }

        [TestMethod]
        public void EllipticEExtremeValues () {
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticE(Double.Epsilon), Math.PI / 2.0));
            Assert.IsTrue(Double.IsNaN(AdvancedMath.EllipticE(Double.NaN)));
        }

        [TestMethod]
        public void EllipticLegendreRelation () {
            foreach (double k in TestUtilities.GenerateRealValues(1.0E-4, 1.0, 8)) {
                double kp = Math.Sqrt(1.0 - k * k);
                double E = AdvancedMath.EllipticE(k);
                double EP = AdvancedMath.EllipticE(kp);
                double K = AdvancedMath.EllipticK(k);
                double KP = AdvancedMath.EllipticK(kp);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(K * EP + KP * E, Math.PI / 2.0 + K * KP));
            }
        }

        [TestMethod]
        public void EllipticECompleteAgreement () {
            foreach (double k in TestUtilities.GenerateRealValues(0.01, 1.0, 8)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.EllipticE(Math.PI / 2.0, k),
                    AdvancedMath.EllipticE(k)
                ));
            }
        }

        [TestMethod]
        public void EllipticEDisplacement () {
            foreach (int m in TestUtilities.GenerateUniformIntegerValues(-100, 100, 4)) {
                foreach (double phi in TestUtilities.GenerateUniformRealValues(-Math.PI / 2.0, +Math.PI / 2.0, 4)) {
                    foreach (double k in TestUtilities.GenerateRealValues(1.0E-4, 1.0, 4)) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(
                            AdvancedMath.EllipticE(m * Math.PI + phi, k),
                            2 * m * AdvancedMath.EllipticE(k) + AdvancedMath.EllipticE(phi, k)
                        ));
                    }
                }
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

        [TestMethod]
        public void EllipticPiSpecialCases () {
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticPi(0.0, 0.0), Math.PI / 2.0));
        }

        [TestMethod]
        public void EllipticPiZeroModulus () {
            // DLMF 19.6.3
            foreach (double n in TestUtilities.GenerateUniformRealValues(-4.0, 1.0, 8)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticPi(n, 0.0), Math.PI / 2.0 / Math.Sqrt(1.0 - n)));
            }
        }

        [TestMethod]
        public void EllipticPiSpecialCharacteristic () {
            foreach (double k in TestUtilities.GenerateRealValues(1.0E-4, 1.0, 4)) {
                // DLMF 19.6.1
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticPi(k * k, k), AdvancedMath.EllipticE(k) / (1.0 - k * k)));
                // DLMF 19.6.2
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticPi(-k, k), Math.PI / 4.0 / (1.0 + k) + AdvancedMath.EllipticK(k) / 2.0));
                // DLMF 19.6.3
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticPi(0.0, k), AdvancedMath.EllipticK(k)));
                Assert.IsTrue(Double.IsPositiveInfinity(AdvancedMath.EllipticPi(1.0, k)));
            }
        }

        // Defining integtral for Pi(n, k)
        [TestMethod]
        public void EllipticPiIntegration () {

            Interval i = Interval.FromEndpoints(0.0, Math.PI / 2.0);
            foreach (double k in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 4)) {
                double m = k * k;
                foreach (double n in TestUtilities.GenerateUniformRealValues(-2.0, 1.0, 4)) {
                    Func<double, double> f = delegate (double t) {
                        double s2 = MoreMath.Sqr(Math.Sin(t));
                        return (1.0 / (1.0 - n * s2) / Math.Sqrt(1.0 - m * s2));
                    };

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        FunctionMath.Integrate(f, i), AdvancedMath.EllipticPi(n, k)
                    ));
                }

            }

        }

        [TestMethod]
        public void EllipticNomeAgreement () {

            foreach (double k in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 4)) {

                double K = AdvancedMath.EllipticK(k);

                double k1 = Math.Sqrt((1.0 - k) * (1.0 + k));
                double K1 = AdvancedMath.EllipticK(k1);

                double q = AdvancedMath.EllipticNome(k);

                // For k << 1, this test fails if formulated as q = e^{\pi K' / K}.
                // Problem is that computation of k' looses some digits to cancellation,
                // then K' is near singularity so error is amplified, then
                // exp function amplifies error some more. Investigation indicates
                // that q agrees with Mathematica to nearly all digits, so problem
                // isn't in nome function. Run test in log space and don't let k get
                // extremely small.

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    Math.Log(q),
                    -Math.PI * K1 / K
                ));
            }

        }

    }
}
