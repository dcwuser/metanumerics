using System;
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

    }
}
