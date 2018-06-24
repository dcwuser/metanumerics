using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_Exponential {

        [TestMethod]
        public void IntegralEAtZero () {
            // A & S 5.1.23
            foreach (int n in TestUtilities.GenerateIntegerValues(2, 100, 8)) {
                Assert.IsTrue(AdvancedMath.IntegralE(n, 0.0) == 1.0 / (n - 1.0));
            }
        }

        [TestMethod]
        public void IntegralE0 () {
            // A & S 5.1.24
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E3, 8)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.IntegralE(0, x), Math.Exp(-x) / x));
            }
        }

        [TestMethod]
        public void IntegralEDefinition () {
            // A & S 5.1.4
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 10.0, 8)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.IntegralE(n, x),
                        FunctionMath.Integrate(t => Math.Exp(-x * t) / MoreMath.Pow(t, n), 1.0, Double.PositiveInfinity)
                    ));
                }
            }
        }

        [TestMethod]
        public void IntegralERecurrence () {
            // A & S 5.1.14
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 24)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        n * AdvancedMath.IntegralE(n + 1, x) + x * AdvancedMath.IntegralE(n, x),
                        Math.Exp(-x)
                    ));
                }
            }
        }

        [TestMethod]
        public void IntegralE1Inequality () {
            // A & S 5.1.20
            // Choose maximum x to keep e^x from overflowing.
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, Math.Log(Double.MaxValue / 10.0), 8)) {
                double lower = Math.Log(1.0 + 2.0 / x) / 2.0;
                double upper = Math.Log(1.0 + 1.0 / x);
                double value = Math.Exp(x) * AdvancedMath.IntegralE(1, x);
                Assert.IsTrue((lower < value) && (value < upper));
            }
        }

        [TestMethod]
        public void IntegralEInequality () {
            // A & S 5.1.19
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 8)) {
                // Choose maximum x to keep e^x from overflowing.
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, Math.Log(Double.MaxValue / 10.0), 8)) {
                    double lower = 1.0 / (x + n);
                    double upper = 1.0 / (x + (n - 1));
                    double value = Math.Exp(x) * AdvancedMath.IntegralE(n, x);
                    Assert.IsTrue((lower < value) && (value <= upper));
                }
            }
        }

        [TestMethod]
        public void IntegralE1Integral () {
            // A & S 5.1.33
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate(t => MoreMath.Sqr(AdvancedMath.IntegralE(1, t)), 0.0, Double.PositiveInfinity),
                2.0 * Math.Log(2.0)
            ));
        }

        [TestMethod]
        public void IntegralSiSpecialCases () {
            Assert.IsTrue(AdvancedMath.IntegralSi(Double.NegativeInfinity) == -Math.PI / 2.0);
            Assert.IsTrue(AdvancedMath.IntegralSi(0.0) == 0.0);
            Assert.IsTrue(AdvancedMath.IntegralSi(Double.PositiveInfinity) == Math.PI / 2.0);
            Assert.IsTrue(Double.IsNaN(AdvancedMath.IntegralSi(Double.NaN)));
        }

        [TestMethod]
        public void IntegralSiReflection () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 8)) {
                Assert.IsTrue(AdvancedMath.IntegralSi(-x) == -AdvancedMath.IntegralSi(x));
            }
        }

        [TestMethod]
        public void IntegralSiDefinition () {
            Func<double, double> f = t => Math.Sin(t) / t;
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 8)) {
                Interval r = Interval.FromEndpoints(0.0, x);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.IntegralSi(x), FunctionMath.Integrate(f, r)));
            }
        }

        [TestMethod]
        public void IntegralCiSpecialCases () {
            Assert.IsTrue(Double.IsNegativeInfinity(AdvancedMath.IntegralCi(0.0)));
            Assert.IsTrue(Double.IsNaN(AdvancedMath.IntegralCi(Double.NaN)));
        }

        // Tried to also test Ci definition, but the integral is too close to divergent
        // for us to get good numerical accuracy.

        /*
        [TestMethod]
        public void IntegralCiDefinition () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E1, 8)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    FunctionMath.Integrate(t => MoreMath.Cos(t) / t, 0.0, x),
                    AdvancedMath.IntegralCi(x)
                ));
            }
        }
        */

        [TestMethod]
        public void IntegralCinSpecialCases () {
            Assert.IsTrue(AdvancedMath.IntegralCin(0.0) == 0.0);
        }

        [TestMethod]
        public void IntegralCinReflection () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 8)) {
                Assert.IsTrue(AdvancedMath.IntegralCin(-x) == AdvancedMath.IntegralCin(x));
            }
        }

        [TestMethod]
        public void IntegralCinDefinition () {
            // To avoid cancellation error, near cos(x) ~ 1, use half-angle formula
            //   1 - \cos(x) = 2 \sin^2(x/2)
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 8)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    FunctionMath.Integrate(t => 2.0 * MoreMath.Sqr(MoreMath.Sin(t / 2.0)) / t, 0.0, x),
                    AdvancedMath.IntegralCin(x)
                ));
            }
        }

        [TestMethod]
        public void IntegralCiCinRelation () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0, 4)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.EulerGamma + Math.Log(x) - AdvancedMath.IntegralCin(x),
                    AdvancedMath.IntegralCi(x)
                ));
            }
            foreach (double x in TestUtilities.GenerateRealValues(1.0, 1.0E4, 4)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.EulerGamma + Math.Log(x) - AdvancedMath.IntegralCi(x),
                    AdvancedMath.IntegralCin(x)
                ));
            }
        }

        [TestMethod]
        public void IntegralCiSiIntegrals () {

            // these integrals are from Oldham et al, An Atlas of Functions

            // A & S 5.2.28
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate(t => AdvancedMath.IntegralCi(t) * Math.Exp(-t), Interval.FromEndpoints(0.0, Double.PositiveInfinity)),
                -Math.Log(2.0) / 2.0
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate(t => AdvancedMath.IntegralSi(t) * Math.Exp(-t), Interval.FromEndpoints(0.0, Double.PositiveInfinity)),
                Math.PI / 4.0
            ));

            /*
            // the integral of [Ci(t)]^2 does not converge numerically 
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate(t => MoreMath.Sqr(AdvancedMath.IntegralCi(t)), Interval.FromEndpoints(0.0, Double.PositiveInfinity)),
                Math.PI / 2.0
            ));
            */
        }

        [TestMethod]
        public void IntegralShiSpecialCases () {
            Assert.IsTrue(AdvancedMath.IntegralShi(0.0) == 0.0);
        }

        [TestMethod]
        public void IntegralShiReflection () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 4)) {
                Assert.IsTrue(AdvancedMath.IntegralShi(-x) == -AdvancedMath.IntegralShi(x));
            }
        }

        [TestMethod]
        public void IntegralShiDefintion () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E1, 4)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    FunctionMath.Integrate(t => Math.Sinh(t) / t, 0.0, x),
                    AdvancedMath.IntegralShi(x)
                ));
            }
        }

    }
}
