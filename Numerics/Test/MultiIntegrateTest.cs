using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public class MultiIntegrateTest {

        public Interval[] UnitCube (int d) {
            Interval[] box = new Interval[d];
            for (int j = 0; j < d; j++) {
                box[j] = Interval.FromEndpoints(0.0, 1.0);
            }
            return (box);
        }

        [TestMethod]
        public void SeperableIntegrals () {

            Func<double[], double> f = delegate(double[] x) {
                double y = 1.0;
                for (int j = 0; j < x.Length; j++) {
                    y *= x[j];
                }
                return (y);
            };

            for (int d = 1; d <= 8; d++) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    FunctionMath.Integrate(f, UnitCube(d)),
                    MoreMath.Pow(2.0, -d),
                    1.0E-3 * Math.Pow(2.0, d / 2.0)
                ));
            }

        }

        [TestMethod]
        public void BallVolumeIntegrals () {

            // The volume of a d-sphere is \frac{\pi^{d/2}}{\Gamma(d/2 + 1)}
            // and the fraction in the first quadrant is 1/2^d of that
            
            // This is a simple test of the integration of a discontinuous function

            Func<double[], double> f = delegate (double[] x) {
                double r2 = 0.0;
                for (int j = 0; j < x.Length; j++) {
                    r2 += x[j] * x[j];
                }
                if (r2 <= 1.0) {
                    return(1.0);
                } else {
                    return(0.0);
                }
            };

            for (int d = 1; d <= 8; d++) {
                Console.WriteLine(d);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    FunctionMath.Integrate(f, UnitCube(d)),
                    Math.Pow(Math.PI, d / 2.0) / AdvancedMath.Gamma(d / 2.0 + 1.0) * MoreMath.Pow(2.0, -d),
                    1.0E-3 * Math.Pow(2.0, d / 2.0)
                ));

            }

        }

        [TestMethod]
        public void ZetaIntegrals () {

            // By expanding 1/(1-x) = 1 + x + x^2 + ... and integrating term-by-term
            // it's easy to show that this integral is \zeta(d)

            Func<double[], double> f = delegate (double[] x) {
                double p = 1.0;
                for (int i = 0; i < x.Length; i++) {
                    p *= x[i];
                }
                return (1.0 / (1.0 - p));
            };

            // \zeta(1) is infinite, so skip d=1

            for (int d = 2; d <= 8; d++) {
                Console.WriteLine(d);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    FunctionMath.Integrate(f, UnitCube(d)),
                    AdvancedMath.RiemannZeta(d),
                    1.0E-3 * Math.Pow(2.0, d / 2.0)
                ));
            }

        }

        [TestMethod]
        public void DoubleIntegrals () {

            // At http://mathworld.wolfram.com/DoubleIntegral.html, Mathworld
            // lists a few double integrals with known values that we take as
            // tests.
            
            // One of the integrals listed there is just the zeta integral for d=2
            // which we do above, so we omit it here.

            // Because these are non-ocsilatory and have relatively low
            // dimension, we can demand fairly high accuracy.

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate((double[] x) => 1.0 / (1.0 - x[0] * x[0] * x[1] * x[1]), UnitCube(2)),
                Math.PI * Math.PI / 8.0,
                5.0E-3
            ));
            /*
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate((double[] x) => 1.0 / (x[0] + x[1]) / Math.Sqrt((1.0 - x[0]) * (1.0 - x[1])), UnitCube(2)),
                4.0 * AdvancedMath.Catalan,
                5.0E-3
            ));
            */
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate((double[] x) => (x[0] - 1.0) / (1.0 - x[0] * x[1]) / Math.Log(x[0] * x[1]), UnitCube(2)),
                AdvancedMath.EulerGamma,
                5.0E-3
            ));
            

        }

        [TestMethod]
        public void WatsonIntegrals () {

            // Watson defined and analytically integrated three complicated tripple integrals related to random walks in three dimension
            // See http://mathworld.wolfram.com/WatsonsTripleIntegrals.html

            // These integrals are oscilatory, so up the budget to about 4,000,000 and reduce the target accuracy to about 1/2%
            EvaluationSettings settings = new EvaluationSettings() { RelativePrecision = 2.0E-3, EvaluationBudget = (1 << 23) };

            Interval watsonWidth = Interval.FromEndpoints(0.0, Math.PI);
            Interval[] watsonBox = new Interval[] { watsonWidth, watsonWidth, watsonWidth };
            
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate((double[] x) => 1.0 / (1.0 - Math.Cos(x[0]) * Math.Cos(x[1]) * Math.Cos(x[2])), watsonBox, settings),
                MoreMath.Pow(AdvancedMath.Gamma(1.0 / 4.0), 4) / 4.0,
                8.0E-3
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate((double[] x) => 1.0 / (3.0 - Math.Cos(x[0]) - Math.Cos(x[1]) - Math.Cos(x[2])), watsonBox, settings),
                Math.Sqrt(6.0) / 96.0 * AdvancedMath.Gamma(1.0 / 24.0) * AdvancedMath.Gamma(5.0 / 24.0) * AdvancedMath.Gamma(7.0 / 24.0) * AdvancedMath.Gamma(11.0 / 24.0),
                8.0E-3
            ));
            
        }

    }
}
