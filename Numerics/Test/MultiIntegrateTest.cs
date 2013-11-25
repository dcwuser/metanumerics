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

        public Interval[] SymmetricUnitCube (int d) {
            Interval[] box = new Interval[d];
            for (int j = 0; j < d; j++) {
                box[j] = Interval.FromEndpoints(-1.0, 1.0);
            }
            return (box);
        }

        [TestMethod]
        public void SeperableIntegrals () {

            Func<IList<double>, double> f = delegate(IList<double> x) {
                double y = 1.0;
                for (int j = 0; j < x.Count; j++) {
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

            Func<IList<double>, double> f = delegate (IList<double> x) {
                double r2 = 0.0;
                for (int j = 0; j < x.Count; j++) {
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

            Func<IList<double>, double> f = delegate (IList<double> x) {
                double p = 1.0;
                for (int i = 0; i < x.Count; i++) {
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
                FunctionMath.Integrate((IList<double> x) => 1.0 / (1.0 - x[0] * x[0] * x[1] * x[1]), UnitCube(2)),
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
                FunctionMath.Integrate((IList<double> x) => (x[0] - 1.0) / (1.0 - x[0] * x[1]) / Math.Log(x[0] * x[1]), UnitCube(2)),
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
                FunctionMath.Integrate((IList<double> x) => 1.0 / (1.0 - Math.Cos(x[0]) * Math.Cos(x[1]) * Math.Cos(x[2])), watsonBox, settings),
                MoreMath.Pow(AdvancedMath.Gamma(1.0 / 4.0), 4) / 4.0,
                8.0E-3
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate((IList<double> x) => 1.0 / (3.0 - Math.Cos(x[0]) - Math.Cos(x[1]) - Math.Cos(x[2])), watsonBox, settings),
                Math.Sqrt(6.0) / 96.0 * AdvancedMath.Gamma(1.0 / 24.0) * AdvancedMath.Gamma(5.0 / 24.0) * AdvancedMath.Gamma(7.0 / 24.0) * AdvancedMath.Gamma(11.0 / 24.0),
                8.0E-3
            ));
            
        }

        [TestMethod]
        public void SteinmetzVolume () {

            // Steinmetz solid is intersection of unit cylinders along all axes. This is another hard-edged integral. Analytic values are known for d=2-5.
            // http://www.math.illinois.edu/~hildebr/ugresearch/cylinder-spring2013report.pdf
            // http://www.math.uiuc.edu/~hildebr/igl/nvolumes-fall2012report.pdf

            EvaluationSettings settings = new EvaluationSettings() { EvaluationBudget = 10000000, RelativePrecision = 1.0E-3 };

            for (int d = 2; d <= 5; d++) {

                double v1 = FunctionMath.Integrate((IList<double> x) => {
                    for (int i = 0; i < d; i++) {
                        double s = 0.0;
                        for (int j = 0; j < d; j++) {
                            if (j != i) s += x[j] * x[j];
                        }
                        if (s > 1.0) return (0.0);
                    }
                    return (1.0);
                }, SymmetricUnitCube(d), settings);

                double v2 = 0.0;
                switch (d) {
                    case 2:
                        // trivial square
                        v2 = 4.0;
                        break;
                    case 3:
                        v2 = 16.0 - 8.0 * Math.Sqrt(2.0);
                        break;
                    case 4:
                        v2 = 48.0 * (Math.PI / 4.0 - Math.Atan(Math.Sqrt(2.0)) / Math.Sqrt(2.0));
                        break;
                    case 5:
                        v2 = 256.0 * (Math.PI / 12.0 - Math.Atan(1.0 / (2.0 * Math.Sqrt(2.0))) / Math.Sqrt(2.0));
                        break;
                }

                Console.WriteLine("{0} {1} {2}", d, v1, v2);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(v1, v2, settings.RelativePrecision));

            }

        }

        [TestMethod]
        public void RambleIntegral () {

            // W_{n}(s) = \int_{0}^{1} dx_1 \cdots dx_n \vert \sum_k e^{2 \pi i x_k} \vert^{s}
            // appears in problems of random walk in n dimensions. This is an oscilatory integral. Analytic values are known for W_3(1) and W_3(-1).

            EvaluationSettings settings = new EvaluationSettings() { EvaluationBudget = 10000000, RelativePrecision = 1.0E-3 };

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                RambleIntegral(3, -1, settings),
                3.0 / 16.0 * Math.Pow(2.0, 1.0 / 3.0) / MoreMath.Pow(Math.PI, 4) * MoreMath.Pow(AdvancedMath.Gamma(1.0 / 3.0), 6),
                settings.RelativePrecision
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                RambleIntegral(3, 1, settings),
                3.0 / 16.0 * Math.Pow(2.0, 1.0 / 3.0) / MoreMath.Pow(Math.PI, 4) * MoreMath.Pow(AdvancedMath.Gamma(1.0 / 3.0), 6) +
                27.0 / 4.0 * Math.Pow(2.0, 2.0 / 3.0) / MoreMath.Pow(Math.PI, 4) * MoreMath.Pow(AdvancedMath.Gamma(2.0 / 3.0), 6),
                settings.RelativePrecision
            ));

        }

        private double RambleIntegral (int d, int s, EvaluationSettings settings) {
            return (FunctionMath.Integrate((IList<double> x) => {
                Complex z = 0.0;
                for (int k = 0; k < d; k++) {
                    z += ComplexMath.Exp(2.0 * Math.PI * ComplexMath.I * x[k]);
                }
                return (MoreMath.Pow(ComplexMath.Abs(z), s));
            }, UnitCube(d), settings));
        }


        [TestMethod]
        public void BoxIntegrals () {

            // The box integrals
            //   B_n(r) = \int{0}^{1} dx_1 \cdots dx_n ( x_1^2 + \cdots x_n^2 )^{r/2}
            //   D_n(r) = \int{0}^{1} dx_1 \cdots dx_n dy_1 \cdots dy_n \left[ (x_1 - y_1)^2 + \cdots (x_n - y_n)^2 \right]^{r/2}
            // Give the mean distance of a point in a box from its center or from other points. Various of these are known analytically.
            // see Bailey, Borwein, Crandall, "Box Integrals", Journal of Computational and Applied Mathematics 206 (2007) 196
            // http://www.davidhbailey.com/dhbpapers/boxintegrals.pdf

            EvaluationSettings settings = new EvaluationSettings() { EvaluationBudget = 10000000, RelativePrecision = 1.0E-3 };

            // 2D integrals

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                BoxIntegralB(2, -1, settings), Math.Log(3.0 + 2.0 * Math.Sqrt(2.0)), settings.RelativePrecision * 2
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                BoxIntegralB(2, 1, settings), (Math.Sqrt(2.0) + Math.Log(Math.Sqrt(2.0) + 1.0)) / 3.0, settings.RelativePrecision * 2
            ));

            //Assert.IsTrue(TestUtilities.IsNearlyEqual(
            //    BoxIntegralD(1, -1, settings), (2.0 - 4.0 * Math.Sqrt(2.0)) / 3.0 + 4.0 * Math.Log(1.0 + Math.Sqrt(2.0)) , settings.RelativePrecision * 2
            //));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                BoxIntegralD(1, 1, settings), 1.0 / 3.0, settings.RelativePrecision * 2
            ));

            // 3D integrals

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                BoxIntegralB(3, -1, settings), Math.Log(5.0 + 3.0 * Math.Sqrt(3.0)) - Math.Log(2.0) / 2.0 - Math.PI / 4.0, settings.RelativePrecision * 2
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                BoxIntegralB(3, 1, settings), Math.Log(2.0 + Math.Sqrt(3.0)) / 2.0 + Math.Sqrt(3.0) / 4.0 - Math.PI / 24.0, settings.RelativePrecision * 2
            ));

            // 4D integrals

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                BoxIntegralD(2, 1, settings), (2.0 + Math.Sqrt(2.0) + 5.0 * Math.Log(1.0 + Math.Sqrt(2.0))) / 15.0, settings.RelativePrecision * 2
            ));

            // 6D integrals

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                BoxIntegralD(3, 1, settings), 
                4.0 / 105.0 + 17.0 * Math.Sqrt(2.0) / 105.0 - 2.0 * Math.Sqrt(3.0) / 35.0 + Math.Log(1.0 + Math.Sqrt(2.0)) / 5.0 + 2.0 * Math.Log(2.0 + Math.Sqrt(3.0)) / 5.0 - Math.PI / 15.0,
                settings.RelativePrecision * 2
            ));

        }

        public double BoxIntegralB (int d, int r, EvaluationSettings settings) {
            return (FunctionMath.Integrate((IList<double> x) => {
                double s = 0.0;
                for (int k = 0; k < d; k++) {
                    s += x[k] * x[k];
                }
                return (Math.Pow(s, r / 2.0));
            }, UnitCube(d), settings));
        }

        public double BoxIntegralD (int d, int r, EvaluationSettings settings) {
            return (FunctionMath.Integrate((IList<double> x) => {
                double s = 0.0;
                for (int k = 0; k < d; k++) {
                    double z = x[k] - x[k + d];
                    s += z * z;
                }
                return (Math.Pow(s, r / 2.0));
            }, UnitCube(2 * d), settings));
        }

    }
}
