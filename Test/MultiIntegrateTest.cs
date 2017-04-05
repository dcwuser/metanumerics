using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Analysis;
using Meta.Numerics.Matrices;

namespace Test {

    [TestClass]
    public class MultiIntegrateTest {


        // Returns a d-dimensional unit cube [0,1]^d

        public Interval[] UnitCube (int d) {
            Interval[] box = new Interval[d];
            for (int j = 0; j < d; j++) {
                box[j] = Interval.FromEndpoints(0.0, 1.0);
            }
            return (box);
        }

        // Returns a d-dimensional symmetric unit cube [-1,+1]^d

        public Interval[] SymmetricUnitCube (int d) {
            Interval[] box = new Interval[d];
            for (int j = 0; j < d; j++) {
                box[j] = Interval.FromEndpoints(-1.0, 1.0);
            }
            return (box);
        }

        [TestMethod]
        public void SeperableIntegrals () {

            // Integrates \Pi_{j=0}^{d-1} \int_0^1 \! dx \, x_j^j = \Pi_{j=0}^{d-1} \frac{1}{j+1} = \frac{1}{d!}

            // This is a simple test of a seperable integral

            Func<IList<double>, double> f = delegate (IList<double> x) {
                double y = 1.0;
                for (int j = 0; j < x.Count; j++) {
                    y *= MoreMath.Pow(x[j], j);
                }
                return (y);
            };

            for (int d = 1; d <= 10; d++) {
                Console.WriteLine(d);
                // Result gets very small, so rely on relative rather than absolute precision.
                IntegrationSettings s = new IntegrationSettings() { AbsolutePrecision = 0, RelativePrecision = Math.Pow(10.0, -(6.0 - d / 2.0)) };
                IntegrationResult r = MultiFunctionMath.Integrate(f, UnitCube(d), s);
                Assert.IsTrue(r.Estimate.ConfidenceInterval(0.95).ClosedContains(1.0 / AdvancedIntegerMath.Factorial(d)));
            }

        }

        [TestMethod]
        public void BallVolumeIntegrals () {

            // The volume of a d-sphere is \frac{\pi^{d/2}}{\Gamma(d/2 + 1)}
            // and the fraction in the first quadrant is 1/2^d of that.

            // This is a simple test of the integration of a discontinuous function.

            Func<IList<double>, double> f = delegate (IList<double> x) {
                double r2 = 0.0;
                for (int j = 0; j < x.Count; j++) {
                    r2 += x[j] * x[j];
                }
                if (r2 <= 1.0) {
                    return (1.0);
                } else {
                    return (0.0);
                }
            };

            for (int d = 1; d <= 8; d++) {
                if (d == 6) continue; //. For d=6, integral returns after just ~300 evaluation with an underestimated error; look into this
                Console.WriteLine(d);
                IntegrationSettings settings = new IntegrationSettings() { AbsolutePrecision = 0.0, RelativePrecision = 1.0E-2 };
                IntegrationResult result = MultiFunctionMath.Integrate(f, UnitCube(d), settings);
                double V = Math.Pow(Math.PI, d / 2.0) / AdvancedMath.Gamma(d / 2.0 + 1.0) * MoreMath.Pow(2.0, -d);
                Assert.IsTrue(result.Estimate.ConfidenceInterval(0.95).ClosedContains(V));
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
            for (int d = 2; d <= 15; d++) {
                Console.WriteLine(d);
                IntegrationResult result = MultiFunctionMath.Integrate(f, UnitCube(d));
                Assert.IsTrue(result.Estimate.ConfidenceInterval(0.95).ClosedContains(AdvancedMath.RiemannZeta(d)));
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

            IntegrationResult i1 = MultiFunctionMath.Integrate(
                (IList<double> x) => 1.0 / (1.0 - x[0] * x[0] * x[1] * x[1]),
                UnitCube(2)
            );
            Assert.IsTrue(i1.Estimate.ConfidenceInterval(0.95).ClosedContains(Math.PI * Math.PI / 8.0));

            IntegrationResult i2 = MultiFunctionMath.Integrate(
                (IList<double> x) => 1.0 / (x[0] + x[1]) / Math.Sqrt((1.0 - x[0]) * (1.0 - x[1])),
                UnitCube(2),
                new IntegrationSettings() { RelativePrecision = 1.0E-6 }
            );
            Assert.IsTrue(i2.Estimate.ConfidenceInterval(0.95).ClosedContains(4.0 * AdvancedMath.Catalan));
            // For higher precision demands on this integral, we start getting Infinity +/- NaN for estimate and never terminate, look into this.

            IntegrationResult i3 = MultiFunctionMath.Integrate(
                (IList<double> x) => (x[0] - 1.0) / (1.0 - x[0] * x[1]) / Math.Log(x[0] * x[1]),
                UnitCube(2)
            );
            Assert.IsTrue(i3.Estimate.ConfidenceInterval(0.95).ClosedContains(AdvancedMath.EulerGamma));

        }

        [TestMethod]
        public void UnitSquareIntegrals () {

            // http://mathworld.wolfram.com/UnitSquareIntegral.html has a long list of integrals over the unit square.

            // Many are taken from Guillera and Sondow,
            // "Double Integrals and Infinite Products for Some Classical Constants Via Analytic Continuations of Lerch's Transcendent.",
            // 16 June 2005. (http://arxiv.org/abs/math.NT/0506319)

            Assert.IsTrue(
                MultiFunctionMath.Integrate(
                    (IList<double> x) => 1.0 / (1.0 - x[0] * x[1]),
                    UnitCube(2)
                ).Estimate.ConfidenceInterval(0.95).ClosedContains(
                    AdvancedMath.RiemannZeta(2.0)
                )
            );

            Assert.IsTrue(
                MultiFunctionMath.Integrate(
                    (IList<double> x) => -Math.Log(x[0] * x[1]) / (1.0 - x[0] * x[1]),
                    UnitCube(2)
                ).Estimate.ConfidenceInterval(0.95).ClosedContains(
                    2.0 * AdvancedMath.RiemannZeta(3.0)
                )
            );

            Assert.IsTrue(
                MultiFunctionMath.Integrate(
                    (IList<double> x) => (x[0] - 1.0) / (1.0 - x[0] * x[1]) / Math.Log(x[0] * x[1]),
                    UnitCube(2)
                ).Estimate.ConfidenceInterval(0.95).ClosedContains(
                    AdvancedMath.EulerGamma
                )
            );

            Assert.IsTrue(
                MultiFunctionMath.Integrate(
                    (IList<double> x) => (x[0] - 1.0) / (1.0 + x[0] * x[1]) / Math.Log(x[0] * x[1]),
                    UnitCube(2)
                ).Estimate.ConfidenceInterval(0.95).ClosedContains(
                    Math.Log(4.0 / Math.PI)
                )
            );

            Assert.IsTrue(
                MultiFunctionMath.Integrate(
                    (IList<double> x) => 1.0 / (1.0 + MoreMath.Sqr(x[0] * x[1])),
                    UnitCube(2)
                ).Estimate.ConfidenceInterval(0.95).ClosedContains(
                    AdvancedMath.Catalan
                )
            );

        }

        [TestMethod]
        public void WatsonIntegrals () {

            // Watson defined and analytically integrated three complicated triple integrals related to random walks in three dimension
            // See http://mathworld.wolfram.com/WatsonsTripleIntegrals.html

            Interval watsonWidth = Interval.FromEndpoints(0.0, Math.PI);
            Interval[] watsonBox = new Interval[] { watsonWidth, watsonWidth, watsonWidth };
            
            Assert.IsTrue(
                MultiFunctionMath.Integrate(
                    (IList<double> x) => 1.0 / (1.0 - Math.Cos(x[0]) * Math.Cos(x[1]) * Math.Cos(x[2])), watsonBox
                ).Estimate.ConfidenceInterval(0.99).ClosedContains(
                    MoreMath.Pow(AdvancedMath.Gamma(1.0 / 4.0), 4) / 4.0
                )
            );

            Assert.IsTrue(
                MultiFunctionMath.Integrate(
                    (IList<double> x) => 1.0 / (3.0 - Math.Cos(x[0]) * Math.Cos(x[1]) - Math.Cos(x[1]) * Math.Cos(x[2]) - Math.Cos(x[0]) * Math.Cos(x[2])), watsonBox
                ).Estimate.ConfidenceInterval(0.99).ClosedContains(
                    3.0 * MoreMath.Pow(AdvancedMath.Gamma(1.0 / 3.0), 6) / Math.Pow(2.0, 14.0 / 3.0) / Math.PI
                )
            );

            Assert.IsTrue(
                MultiFunctionMath.Integrate(
                    (IList<double> x) => 1.0 / (3.0 - Math.Cos(x[0]) - Math.Cos(x[1]) - Math.Cos(x[2])), watsonBox
                ).Estimate.ConfidenceInterval(0.99).ClosedContains(
                    Math.Sqrt(6.0) / 96.0 * AdvancedMath.Gamma(1.0 / 24.0) * AdvancedMath.Gamma(5.0 / 24.0) * AdvancedMath.Gamma(7.0 / 24.0) * AdvancedMath.Gamma(11.0 / 24.0)
                )
            );
            
        }

        [TestMethod]
        public void SteinmetzVolume () {

            // Steinmetz solid is intersection of unit cylinders along all axes. This is another hard-edged integral.
            // Analytic values are known for d=2-5.
            // http://www.math.illinois.edu/~hildebr/ugresearch/cylinder-spring2013report.pdf
            // http://www.math.uiuc.edu/~hildebr/igl/nvolumes-fall2012report.pdf

            IntegrationSettings settings = new IntegrationSettings() { RelativePrecision = 1.0E-2 };

            for (int d = 2; d <= 5; d++) {

                IntegrationResult v1 = MultiFunctionMath.Integrate((IList<double> x) => {
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

                Console.WriteLine("{0} {1} {2}", d, v1.Value, v2);
                Assert.IsTrue(v1.Estimate.ConfidenceInterval(0.99).ClosedContains(v2));
            }

        }

        [TestMethod]
        public void RambleIntegrals () {

            // W_{n}(s) = \int_{0}^{1} dx_1 \cdots dx_n \vert \sum_k e^{2 \pi i x_k} \vert^{s}
            // appears in problems of random walk in n dimensions. This is an oscilatory integral.
            // Analytic values are known for W_3(1) and W_3(-1) (https://www.carma.newcastle.edu.au/jon/walks.pdf).
            // Also W_2(1) = 4 / \pi and even arguments (http://www.carma.newcastle.edu.au/jon/walkstalk.pdf).
            // W_4(-1) can be related to a one-dimensional integral involving elliptic functions
            // (https://www.carma.newcastle.edu.au/jon/walks2.pdf and http://crd-legacy.lbl.gov/~dhbailey/dhbpapers/combat.pdf).

            // This is an oscilatory integral, so drop the required precision.
            IntegrationSettings settings2 = new IntegrationSettings() { RelativePrecision = 1.0E-4 };
            IntegrationSettings settings3 = new IntegrationSettings() { RelativePrecision = 1.0E-3 };
            IntegrationSettings settings4 = new IntegrationSettings() { RelativePrecision = 1.0E-2 };

            Assert.IsTrue(RambleIntegral(2, 1, settings3).Estimate.ConfidenceInterval(0.95).ClosedContains(4.0 / Math.PI));

            Assert.IsTrue(RambleIntegral(3, -1, settings3).Estimate.ConfidenceInterval(0.95).ClosedContains(
                3.0 / 16.0 * Math.Pow(2.0, 1.0 / 3.0) / MoreMath.Pow(Math.PI, 4) * MoreMath.Pow(AdvancedMath.Gamma(1.0 / 3.0), 6)
            ));

            Assert.IsTrue(RambleIntegral(3, 1, settings3).Estimate.ConfidenceInterval(0.95).ClosedContains(
                3.0 / 16.0 * Math.Pow(2.0, 1.0 / 3.0) / MoreMath.Pow(Math.PI, 4) * MoreMath.Pow(AdvancedMath.Gamma(1.0 / 3.0), 6) +
                27.0 / 4.0 * Math.Pow(2.0, 2.0 / 3.0) / MoreMath.Pow(Math.PI, 4) * MoreMath.Pow(AdvancedMath.Gamma(2.0 / 3.0), 6)
            ));

            Assert.IsTrue(RambleIntegral(4, -1, settings4).Estimate.ConfidenceInterval(0.95).ClosedContains(
                8.0 / MoreMath.Pow(Math.PI, 3) * FunctionMath.Integrate(k => MoreMath.Sqr(AdvancedMath.EllipticK(k)), Interval.FromEndpoints(0.0, 1.0))
            ));

        }

        private IntegrationResult RambleIntegral (int d, int s, IntegrationSettings settings) {
            return (MultiFunctionMath.Integrate((IList<double> x) => {
                Complex z = 0.0;
                for (int k = 0; k < d; k++) {
                    z += ComplexMath.Exp(2.0 * Math.PI * Complex.I * x[k]);
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
            // http://www.davidhbailey.com/dhbpapers/boxintegrals.pdf and http://www.davidhbailey.com/dhbpapers/bbbz-conmath.pdf

            // More results are in http://crd-legacy.lbl.gov/~dhbailey/dhbpapers/BoxII.pdf

            IntegrationSettings settings = new IntegrationSettings() { EvaluationBudget = 1000000, RelativePrecision = 1.0E-3 };

            // 2D integrals

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                BoxIntegralB(2, -1, settings), 2.0 * Math.Log(1.0 + Math.Sqrt(2.0)), settings.RelativePrecision * 2
            ));
            // Note 2 \ln(1 + \sqrt{2}) = \ln(3 + 2 \sqrt{2}) because (1 + \sqrt{2})^2 = 3 + 2 \sqrt{2}

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                BoxIntegralB(2, 1, settings), (Math.Sqrt(2.0) + Math.Log(Math.Sqrt(2.0) + 1.0)) / 3.0, settings.RelativePrecision * 2
            ));

            /*
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                BoxIntegralD(1, -1, new EvaluationSettings() { EvaluationBudget = 100000, RelativePrecision = 1.0E-2 }),
                (2.0 - 4.0 * Math.Sqrt(2.0)) / 3.0 + 4.0 * Math.Log(1.0 + Math.Sqrt(2.0)), settings.RelativePrecision * 2
            ));
            */

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                BoxIntegralD(1, 1, settings), 1.0 / 3.0, settings.RelativePrecision * 2
            ));

            // 3D integrals

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                BoxIntegralB(3, -1, settings), Math.Log(5.0 + 3.0 * Math.Sqrt(3.0)) - Math.Log(2.0) / 2.0 - Math.PI / 4.0, settings.RelativePrecision * 4
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                BoxIntegralB(3, 1, settings), Math.Log(2.0 + Math.Sqrt(3.0)) / 2.0 + Math.Sqrt(3.0) / 4.0 - Math.PI / 24.0, settings.RelativePrecision * 2
            ));

            // 4D integrals

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                BoxIntegralB(4, -2, settings), Math.PI * Math.Log(2.0 + Math.Sqrt(3.0)) - 2.0 * AdvancedMath.Catalan - Math.PI * Math.PI / 8.0, settings.RelativePrecision * 4
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                BoxIntegralD(2, 1, settings), (2.0 + Math.Sqrt(2.0) + 5.0 * Math.Log(1.0 + Math.Sqrt(2.0))) / 15.0, settings.RelativePrecision * 4
            ));

            // 6D integrals
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                BoxIntegralD(3, 1, settings), 
                4.0 / 105.0 + 17.0 * Math.Sqrt(2.0) / 105.0 - 2.0 * Math.Sqrt(3.0) / 35.0 + Math.Log(1.0 + Math.Sqrt(2.0)) / 5.0 + 2.0 * Math.Log(2.0 + Math.Sqrt(3.0)) / 5.0 - Math.PI / 15.0,
                settings.RelativePrecision * 2
            ));


        }

        [TestMethod]
        public void MoreBox () {
            // 8-D integral!
            IntegrationSettings settings = new IntegrationSettings() { RelativePrecision = 1.0E-3 };
            Console.WriteLine(BoxIntegralD(4, -1, settings));
            /*
            double K1 = FunctionMath.Integrate(Math.Asec)
            double alpha = Math.Asin(2.0 / 3.0 - 1.0 / 6.0 / Math.Sqrt(2.0));
            double d41 = 26.0 / 15.0 * AdvancedMath.Catalan - 380.0 / 6237.0 * Math.Sqrt(5.0) +
                568.0 / 3465.0 * Math.Sqrt(3.0) - 4.0 / 189.0 * Math.PI - 449.0 / 3465.0 -
                73.0 / 63.0 * Math.Sqrt(2.0) * Math.Atan(Math.Sqrt(2.0) / 4.0) -
                184.0 / 189.0 * Math.Log(2.0) + 64.0 / 189.0 * Math.Log(Math.Sqrt(5.0) + 1.0) +
                1.0 / 54.0 * Math.Log(1.0 + Math.Sqrt(2.0)) + 40.0 / 63.0 * Math.Log(Math.Sqrt(2.0) + Math.Sqrt(6.0)) -
                5.0 / 28.0 * Math.PI * Math.Log(1.0 + Math.Sqrt(2.0)) + 52.0 / 63.0 * Math.PI * Math.Log(2.0) +
                295.0 / 252.0 * Math.Log(3.0) + 4.0 / 315.0 * Math.PI * Math.PI + 3239.0 / 62370.0 * Math.Sqrt(2.0) -
                8.0 / 21.0 * Math.Sqrt(3.0) * Math.Atan(1.0 / Math.Sqrt(15.0)) -
                52.0 / 63.0 * Math.PI * Math.Log(Math.Sqrt(2.0) + Math.Sqrt(6.0)) +
                5.0 / 7.0 * alpha * Math.Log(1.0 + Math.Sqrt(2.0)) - 5.0 / 7.0 * AdvancedMath.Clausen(alpha)
                - 5.0 / 7.0 * AdvancedMath.Clausen(alpha + Math.PI / 2.0) + 52.0 / 63.0 * K1;
            */
        }

        
        [TestMethod]
        public void ExponentialBox () {

            // Mentioned in passing in http://crd-legacy.lbl.gov/~dhbailey/dhbpapers/BoxII.pdf
            // This is easy because the integrand is seperable and smooth.

            for (int d = 2; d <= 10; d++) {

                IntegrationResult result = MultiFunctionMath.Integrate((IList<double> r) => {
                    double s = 0.0;
                    foreach (double x in r) s += x * x;
                    return (Math.Exp(-s));
                }, UnitCube(d));

                Assert.IsTrue(result.Estimate.ConfidenceInterval(0.95).ClosedContains(
                    MoreMath.Pow(Math.Sqrt(Math.PI) / 2.0 * AdvancedMath.Erf(1.0), d)
                ));

            }

        }
        

        [TestMethod]
        public void IsingIntegrals () {

            // See http://www.davidhbailey.com/dhbpapers/ising.pdf.

            IntegrationSettings settings = new IntegrationSettings() { RelativePrecision = 1.0E-3, EvaluationBudget = 1000000 };

            int d = 4;
            Interval[] volume = new Interval[d];
            for (int i = 0; i < volume.Length; i++) {
                volume[i] = Interval.FromEndpoints(0.0, Double.PositiveInfinity);
            }

            IntegrationResult r = MultiFunctionMath.Integrate((IList<double> x) => {
                double p = 1.0;
                double q = 0.0;
                for (int i = 0; i < x.Count; i++) {
                    double u = x[i];
                    double v = 1.0 / u;
                    q += (u + v);
                    p *= v;
                }
                return (p / MoreMath.Sqr(q));
            }, volume, settings);
            double c = 4.0 / AdvancedIntegerMath.Factorial(d);
            Console.WriteLine("{0} {1}", c * r.Estimate, r.EvaluationCount);

            Console.WriteLine(7.0 * AdvancedMath.RiemannZeta(3.0) / 12.0);

            Assert.IsTrue((c * r.Estimate).ConfidenceInterval(0.99).ClosedContains(7.0 * AdvancedMath.RiemannZeta(3.0) / 12.0));

        }

        public double BoxIntegralB (int d, int r, IntegrationSettings settings) {
            return (MultiFunctionMath.Integrate((IList<double> x) => {
                double s = 0.0;
                for (int k = 0; k < d; k++) {
                    s += x[k] * x[k];
                }
                return (Math.Pow(s, r / 2.0));
            }, UnitCube(d), settings).Value);
        }

        public double BoxIntegralD (int d, int r, IntegrationSettings settings) {
            return (MultiFunctionMath.Integrate((IList<double> x) => {
                double s = 0.0;
                for (int k = 0; k < d; k++) {
                    double z = x[k] - x[k + d];
                    s += z * z;
                }
                return (Math.Pow(s, r / 2.0));
            }, UnitCube(2 * d), settings).Value);
        }

        [TestMethod]
        public void GaussianIntegrals () {

            Random rng = new Random(1);
            for (int d = 2; d < 8; d++) {

                Console.WriteLine(d);

                // Create a symmetric matrix
                SymmetricMatrix A = new SymmetricMatrix(d);
                for (int r = 0; r < d; r++) {
                    for (int c = 0; c < r; c++) {
                        A[r, c] = rng.NextDouble();
                    }
                    // Ensure it is positive definite by diagonal dominance
                    A[r, r] = r + 1.0;
                }

                // Compute its determinant, which appears in the analytic value of the integral
                CholeskyDecomposition CD = A.CholeskyDecomposition();
                double detA = CD.Determinant();

                // Compute the integral
                Func<IList<double>, double> f = (IList<double> x) => {
                    ColumnVector v = new ColumnVector(x);
                    double s = v.Transpose() * (A * v);
                    return (Math.Exp(-s));
                };

                Interval[] volume = new Interval[d];
                for (int i = 0; i < d; i++) volume[i] = Interval.FromEndpoints(Double.NegativeInfinity, Double.PositiveInfinity);

                // These are difficult integrals; demand reduced precision.
                IntegrationSettings settings = new IntegrationSettings() { RelativePrecision = Math.Pow(10.0, -(4.0 - d / 2.0)) };

                IntegrationResult I = MultiFunctionMath.Integrate(f, volume, settings);

                // Compare to the analytic result
                Assert.IsTrue(I.Estimate.ConfidenceInterval(0.95).ClosedContains(Math.Sqrt(MoreMath.Pow(Math.PI, d) / detA)));
            }

        }

    }
}
