using System;
using System.Collections.Generic;
using System.Diagnostics;

using TestClassAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestClassAttribute;
using TestMethodAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestMethodAttribute;
using ExpectedExceptionAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.ExpectedExceptionAttribute;
using Assert = Microsoft.VisualStudio.TestTools.UnitTesting.Assert;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace Test {

    [TestClass]
    public class OdeTest {

        [TestMethod]
        public void TestMultivariateRegression () {

            // Collect r^2 values from multivariate linear regressions.

            double cz = 1.0;
            double cx = 0.0;
            double cy = 0.0;

            Random rng = new Random(1001110000);
            ContinuousDistribution xDistribution = new UniformDistribution(Interval.FromEndpoints(-4.0, 8.0));
            ContinuousDistribution yDistribution = new UniformDistribution(Interval.FromEndpoints(-8.0, 4.0));
            ContinuousDistribution eDistribution = new NormalDistribution();

            List<double> r2Sample = new List<double>();

            for (int i = 0; i < 500; i++) {

                MultivariateSample xyzSample = new MultivariateSample(3);
                for (int k = 0; k < 12; k++) {
                    double x = xDistribution.GetRandomValue(rng);
                    double y = yDistribution.GetRandomValue(rng);
                    double z = cx * x + cy * y + cz + eDistribution.GetRandomValue(rng);
                    xyzSample.Add(x, y, z);
                }
                MultiLinearRegressionResult fit = xyzSample.LinearRegression(2);
                double fcx = fit.Parameters.ValuesVector[0];
                double fcy = fit.Parameters.ValuesVector[1];
                double fcz = fit.Parameters.ValuesVector[2];

                r2Sample.Add(fit.RSquared);
            }

            // r^2 values should be distributed as expected.
            ContinuousDistribution r2Distribution = new BetaDistribution((3 - 1) / 2.0, (12 - 3) / 2.0);

            TestResult ks = r2Sample.KolmogorovSmirnovTest(r2Distribution);
            Assert.IsTrue(ks.Probability > 0.05);
        }


        [TestMethod]
        public void TestBivariateRegression () {

            // Do a bunch of linear regressions. r^2 should be distributed as expected.

            double a0 = 1.0;
            double b0 = 0.0;

            Random rng = new Random(1001110000);
            ContinuousDistribution xDistribution = new UniformDistribution(Interval.FromEndpoints(-2.0, 4.0));
            ContinuousDistribution eDistribution = new NormalDistribution();

            List<double> r2Sample = new List<double>();

            for (int i = 0; i < 500; i++) {

                BivariateSample xySample = new BivariateSample();
                for (int k = 0; k < 10; k++) {
                    double x = xDistribution.GetRandomValue(rng);
                    double y = a0 + b0 * x + eDistribution.GetRandomValue(rng);
                    xySample.Add(x, y);
                }
                LinearRegressionResult fit = xySample.LinearRegression();
                double a = fit.Intercept.Value;
                double b = fit.Slope.Value;

                r2Sample.Add(fit.RSquared);
            }

            ContinuousDistribution r2Distribution = new BetaDistribution((2 - 1) / 2.0, (10 - 2) / 2.0);
            TestResult ks = r2Sample.KolmogorovSmirnovTest(r2Distribution);
            Assert.IsTrue(ks.Probability > 0.05);
        }

        [TestMethod]
        public void OdeExponential () {

            // The exponential function y = e^x satisfies
            //    y' = y
            // This is perhaps the simplest differential equation.

            Func<double, double, double> f = (double x, double y) => y;

            OdeResult r = FunctionMath.IntegrateOde(f, 0.0, 1.0, 2.0);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(r.Y, MoreMath.Sqr(Math.E)));

            Console.WriteLine(r.EvaluationCount);
        }

        [TestMethod]
        public void OdeNonlinear () {

            // y = \frac{y_0}{1 - y_0 (x - x_0)}
            Func<double, double, double> f = (double x, double y) => MoreMath.Sqr(y);

            int count = 0;
            OdeSettings settings = new OdeSettings() {
                RelativePrecision = 1.0E-8,
                EvaluationBudget = 1024,
                Listener = (OdeResult) => count++
            };
            OdeResult result = FunctionMath.IntegrateOde(f, 0.0, 1.0, 0.99, settings);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Y, 1.0 / (1.0 - 1.0 * (0.99 - 0.0)), result.Settings));

            Assert.IsTrue(count > 0);

            Console.WriteLine(result.EvaluationCount);
        }

        [TestMethod]
        public void OdeLogistic () {

            // y = \frac{y_0}{y_0 + (1 - y_0) e^{-(x - x_0)}
            Func<double, double, double> rhs = (double x, double y) => y * (1.0 - y);
            Func<double, double, double> solution = (double y0, double x) => y0 / (y0 + (1.0 - y0) * Math.Exp(-x));

            int count = 0;
            foreach (double y0 in new double[] { -0.1, 0.0, 0.4, 1.0, 1.6 }) {
                Console.WriteLine(y0);
                Interval r = Interval.FromEndpoints(0.0, 2.0);
                OdeResult s = FunctionMath.IntegrateOde(rhs, 0.0, y0, 2.0);
                double y1 = s.Y;
                Console.WriteLine(y1);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(y1, solution(y0, 2.0)));
                count += s.EvaluationCount;
            }
            Console.WriteLine(count);
        }

        [TestMethod]
        public void OdeExample () {

            Func<double, double, double> rhs = (double t, double u) => (1.0 - 2.0 * t) * u;
            Func<double, double, double> solution = (double u0, double t) => u0 * Math.Exp(t - t * t);

            int count = 0;
            foreach (double t in new double[] { 0.5, 0.75, 1.50, 2.25, 3.25 }) {
                OdeResult r = FunctionMath.IntegrateOde(rhs, 0.0, 1.0, t);
                double y1 = r.Y;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(y1, solution(1.0, t)));
                count += r.EvaluationCount;
            }
            Console.WriteLine(count);
        }

        [TestMethod]
        public void OdeSine () {

            // The sine and cosine functions satisfy
            //   y'' = - y
            // This is perhaps the simplest conservative differential equation.
            // (i.e. right hand side depends only on y, not y')

            Func<double, double, double> f = (double x, double y) => -y;

            int count = 0;
            OdeSettings settings = new OdeSettings() {
                Listener = (OdeResult r) => {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        MoreMath.Sqr(r.Y) + MoreMath.Sqr(r.YPrime), 1.0, r.Settings
                    ));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        r.Y, MoreMath.Sin(r.X), r.Settings
                    ));
                    count++;
                }
            };
            OdeResult result = FunctionMath.IntegrateConservativeOde(f, 0.0, 0.0, 1.0, 5.0, settings);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Y, MoreMath.Sin(5.0)));
            Assert.IsTrue(count > 0);
        }

        private static double EulerKineticEnergy (IList<double> I, IList<double> w) {
            return (I[0] * MoreMath.Sqr(w[0]) + I[1] * MoreMath.Sqr(w[1]) + I[2] * MoreMath.Sqr(w[2]));
        }

        private static double EulerAngularMomentum (IList<double> I, IList<double> w) {
            return (MoreMath.Sqr(I[0] * w[0]) + MoreMath.Sqr(I[1] * w[1]) + MoreMath.Sqr(I[2] * w[2]));

        }

        [TestMethod]
        public void OdeEuler () {

            // Euler's equations for rigid body motion (without external forces) are
            //   I_1 \dot{\omega}_1  = (I_2 - I_3) \omega_2 \omega_3
            //   I_2 \dot{\omega}_2  = (I_3 - I_1) \omega_3 \omega_1
            //   I_3 \dot{\omega}_3  = (I_1 - I_2) \omega_1 \omega_2
            // where \omega's are rotations about the principal axes and I's are the moments
            // of inertia about those axes. The rotational kinetic energy
            //   I_1 \omega_1^2 + I_2 \omega_2^2 + I_3 \omega_3^2 = 2 E
            // and angular momentum
            //   I_1^2 \omega_1^2 + I_2^2 \omega_2^2 + I_3^2 \omega_3^2 = L^2
            // are conserved quantities.

            // Landau & Lifshitz, Mechanics (3rd edition) solves these equations.

            // Note
            //   (\cn u)' = - (\sn u)(\dn u)
            //   (\sn u)' = (\dn u)(\cn u)
            //   (\dn u)' = -k^2 (\cn u)(\sn u)
            // So a solution with \omega_1 ~ \cn, \omega_2 ~ \sn, \omega_3 ~ \dn
            // would seem to work, if we can re-scale variables to eliminate factors.
            // In fact, this turns out to be the case.
            
            // Label axis so that I_3 > I_2 > I_1. If L^2 > 2 E I_2, define
            //   v^2 = \frac{(I_3 - I_2)(L^2 - 2E I_1)}{I_1 I_2 I_3}
            //   k^2 = \frac{(I_2 - I_1)(2E I_3 - L^2)}{(I_3 - I_2)(L^2 - 2E I_1)}
            //   A_1^2 = \frac{2E I_3 - L^2}{I_1 (I_3 - I_1)}
            //   A_2^2 = \frac{2E I_3 - L^2}{I_2 (I_3 - I_2)}
            //   A_3^2 = \frac{L^2 - 2E I_1}{I_3 (I_3 - I_1)}
            // (If L^2 < 2 E I_2, just switch subscripts 1 <-> 3). Then
            //   \omega_1 = A_1 \cn (c t, k)
            //   \omega_2 = A_2 \sn (c t, k)
            //   \omega_3 = A_3 \dn (c t, k)
            // Period is complete when v T = 4 K.

            ColumnVector I = new ColumnVector(3.0, 4.0, 5.0);

            Func<double, IReadOnlyList<double>, IReadOnlyList<double>> rhs = (double t, IReadOnlyList<double> w) => {
                return (new ColumnVector((I[1] - I[2]) * w[1] * w[2] / I[0], (I[2] - I[0]) * w[2] * w[0] / I[1], (I[0] - I[1]) * w[0] * w[1] / I[2]));
            };

            ColumnVector w0 = new ColumnVector(6.0, 0.0, 1.0);

            // Determine L^2 and 2 E
            double L2 = EulerAngularMomentum(I, w0);
            double E2 = EulerKineticEnergy(I, w0);

            // As Landau points out, these inequalities should be ensured by the definitions of E and L.
            Debug.Assert(E2 * I[0] < L2);
            Debug.Assert(L2 < E2 * I[2]);

            double v, k;
            Func<double, ColumnVector> sln;
            if (L2 > E2 * I[1]) {
                v = Math.Sqrt((I[2] - I[1]) * (L2 - E2 * I[0]) / I[0] / I[1] / I[2]);
                k = Math.Sqrt((I[1] - I[0]) * (E2 * I[2] - L2) / (I[2] - I[1]) / (L2 - E2 * I[0]));
                ColumnVector A = new ColumnVector(
                    Math.Sqrt((E2 * I[2] - L2) / I[0] / (I[2] - I[0])),
                    Math.Sqrt((E2 * I[2] - L2) / I[1] / (I[2] - I[1])),
                    Math.Sqrt((L2 - E2 * I[0]) / I[2] / (I[2] - I[0]))
                );
                sln = t => new ColumnVector(
                    A[0] * AdvancedMath.JacobiCn(v * t, k),
                    A[1] * AdvancedMath.JacobiSn(v * t, k),
                    A[2] * AdvancedMath.JacobiDn(v * t, k)
                );

            } else {
                v = Math.Sqrt((I[0] - I[1]) * (L2 - E2 * I[2]) / I[0] / I[1] / I[2]);
                k = Math.Sqrt((I[1] - I[2]) * (E2 * I[0] - L2) / (I[0] - I[1]) / (L2 - E2 * I[2]));
                ColumnVector A = new ColumnVector(
                    Math.Sqrt((L2 - E2 * I[2]) / I[0] / (I[0] - I[2])),
                    Math.Sqrt((E2 * I[0] - L2) / I[1] / (I[0] - I[1])),
                    Math.Sqrt((E2 * I[0] - L2) / I[2] / (I[0] - I[2]))
                );
                sln = t => new ColumnVector(
                    A[0] * AdvancedMath.JacobiDn(v * t, k),
                    A[1] * AdvancedMath.JacobiSn(v * t, k),
                    A[2] * AdvancedMath.JacobiCn(v * t, k)
                );
            }

            Debug.Assert(k < 1.0);

            double T = 4.0 * AdvancedMath.EllipticK(k) / v;

            int listenerCount = 0;
            MultiOdeSettings settings = new MultiOdeSettings() {
                Listener = (MultiOdeResult q) => {
                    listenerCount++;

                    // Verify that energy and angular momentum conservation is respected
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(EulerKineticEnergy(I, q.Y), E2, q.Settings));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(EulerAngularMomentum(I, q.Y), L2, q.Settings));

                    // Verify that the result agrees with the analytic solution
                    //ColumnVector YP = new ColumnVector(
                    //    A[0] * AdvancedMath.JacobiCn(v * q.X, k),
                    //    A[1] * AdvancedMath.JacobiSn(v * q.X, k),
                    //    A[2] * AdvancedMath.JacobiDn(v * q.X, k)
                    //);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(q.Y, sln(q.X), q.Settings));
                }
            };
            MultiOdeResult result = MultiFunctionMath.IntegrateOde(rhs, 0.0, w0, T, settings);
            Console.WriteLine(result.EvaluationCount);
            //MultiOdeResult r2 = NewMethods.IntegrateOde(rhs, 0.0, w0, T, settings);

            // Verify that one period of evolution has brought us back to our initial state
            Assert.IsTrue(TestUtilities.IsNearlyEqual(result.X, T));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Y, w0, settings));
            Assert.IsTrue(listenerCount > 0);

        }

        [TestMethod]
        public void OdeOrbit2 () {

            Func<double, IReadOnlyList<double>, IReadOnlyList<double>> rhs = (double t, IReadOnlyList<double> r) => {
                double d = MoreMath.Hypot(r[0], r[2]);
                double d3 = MoreMath.Pow(d, 3);
                return (new ColumnVector(r[1], -r[0] / d3, r[3], -r[2] / d3));
            };

            foreach (double e in TestUtilities.GenerateUniformRealValues(0.0, 1.0, 8)) {
                Console.WriteLine(e);

                ColumnVector r0 = new ColumnVector(1.0 - e, 0.0, 0.0, Math.Sqrt((1.0 + e) / (1.0 - e)));

                ColumnVector r1 = MultiFunctionMath.IntegrateOde(rhs, 0.0, r0, 2.0 * Math.PI).Y;

                Assert.IsTrue(TestUtilities.IsNearlyEqual(r0, r1, new EvaluationSettings() { RelativePrecision = 1.0E-9 }));
            }

        }


        private double OrbitAngularMomentum (IList<double> r, IList<double> rDot) {
            return (r[0] * rDot[1] - r[1] * rDot[0]);
        }

        private double OrbitEnergy (IList<double> r, IList<double> rDot) {
            return (MoreMath.Sqr(rDot[0]) + MoreMath.Sqr(rDot[1]) - 2.0 / MoreMath.Hypot(r[0], r[1]));
        }

        [TestMethod]
        public void OdeOrbitKepler () {

            // This is a simple Keplerian orbit.
            // Hull (1972) constructed initial conditions that guarantee a given orbital eccentricity
            // and a period of 2 \pi with unit masses and unit gravitational constant.

            Func<double, IReadOnlyList<double>, IReadOnlyList<double>> rhs = (double t, IReadOnlyList<double> r) => {
                double d = MoreMath.Hypot(r[0], r[1]);
                double d3 = MoreMath.Pow(d, 3);
                return (new double[] { -r[0] / d3, -r[1] / d3 });
            };

            foreach (double e in new double[] { 0.0, 0.1, 0.3, 0.5, 0.7, 0.9 }) {

                ColumnVector r0 = new ColumnVector(1.0 - e, 0.0);
                ColumnVector rDot0 = new ColumnVector(0.0, Math.Sqrt((1.0 + e) / (1.0 - e)));

                double E = OrbitEnergy(r0, rDot0);
                double L = OrbitAngularMomentum(r0, rDot0);

                MultiOdeSettings settings = new MultiOdeSettings() {
                    Listener = (MultiOdeResult b) => {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(OrbitAngularMomentum(b.Y, b.YPrime), L, b.Settings));
                        EvaluationSettings relaxed = new EvaluationSettings() { RelativePrecision = 2.0 * b.Settings.RelativePrecision, AbsolutePrecision = 2.0 * b.Settings.AbsolutePrecision };
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(OrbitEnergy(b.Y, b.YPrime), E, relaxed));
                    }
                };
                
               
                MultiOdeResult result = MultiFunctionMath.IntegrateConservativeOde(rhs, 0.0, r0, rDot0, 2.0 * Math.PI, settings);
                ColumnVector r1 = result.Y;
                
                Console.WriteLine(result.EvaluationCount);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(r0, r1, new EvaluationSettings() { RelativePrecision = 512 * result.Settings.RelativePrecision }));

                // For large eccentricities, we loose precision. This is apparently typical behavior.
                // Would be nice if we could quantify expected loss, or correct to do better.

            }
        }

        [TestMethod]
        public void OdeOrbitFigureEight () {

            Func<double, IReadOnlyList<double>, IReadOnlyList<double>> rhs = (double t, IReadOnlyList<double> q) => {

                double x01 = q[0] - q[2];
                double y01 = q[1] - q[3];

                double x02 = q[0] - q[4];
                double y02 = q[1] - q[5];

                double x12 = q[2] - q[4];
                double y12 = q[3] - q[5];

                double r01 = MoreMath.Pow(MoreMath.Hypot(x01, y01), 3);
                double r02 = MoreMath.Pow(MoreMath.Hypot(x02, y02), 3);
                double r12 = MoreMath.Pow(MoreMath.Hypot(x12, y12), 3);

                return (new double[] {
                    - x01 / r01 - x02 / r02,
                    - y01 / r01 - y02 / r02,
                    + x01 / r01 - x12 / r12,
                    + y01 / r01 - y12 / r12,
                    + x02 / r02 + x12 / r12,
                    + y02 / r02 + y12 / r12
                });

            };

            double[] q1 = new double[] { 0.97000435669734, -0.24308753153583 };
            double[] p3 = new double[] { -0.93240737144104, -0.86473146092102 };

            ColumnVector q0 = new ColumnVector( q1[0], q1[1], -q1[0], -q1[1], 0.0, 0.0 );
            ColumnVector p0 = new ColumnVector( -p3[0] / 2.0, -p3[1] / 2.0, -p3[0] / 2.0, -p3[1] / 2.0, p3[0], p3[1] );

            double T = 6.32591398292621;

            MultiOdeResult result = MultiFunctionMath.IntegrateConservativeOde(rhs, 0.0, q0, p0, T);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(q0, result.Y, 1.0E-10));

        }

        [TestMethod]
        public void OdeOrbitTriangle () {

            // Lagrange described three-body systems where each body sits at the vertex
            // of an equilateral triangle. Each body executes an elliptical orbit
            // and the triangle remains equilateral, although its size can change.

            // The case of equal masses on a circle is easily analytically tractable.
            // In this case the orbits are all along the same circle and the 
            // the triangle doesn't change size, it just rotates. 

            Func<double, IReadOnlyList<double>, IReadOnlyList<double>> rhs = (double t, IReadOnlyList<double> r) => {

                double x01 = r[0] - r[1];
                double x02 = r[0] - r[2];
                double x12 = r[1] - r[2];

                double y01 = r[3] - r[4];
                double y02 = r[3] - r[5];
                double y12 = r[4] - r[5];

                double r01 = MoreMath.Pow(MoreMath.Hypot(x01, y01), 3);
                double r02 = MoreMath.Pow(MoreMath.Hypot(x02, y02), 3);
                double r12 = MoreMath.Pow(MoreMath.Hypot(x12, y12), 3);

                return (new double[] {
                    - x01 / r01 - x02 / r02,
                    x01 / r01 - x12 / r12,
                    x02 / r02 + x12 / r12,
                    - y01 / r01 - y02 / r02,
                    y01 / r01 - y12 / r12,
                    y02 / r02 + y12 / r12
                });

            };

            double L = 1.0;
            double R = L / Math.Sqrt(3.0);
            double v = Math.Sqrt(1.0 / L);
            double T = 2.0 * Math.PI * Math.Sqrt(L * L * L / 3.0);

            double c = 1.0 / 2.0;
            double s = Math.Sqrt(3.0) / 2.0;

            ColumnVector r0 = new ColumnVector(R, -c * R, -c * R, 0.0, s * R, -s * R );
            double[] rp0 = new double[] { 0.0, -s * v, s * v, v, -c * v, -c * v };

            int count = 0;
            MultiOdeSettings settings = new MultiOdeSettings() {
                RelativePrecision = 1.0E-12,
                AbsolutePrecision = 1.0E-12,
                Listener = (MultiOdeResult mer) => {

                    count++;
                    ColumnVector r = mer.Y;

                    double x01 = r[0] - r[1];
                    double x02 = r[0] - r[2];
                    double x12 = r[1] - r[2];

                    double y01 = r[3] - r[4];
                    double y02 = r[3] - r[5];
                    double y12 = r[4] - r[5];

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(MoreMath.Hypot(x01, y01), L, mer.Settings));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(MoreMath.Hypot(x02, y02), L, mer.Settings));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(MoreMath.Hypot(x12, y12), L, mer.Settings));

                }
            };

            MultiOdeResult result = MultiFunctionMath.IntegrateConservativeOde(rhs, 0.0, r0, rp0, T, settings);

            // Test that one period brought us back to the same place.
            Assert.IsTrue(TestUtilities.IsNearlyEqual(r0, result.Y, settings));

            Assert.IsTrue(count > 0);
        }


        [TestMethod]
        public void OdeAiry () {

            // This is the airy differential equation
            Func<double, double, double> f = (double x, double y) => x * y;
            
            // Solutions should be of the form f(x) = a Ai(x) + b Bi(x).
            // Given initial value of f and f', this equation plus the Wronskian can be solved to give
            //   a = \pi ( f Bi' - f' Bi )    b = \pi ( f' Ai - f Ai' )

            // Start with some initial conditions
            double x0 = 0.0;
            double y0 = 0.0;
            double yp0 = 1.0;

            // Find the a and b coefficients consistent with those values
            SolutionPair s0 = AdvancedMath.Airy(x0);
            double a = Math.PI * (y0 * s0.SecondSolutionDerivative - yp0 * s0.SecondSolutionValue);
            double b = Math.PI * (yp0 * s0.FirstSolutionValue - y0 * s0.FirstSolutionDerivative);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(y0, a * s0.FirstSolutionValue + b * s0.SecondSolutionValue));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(yp0, a * s0.FirstSolutionDerivative + b * s0.SecondSolutionDerivative));

            // Integrate to a new point (pick a negative one so we test left integration)
            double x1 = -5.0;
            OdeResult result = FunctionMath.IntegrateConservativeOde(f, x0, y0, yp0, x1);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(result.X, x1));
            Console.WriteLine(result.EvaluationCount);

            // The solution should still hold
            SolutionPair s1 = AdvancedMath.Airy(x1);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Y, a * s1.FirstSolutionValue + b * s1.SecondSolutionValue, result.Settings));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(result.YPrime, a * s1.FirstSolutionDerivative + b * s1.SecondSolutionDerivative, result.Settings));

        }

        [TestMethod]
        public void OdeLRC () {

            // Demanding zero net voltage across L, R, and C elements in series gives
            //   Q / C + \dot{Q} R + \ddot{Q} L = 0
            // This is a second order linear ODE with constant coefficients, i.e.
            // universal damped harmonic oscillator.

            // Characteristic equation r^2 L + R r + r / C = 0 with solutions
            //   r = \frac{-R \pm \sqrt{R^2 - 4 L / C}{2L}. Define
            // t_0 = \sqrt{RC} = 1 / w_0 and t_1 = 2 L / R = 1 / w_1.

            // If t_0 < t_1 the discriminant is negative and
            //   r = - w_0 \pm i w
            // where w = \sqrt{ w_0^2 - w_1^2 } so
            //   Q = e^{-w_0 t} \left[ A cos(w t) + B sin(w t) \right]
            // We get damped oscillatory behavior.

            // If t_0 > t_1 the discriminant is positive and
            //   r = \sqrt{ w_1^2 - w_0^2 } \pm w_1
            // so
            //   Q = A e^{-w_a t} + B^{-w_b t}
            // We get purely damped behavior. 

            double q0 = 1.0;
            double qp0 = 0.0;

            foreach(double L in TestUtilities.GenerateRealValues(0.1, 10.0, 4, 1)) {
                foreach (double R in TestUtilities.GenerateRealValues(0.1, 1.0, 4, 2)) {
                    foreach (double C in TestUtilities.GenerateRealValues(0.1, 1.0, 4, 3)) {

                        Func<double, IReadOnlyList<double>, IReadOnlyList<double>> rhs = (double x, IReadOnlyList<double> y) => {
                            double q = y[0]; double qp = y[1];
                            return (new double[] {
                                qp,
                                (0.0 - q / C - R * qp) / L
                            });
                        };

                        double t0 = Math.Sqrt(L * C);
                        double w0 = 1.0 / t0;

                        double t1 = 2.0 * L / R;
                        double w1 = 1.0 / t1;

                        double t = 4.0;
                        double qt;
                        if (t0 < t1) {

                            double w = Math.Sqrt(w0 * w0 - w1 * w1);

                            double A = q0;
                            double B = (qp0 + w1 * q0) / w;

                            qt = (A * Math.Cos(w * t) + B * Math.Sin(w * t)) * Math.Exp(-w1 * t);

                        } else {

                            double w = Math.Sqrt(w1 * w1 - w0 * w0);
                            double wa = w1 + w;
                            double wb = w0 * w0 / wa;

                            double A = (wb * q0 + qp0) / (2.0 * w);
                            double B = (wa * q0 + qp0) / (2.0 * w);

                            qt = -A * Math.Exp(-wa * t) + B * Math.Exp(-wb * t);
                        }
                        Console.WriteLine(qt);

                        MultiOdeResult result = MultiFunctionMath.IntegrateOde(rhs, 0.0, new double[] { q0, qp0 }, t);
                        Console.WriteLine(result.Y[0]);

                        Console.WriteLine(result.EvaluationCount);

                        Assert.IsTrue(TestUtilities.IsNearlyEqual(qt, result.Y[0], result.Settings));

                    }
                }
            }

        }

        [TestMethod]
        public void OdeCatenary () {

            // The equation of a catenary is
            //  \frac{d^2 y}{dx^2} = \sqrt{1.0 + \left(\frac{dy}{dx}\right)^2}
            // with y(0) = 1 and y'(0) = 0 and solution y = \cosh(x)

            Func<double, IReadOnlyList<double>, IReadOnlyList<double>> rhs = (double t, IReadOnlyList<double> y) =>
                new double[] { y[1], MoreMath.Hypot(1.0, y[1]) };

            MultiOdeResult result = MultiFunctionMath.IntegrateOde(rhs, 0.0, new double[] { 1.0, 0.0 }, 2.0);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Y[0], Math.Cosh(2.0)));

        }

        [TestMethod]
        public void OdeLorenz () {

            double sigma = 10.0;
            double beta = 8.0 / 3.0;
            double rho = 28.0;

            Func<double, IReadOnlyList<double>, IReadOnlyList<double>> rhs = (double t, IReadOnlyList<double> u) =>
                new ColumnVector(
                    sigma * (u[1] - u[0]),
                    u[0] * (rho - u[2]) - u[1],
                    u[0] * u[1] - beta * u[2]
                );

            MultiOdeSettings settings = new MultiOdeSettings() {
                RelativePrecision = 1.0E-8,
                EvaluationBudget = 10000
            };

            ColumnVector u0 = new ColumnVector(1.0, 1.0, 1.0);

            MultiOdeResult result = MultiFunctionMath.IntegrateOde(rhs, 0.0, u0, 10.0, settings);

            Console.WriteLine(result.EvaluationCount);

            // Is there anything we can assert? There is no analytic solution or conserved quantity.

        }

        [TestMethod]
        public void OdeLaneEmden () {

            // The Lane-Emden equations describe a simplified model of stellar structure.
            // See http://mathworld.wolfram.com/Lane-EmdenDifferentialEquation.html

            // Analytic solutions are known for the n=0, 1, and 5 cases.

            int n = 0;

            Func<double, IReadOnlyList<double>, IReadOnlyList<double>> rhs = (double x, IReadOnlyList<double> t) =>
                new ColumnVector(
                    t[1], -MoreMath.Pow(t[0], n) - 2.0 * (x == 0.0 ? 0.0 : t[1] / x)
                );

            ColumnVector t0 = new ColumnVector(1.0, 0.0);

            double x1 = 2.0;

            n = 0;
            MultiOdeResult result = MultiFunctionMath.IntegrateOde(rhs, 0.0, t0, x1);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Y[0], 1.0 - MoreMath.Sqr(x1) / 6.0));
            Console.WriteLine(result.EvaluationCount);

            n = 1;
            result = MultiFunctionMath.IntegrateOde(rhs, 0.0, t0, x1);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Y[0], MoreMath.Sin(x1) / x1));
            Console.WriteLine(result.EvaluationCount);

            n = 5;
            result = MultiFunctionMath.IntegrateOde(rhs, 0.0, t0, x1);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Y[0], 1.0 / Math.Sqrt(1.0 + MoreMath.Sqr(x1) / 3.0)));
            Console.WriteLine(result.EvaluationCount);

            // For all of these cases, the initial step fails until it gets very small, then it increases again.
            // Look into why this is.

        }

        [TestMethod]
        public void OdeLotkaVolterra () {

            // Lotka/Volterra equations are a non-linear predator-prey model.
            //   \dot{x} = A x + B x y
            //   \dot{y} = -C y + D x y
            // See http://mathworld.wolfram.com/Lotka-VolterraEquations.html

            // It can be shown solutions are always periodic (but not sinusoidal, and
            // often with phase shift between x and y).

            // Equilibria are x, y = C/D, A /B (and 0, 0).
            // If started positive, x and y never go negative.
            // A conserved quantity is: - D x + C log x - B y + A log y
            // Period of equations linearized around equilibrium is 2 \pi / \sqrt{AC}.

            // Can also be used to model chemical reaction rates.

            double A = 1.5; // Prey growth rate
            double B = 1.0;
            double C = 3.0; // Predator death rate
            double D = 1.0;

            // Try A = 20, B = 1, C = 30, D = 1, X = 8, Y = 12, T = 1
            // Try A = 1.5, B = 1 C = 3, D = 1, X = 10, Y = 5, T = 10
 
            Func<double, IReadOnlyList<double>, IReadOnlyList<double>> rhs = (double t, IReadOnlyList<double> p) =>
                new double[] {
                    A * p[0] - B * p[0] * p[1],
                    -C * p[1] + D * p[0] * p[1]
                };

            Func<IList<double>, double> conservedQuantity = (IList<double> p) =>
                -D * p[0] + C * Math.Log(p[0]) - B * p[1] + A * Math.Log(p[1]);

            ColumnVector p0 = new ColumnVector(10.0, 5.0);

            double L0 = conservedQuantity(p0);

            // Set up a handler that verifies conservation and positivity
            MultiOdeSettings settings = new MultiOdeSettings() {
                Listener = (MultiOdeResult rr) => {
                    double L = conservedQuantity(rr.Y);
                    Assert.IsTrue(rr.Y[0] > 0.0);
                    Assert.IsTrue(rr.Y[1] > 0.0);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(L, L0, rr.Settings));
                }
            };

            // Estimate period
            double T = 2.0 * Math.PI / Math.Sqrt(A * C);

            // Integrate over a few estimated periods
            MultiOdeResult result = MultiFunctionMath.IntegrateOde(rhs, 0.0, p0, 3.0 * T, settings);

            Console.WriteLine(result.EvaluationCount);

        }

        [TestMethod]
        public void OdeDawson () {

            // The Dawson function fulfills a simple ODE.
            //   \frac{dF}{dx} + 2 x F = 1 \qquad F(0) = 0
            // See e.g. https://en.wikipedia.org/wiki/Dawson_function
            // Verify that we get correct values via ODE integration.

            Func<double, double, double> rhs = (double x, double F) => 1.0 - 2.0 * x * F;

            foreach (double x1 in TestUtilities.GenerateRealValues(0.1, 10.0, 8)) {

                EvaluationSettings s = new EvaluationSettings() {
                    RelativePrecision = 1.0E-13,
                    AbsolutePrecision = 0.0
                };
                OdeResult r = FunctionMath.IntegrateOde(rhs, 0.0, 0.0, x1);
                Debug.WriteLine("{0}: {1} {2}: {3}", x1, r.Y, AdvancedMath.Dawson(x1), r.EvaluationCount);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(r.Y, AdvancedMath.Dawson(x1), s));

            }
        }
          
        [TestMethod]
        public void OdePendulum () {

            // Without the small-angle approximation, the period of a pendulum is not
            // independent of angle. Instead, it is given by
            //   P = 4 K(\sin(\phi_0) / 2)
            // where \phi_0 is the release angle and K is the complete elliptic function
            // of the first kind.
            // See https://en.wikipedia.org/wiki/Pendulum_(mathematics)
            // We compute for an initial angle of 45 degrees, far beyond a small angle.

            Func<double, double, double> rhs = (double t, double u) => -MoreMath.Sin(u);

            double u0 = Math.PI / 4.0;
            double p = 4.0 * AdvancedMath.EllipticK(MoreMath.Sin(u0 / 2.0));

            OdeResult r = FunctionMath.IntegrateConservativeOde(rhs, 0.0, u0, 0.0, p);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(r.Y, u0));

            Console.WriteLine(r.EvaluationCount);

        }

        [TestMethod]
        public void OdeErf  () {

            Func<double, double, double> rhs = (double t, double u) =>
                2.0 / Math.Sqrt(Math.PI) * Math.Exp(-t * t);

            OdeResult r = FunctionMath.IntegrateOde(rhs, 0.0, 0.0, 5.0);

             Console.WriteLine(r.Y);

        }

        [TestMethod]
        public void OdeSolarSystem () {

            Planet[] planets = new Planet[] {

                new Planet() {
                    Name = "Sun",
                    Mass = 1.3271244004193938E11 * 2.22972472E-15,
                    Position = new ColumnVector(3.700509269632818E-03, 2.827000367199164E-03, -1.623212133858169E-04),
                    Velocity = new ColumnVector(-1.181051079944745E-06, 7.011580463060376E-06, 1.791336618633265E-08),
                    Year = 0.0
                },
                new Planet() {
                    Name = "Earth",
                    Mass = 398600.440 * 2.22972472E-15,
                    Position = new ColumnVector(5.170309635282939E-01, -8.738395510520275E-01, -1.323433043109283E-04),
                    Velocity = new ColumnVector(1.456221287512929E-02, 8.629625079574064E-03, 2.922661879104068E-07),
                    Year = 365.25636
                },
                new Planet() {
                    Name = "Jupiter",
                    Mass = 126686511 * 2.22972472E-15,
                    Position = new ColumnVector(-5.438893557878444E+00, 1.497713688978628E-01, 1.210124167423688E-01),
                    Velocity = new ColumnVector(-2.955588275385831E-04, -7.186425047188191E-03, 3.648630400553893E-05),
                    Year = 4332.59
                },
                new Planet() {
                    Name = "Saturn",
                    Mass = 37931207.8 * 2.22972472E-15,
                    Position = new ColumnVector(-2.695377264613252E+00, -9.657410990313211E+00, 2.751886819002198E-01),
                    Velocity = new ColumnVector(5.067197776417300E-03, -1.516587142748002E-03, -1.754902991309407E-04),
                    Year = 10759.22
                }
            };

            Func<double, IReadOnlyList<double>, IReadOnlyList<double>> rhs = (double t, IReadOnlyList<double> p) => {

                ColumnVector a = new ColumnVector(3 * planets.Length);

                for (int i = 0; i < planets.Length; i++) {

                    for (int j = 0; j < planets.Length; j++) {

                        if (i == j) continue;

                        double x = p[3 * i] - p[3 * j];
                        double y = p[3 * i + 1] - p[3 * j + 1];
                        double z = p[3 * i + 2] - p[3 * j + 2];
                        double d3 = Math.Pow(x * x + y * y + z * z, 3.0 / 2.0);

                        a[3 * i] -= planets[j].Mass * x / d3;
                        a[3 * i + 1] -= planets[j].Mass * y / d3;
                        a[3 * i + 2] -= planets[j].Mass * z / d3;

                    }

                }

                return (a);
            };

            ColumnVector p0 = new ColumnVector(3 * planets.Length);
            for(int i = 0; i < planets.Length; i++) {
                p0[3 * i] = planets[i].Position[0];
                p0[3 * i + 1] = planets[i].Position[1];
                p0[3 * i + 2] = planets[i].Position[2];
            }

            ColumnVector q0 = new ColumnVector(3 * planets.Length);
            for (int i = 0; i < planets.Length; i++) {
                q0[3 * i] = planets[i].Velocity[0];
                q0[3 * i + 1] = planets[i].Velocity[1];
                q0[3 * i + 2] = planets[i].Velocity[2];
            }

            MultiOdeSettings settings = new MultiOdeSettings() {
                RelativePrecision = 1.0E-8,
                AbsolutePrecision = 1.0E-16
            };

            MultiOdeResult result = MultiFunctionMath.IntegrateConservativeOde(rhs, 0.0, p0, q0, 4332.59, settings);

        }

    }

    // Positions and velocities for all solar system bodies at any time
    // are available at http://ssd.jpl.nasa.gov/horizons.cgi
    // Position units are AU, velocity units are AU/day
    // Mass is given as the product GM in km^3/s^2, converted to AU^3/day^2.


    internal class Planet {

        public string Name { get; set; }

        public double Mass { get; set; }

        public IList<double> Position { get; set; }

        public IList<double> Velocity { get; set; }

        public double Year { get; set; }

    }
}
