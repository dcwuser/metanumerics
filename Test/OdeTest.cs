using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;

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

            double cz = 1.0;
            double cx = 0.0;
            double cy = 0.0;

            Random rng = new Random(1001110000);
            Distribution xDistribution = new UniformDistribution(Interval.FromEndpoints(-4.0, 8.0));
            Distribution yDistribution = new UniformDistribution(Interval.FromEndpoints(-8.0, 4.0));
            Distribution eDistribution = new NormalDistribution();

            Sample r2Sample = new Sample();

            for (int i = 0; i < 500; i++) {

                MultivariateSample xyzSample = new MultivariateSample(3);
                for (int k = 0; k < 12; k++) {
                    double x = xDistribution.GetRandomValue(rng);
                    double y = yDistribution.GetRandomValue(rng);
                    double z = cx * x + cy * y + cz + eDistribution.GetRandomValue(rng);
                    xyzSample.Add(x, y, z);
                }
                FitResult fit = xyzSample.LinearRegression(2);
                double fcx = fit.Parameters[0];
                double fcy = fit.Parameters[1];
                double fcz = fit.Parameters[2];

                double ss2 = 0.0;
                double ss1 = 0.0;
                foreach (double[] xyz in xyzSample) {
                    ss2 += MoreMath.Sqr(xyz[2] - (fcx * xyz[0] + fcy * xyz[1] + fcz));
                    ss1 += MoreMath.Sqr(xyz[2] - xyzSample.Column(2).Mean);
                }
                double r2 = 1.0 - ss2 / ss1;
                r2Sample.Add(r2);
            }

            Console.WriteLine("{0} {1} {2} {3} {4}", r2Sample.Count, r2Sample.PopulationMean, r2Sample.StandardDeviation, r2Sample.Minimum, r2Sample.Maximum);

            Distribution r2Distribution = new BetaDistribution((3 - 1) / 2.0, (12 - 3) / 2.0);
            //Distribution r2Distribution = new BetaDistribution((10 - 2) / 2.0, (2 - 1) / 2.0);
            Console.WriteLine("{0} {1}", r2Distribution.Mean, r2Distribution.StandardDeviation);

            TestResult ks = r2Sample.KolmogorovSmirnovTest(r2Distribution);
            Console.WriteLine(ks.RightProbability);
            Console.WriteLine(ks.Probability);

        }


        [TestMethod]
        public void TestBivariateRegression () {

            double a0 = 1.0;
            double b0 = 0.0;

            Random rng = new Random(1001110000);
            Distribution xDistribution = new UniformDistribution(Interval.FromEndpoints(-2.0, 4.0));
            Distribution eDistribution = new NormalDistribution();

            Sample r2Sample = new Sample();

            for (int i = 0; i < 500; i++) {

                BivariateSample xySample = new BivariateSample();
                for (int k = 0; k < 10; k++) {
                    double x = xDistribution.GetRandomValue(rng);
                    double y = a0 + b0 * x + eDistribution.GetRandomValue(rng);
                    xySample.Add(x, y);
                }
                FitResult fit = xySample.LinearRegression();
                double a = fit.Parameters[0];
                double b = fit.Parameters[1];

                double ss2 = 0.0;
                double ss1 = 0.0;
                foreach (XY xy in xySample) {
                    ss2 += MoreMath.Sqr(xy.Y - (a + b * xy.X));
                    ss1 += MoreMath.Sqr(xy.Y - xySample.Y.Mean);
                }
                double r2 = 1.0 - ss2 / ss1;
                r2Sample.Add(r2);
            }

            Console.WriteLine("{0} {1} {2} {3} {4}", r2Sample.Count, r2Sample.PopulationMean, r2Sample.StandardDeviation, r2Sample.Minimum, r2Sample.Maximum);

            Distribution r2Distribution = new BetaDistribution((2 - 1) / 2.0, (10 - 2) / 2.0);
            //Distribution r2Distribution = new BetaDistribution((10 - 2) / 2.0, (2 - 1) / 2.0);
            Console.WriteLine("{0} {1}", r2Distribution.Mean, r2Distribution.StandardDeviation);

            TestResult ks = r2Sample.KolmogorovSmirnovTest(r2Distribution);
            Console.WriteLine(ks.RightProbability);
            Console.WriteLine(ks.Probability);

        }

        [TestMethod]
        public void OdeExponential () {

            // Exponential
            // y = y_0 e^{x - x_0}
            Func<double, double, double> f = (double x, double y) => y;

            double y1 = FunctionMath.SolveOde(f, 0.0, 1.0, 2.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(y1, MoreMath.Sqr(Math.E)));
        }

        [TestMethod]
        public void OdeNonlinear () {

            // y = \frac{y_0}{1 - y_0 (x - x_0)}
            Func<double, double, double> f = (double x, double y) => MoreMath.Sqr(y);

            EvaluationSettings settings = new EvaluationSettings() { RelativePrecision = 1.0E-8, EvaluationBudget = 1000 };
            double y1 = FunctionMath.SolveOde(f, 0.0, 1.0, 0.99, settings);
            Console.WriteLine(y1);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(y1, 1.0 / (1.0 - 1.0 * (0.99 - 0.0)), settings));

        }

        [TestMethod]
        public void OdeLogistic () {

            // y = \frac{y_0}{y_0 + (1 - y_0) e^{-(x - x_0)}
            Func<double, double, double> rhs = (double x, double y) => y * (1.0 - y);
            Func<double, double, double> solution = (double y0, double x) => y0 / (y0 + (1.0 - y0) * Math.Exp(-x));

            foreach (double y0 in new double[] { -0.1, 0.0, 0.4, 1.0, 1.6 }) {
                Console.WriteLine(y0);
                Interval r = Interval.FromEndpoints(0.0, 2.0);
                double y1 = FunctionMath.SolveOde(rhs, 0.0, y0, 2.0);
                Console.WriteLine(y1);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(y1, solution(y0, 2.0)));
            }

        }

        [TestMethod]
        public void OdeExample () {

            Func<double, double, double> rhs = (double t, double u) => (1.0 - 2.0 * t) * u;
            Func<double, double, double> solution = (double u0, double t) => u0 * Math.Exp(t - t * t);

            foreach (double t in new double[] { 0.5, 0.75, 1.50, 2.25, 3.25 }) {
                double y1 = FunctionMath.SolveOde(rhs, 0.0, 1.0, t);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(y1, solution(1.0, t)));
            }

        }

        [TestMethod]
        public void OdeSine () {
            Func<double, double, double> f = (double x, double y) => -y;
            double y1 = FunctionMath.SolveConservativeOde(f, 0.0, 0.0, 1.0, 5.0);
            Console.WriteLine(y1);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(y1, MoreMath.Sin(5.0)));
        }

        [TestMethod]
        public void OdeSine2 () {
            Func<double, IList<double>, IList<double>> f2 = (double x, IList<double> y) => new double[] { -y[0] };
            double[] y0 = new double[] { 0.0 };
            double[] yp0 = new double[] { 1.0 };
            ColumnVector z = MultiFunctionMath.SolveConservativeOde(f2, 0.0, y0, yp0, 5.0);
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
            //   I_1 \omega_1^2 + I_2 \omega_2^2 + I_3 \omega_3^2 = 2E
            // and angular momentum
            //   I_1^2 \omega_1^2 + I_2^2 \omega_2^2 + I_3^2 \omega_3^2 = M^2
            // are conserved quantities.

            ColumnVector I = new ColumnVector(1.0, 2.0, 3.0);

            Func<double, IList<double>, IList<double>> rhs = (double t, IList<double> w) => {
                return (new ColumnVector((I[1] - I[2]) * w[1] * w[2] / I[0], (I[2] - I[0]) * w[2] * w[0] / I[1], (I[0] - I[1]) * w[0] * w[1] / I[2]));
            };

            ColumnVector w0 = new ColumnVector(1.0, 1.0, 1.0);

            double eps = 1.0E-12;
            EvaluationSettings settings = new EvaluationSettings() {
                RelativePrecision = eps,
                AbsolutePrecision = eps * eps,
                EvaluationBudget = 10000
            };
            settings.UpdateHandler = (EvaluationResult r) => {
                OdeResult<IList<double>> q = (OdeResult<IList<double>>) r;
                Console.WriteLine("{0} {1}", q.EvaluationCount, q.X);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(EulerKineticEnergy(I, q.Y), EulerKineticEnergy(I, w0), settings));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(EulerAngularMomentum(I, q.Y), EulerAngularMomentum(I, w0), settings));
            };

            MultiFunctionMath.SolveOde(rhs, 0.0, w0, 10.0, settings);

        }

        [TestMethod]
        public void OdeOrbit2 () {

            Func<double, IList<double>, IList<double>> rhs = (double t, IList<double> r) => {
                double d = MoreMath.Hypot(r[0], r[2]);
                double d3 = MoreMath.Pow(d, 3);
                return (new ColumnVector(r[1], -r[0] / d3, r[3], -r[2] / d3));
            };

            //double e = 0.5;
            foreach (double e in TestUtilities.GenerateUniformRealValues(0.0, 1.0, 8)) {
                Console.WriteLine(e);

                ColumnVector r0 = new ColumnVector(1.0 - e, 0.0, 0.0, Math.Sqrt((1.0 + e) / (1.0 - e)));

                ColumnVector r1 = new ColumnVector(MultiFunctionMath.SolveOde(rhs, 0.0, r0, 2.0 * Math.PI).Y);

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
        public void OdeOrbit () {

            // This is a simple Keplerian orbit.
            // Hull (1972) constructed initial conditions that guarantee a given orbital eccentricity
            // and a period of 2 \pi with unit masses and unit gravitational constant.

            Func<double, IList<double>, IList<double>> rhs = (double t, IList<double> r) => {
                double d = MoreMath.Hypot(r[0], r[1]);
                double d3 = MoreMath.Pow(d, 3);
                return (new double[] { -r[0] / d3, -r[1] / d3 });
            };

            EvaluationSettings settings = new EvaluationSettings() { RelativePrecision = 1.0E-12, AbsolutePrecision = 1.0E-24, EvaluationBudget = 10000 };
            settings.UpdateHandler = (EvaluationResult a) => {
                OdeResult<IList<double>> b = (OdeResult<IList<double>>) a;
                Console.WriteLine("{0} {1}", b.EvaluationCount, b.X);
            };

            //double e = 0.5;
            foreach (double e in TestUtilities.GenerateUniformRealValues(0.0, 1.0, 8)) {
                Console.WriteLine("e = {0}", e);

                ColumnVector r0 = new ColumnVector(1.0 - e, 0.0);
                ColumnVector rp0 = new ColumnVector(0.0, Math.Sqrt((1.0 + e) / (1.0 - e)));

                settings.UpdateHandler = (EvaluationResult a) => {
                    OdeResult<IList<double>> b = (OdeResult<IList<double>>) a;
                    Console.WriteLine("  {0} {1}", b.EvaluationCount, b.X);
                    Console.WriteLine("  {0} ?= {1}", OrbitEnergy(b.Y, b.YPrime), OrbitEnergy(r0, rp0)); 
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(OrbitAngularMomentum(b.Y, b.YPrime), OrbitAngularMomentum(r0, rp0), settings));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(OrbitEnergy(b.Y, b.YPrime), OrbitEnergy(r0, rp0), new EvaluationSettings() { RelativePrecision = 2.0E-12 }));
                };

                ColumnVector r1 = MultiFunctionMath.SolveConservativeOde(rhs, 0.0, r0, rp0, 2.0 * Math.PI, settings);

                // The 

                Assert.IsTrue(TestUtilities.IsNearlyEqual(r0, r1, new EvaluationSettings() { RelativePrecision = 1.0E-9 }));

            }
        }

        [TestMethod]
        public void OdeAiry () {
            Func<double, double, double> f = (double x, double y) => x * y;
            FunctionMath.SolveConservativeOde(f, 0.0, 1.0 / (Math.Pow(3.0, 2.0 / 3.0) * AdvancedMath.Gamma(2.0 / 3.0)), - 1.0 / (Math.Pow(3.0, 1.0 / 3.0) * AdvancedMath.Gamma(1.0 / 3.0)), -5.0);
        }

        [TestMethod]
        public void OdeLRC () {

            double L = 1.0;
            double R = 1.0;
            double C = 1.0;

            Func<double, IList<double>, IList<double>> rhs = (double x, IList<double> y) => {
                double q = y[0]; double qp = y[1];
                return (new double[] {
                    qp,
                    (0.0 - q / C - R * qp) / L
                });
            };

            double q0 = 1.0;
            double qp0 = 0.0;

            double t0 = 2.0 * L / R;
            double w0 = Math.Sqrt(1.0 / (L * C) - 1.0 / (t0 * t0));

            double t = 1.0;
            double sc = (q0 * Math.Cos(w0 * t) + (qp0 + q0 / t0) / w0 * Math.Sin(w0 * t)) * Math.Exp(-t / t0);
            Console.WriteLine(sc);

            IList<double> s = MultiFunctionMath.SolveOde(rhs, 0.0, new double[] { q0, qp0 }, t).Y;
            Console.WriteLine(s[0]);

        }

        [TestMethod]
        public void OdeLorenz () {

            double sigma = 10.0;
            double beta = 8.0 / 3.0;
            double rho = 28.0;

            Func<double, IList<double>, IList<double>> rhs = (double t, IList<double> u) => new ColumnVector(
                sigma * (u[1] - u[0]),
                u[0] * (rho - u[2]) - u[1],
                u[0] * u[1] - beta * u[2]
            );

            ColumnVector u0 = new ColumnVector(1.0, 1.0, 1.0);

            MultiFunctionMath.SolveOde(rhs, 0.0, u0, 10.0, new EvaluationSettings() { RelativePrecision = 1.0E-8, EvaluationBudget = 100000 });

        }
    

    }
}
