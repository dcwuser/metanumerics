using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;


namespace Test {
    
    [TestClass]
    public class RootsTest {

        [TestMethod]
        public void RootOfEi () {
            double x = FunctionMath.FindZero(AdvancedMath.IntegralEi, Interval.FromEndpoints(0.1, 1.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(x, 0.37250741078136663446));
            Assert.IsTrue(Math.Abs(AdvancedMath.IntegralEi(x)) < TestUtilities.TargetPrecision);
        }

        [TestMethod]
        public void RootOfPsi () {
            double x = FunctionMath.FindZero(AdvancedMath.Psi, 1.5);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(x, 1.46163214496836234126));
            Assert.IsTrue(AdvancedMath.Psi(x) < TestUtilities.TargetPrecision);
        }

        [TestMethod]
        public void RootOfJ0 () {
            Func<double, double> f = delegate(double x) {
                return (AdvancedMath.BesselJ(0, x));
            };
            double y = FunctionMath.FindZero(f, Interval.FromEndpoints(2.0, 4.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(y, 2.40482555769577276862)); 
        }

        [TestMethod]
        public void MultiRootMathworksExample () {

            Func<double[], double[]> f = delegate(double[] u) {
                double x = u[0]; double y = u[1];
                return (new double[] { 2.0 * x - y - Math.Exp(-x), 2.0 * y - x - Math.Exp(-y) });
            };

            double[] x0 = GetRandomStartingVector(new Random(1), 2);

            double[] xz = FunctionMath.FindZero(f, x0);

            double[] fz = f(xz);
            Assert.IsTrue(Math.Abs(fz[0]) < TestUtilities.TargetPrecision);
            Assert.IsTrue(Math.Abs(fz[1]) < TestUtilities.TargetPrecision);

        }

        [TestMethod]
        public void MultiRootExtendedRosenbock () {

            Random rng = new Random(1);

            // dimension must be even for this one
            for (int d = 2; d <= 10; d += 2) {
 
                Func<double[], double[]> f = delegate(double[] u) {
                    double[] v = new double[d];
                    for (int k = 0; k < d / 2; k++) {
                        v[2 * k] = 1.0 - u[2 * k + 1];
                        v[2 * k + 1] = 10.0 * (u[2 * k] - MoreMath.Pow(u[2 * k + 1], 2));
                    }
                    return (v);
                };

                double[] x = FunctionMath.FindZero(f, GetRandomStartingVector(rng, d));

                // the solution is (1, 1, ..., 1)
                for (int i = 0; i < d; i++) Assert.IsTrue(TestUtilities.IsNearlyEqual(x[i], 1.0));

                // and of course the function should evaluate to zero
                double[] y = f(x);
                for (int i = 0; i < d; i++) Assert.IsTrue(Math.Abs(y[i]) < TestUtilities.TargetPrecision);
            }

        }

        [TestMethod]
        public void MultiRootTrigExp () {

            Random rng = new Random(1);

            foreach (int d in new int[] { 3, 5, 8 }) {

                Func<double[], double[]> f = delegate(double[] x) {
                    double[] y = new double[d];
                    y[0] = 3.0 * MoreMath.Pow(x[0], 3) + 2.0 * x[1] - 5.0 + Math.Sin(x[0] - x[1]) * Math.Sin(x[0] + x[1]);
                    for (int k = 1; k < d - 1; k++) {
                        y[k] = -x[k - 1] * Math.Exp(x[k - 1] - x[k]) + x[k] * (4.0 + 3.0 * MoreMath.Pow(x[k], 2)) + 2.0 * x[k + 1] +
                            Math.Sin(x[k] - x[k + 1]) * Math.Sin(x[k] + x[k + 1]) - 8.0;
                    }
                    y[d - 1] = -x[d - 2] * Math.Exp(x[d - 2] - x[d - 1]) + 4.0 * x[d - 1] - 3.0;
                    return (y);
                };

                double[] x0 = GetRandomStartingVector(rng, d);

                double[] xz = FunctionMath.FindZero(f, x0);

                // solution is (1, 1, ..., 1)
                for (int i = 0; i < d; i++) Assert.IsTrue(TestUtilities.IsNearlyEqual(xz[i], 1.0));
                
                // and function there is of course zero
                double[] fz = f(xz);
                for (int i = 0; i < d; i++) Assert.IsTrue(Math.Abs(fz[i]) < TestUtilities.TargetPrecision);

            }

        }

        // create a random starting vector with components in the box -10 < x_k < +10

        private double[] GetRandomStartingVector (Random rng, int d) {
            double[] x = new double[d];
            for (int k = 0; k < d; k++) {
                x[k] = 10.0 * rng.NextDouble() - 5.0;
            }
            return (x);
        }

    }
}
