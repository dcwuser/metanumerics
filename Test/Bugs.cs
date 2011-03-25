using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public class BugTests {

        [TestMethod]
        public void Bug2811 () {

            ChiSquaredDistribution d = new ChiSquaredDistribution(1798);
            double x = d.InverseLeftProbability(0.975);
            Console.WriteLine(x);
        }

        private SquareMatrix CyclicMatrix (int n) {

            SquareMatrix C = new SquareMatrix(n);
            for (int c = 1; c < n; c++) {
                C[c - 1, c] = 1.0;
            }
            C[n - 1, 0] = 1.0;

            return (C);
        }

        [TestMethod]
        public void Bug5504 () {
            // this eigenvalue non-convergence is solved by an ad hoc shift
            // these "cyclic matrices" are good examples of situations that are difficult for the QR algorithm 
            for (int d = 2; d <= 8; d++) {
                Console.WriteLine(d);
                SquareMatrix C = CyclicMatrix(d);
                Complex[] lambdas = C.Eigenvalues();
                foreach (Complex lambda in lambdas) {
                    Console.WriteLine("{0} ({1} {2})", lambda, ComplexMath.Abs(lambda), ComplexMath.Arg(lambda));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Abs(lambda), 1.0));
                }
            }

        }

        [TestMethod]
        public void DifficultEigenvalue () {

            SquareMatrix A = new SquareMatrix(4);
            A[0, 0] = 0.0; A[0, 1] = 2.0; A[0, 2] = 0.0; A[0, 3] = -1.0;
            A[1, 0] = 1.0; A[1, 1] = 0.0; A[1, 2] = 0.0; A[1, 3] = 0.0;
            A[2, 0] = 0.0; A[2, 1] = 1.0; A[2, 2] = 0.0; A[2, 3] = 0.0;
            A[3, 0] = 0.0; A[3, 1] = 0.0; A[3, 2] = 1.0; A[3, 3] = 0.0;

            Complex[] zs = A.Eigenvalues();
            foreach (Complex z in zs) {
                Console.WriteLine("{0} ({1} {2})", z, ComplexMath.Abs(z), ComplexMath.Arg(z));
            }


        }

        [TestMethod]
        public void Bug5505 () {
            // finding the root of x^3 would fail with a NonconvergenceException because our termination criterion tested only
            // for small relative changes in x; adding a test for small absolute changes too solved the problem
            double x0 = FunctionMath.FindZero(delegate(double x) { return (x * x * x); }, 1.0);
            Assert.IsTrue(Math.Abs(x0) < TestUtilities.TargetPrecision);
        }

        [TestMethod]
        public void Bug5705 () {
            // large arguments to Erf and Erfc would cause a NonconvergeException because of the way the continued fraction for GammaQ
            // was evaluated; we added an explicit test before entering the evaluation loop to avoid this
            Assert.IsTrue(AdvancedMath.Erf(Double.NegativeInfinity) == -1.0); Assert.IsTrue(AdvancedMath.Erfc(Double.NegativeInfinity) == 2.0);
            Assert.IsTrue(AdvancedMath.Erf(Double.MinValue) == -1.0); Assert.IsTrue(AdvancedMath.Erfc(Double.MinValue) == 2.0);
            Assert.IsTrue(AdvancedMath.Erf(Double.MaxValue) == 1.0); Assert.IsTrue(AdvancedMath.Erfc(Double.MaxValue) == 0.0);
            Assert.IsTrue(AdvancedMath.Erf(Double.PositiveInfinity) == 1.0); Assert.IsTrue(AdvancedMath.Erfc(Double.PositiveInfinity) == 0.0);
            // we should add tests of the behavior of all our advanced functions with large/infinite arguments
        }

        // 1 sqrt ~ 4 flops
        // 1 log ~ 12 flops

        // Inverse CDF of F-distribution fails for d2 <= 2

        [TestMethod]
        public void Bug5886 () {

            double nu1 = 1.0;
            double nu2 = 0.1;

            FisherDistribution F = new FisherDistribution(nu1, nu2);

            /*
            Console.WriteLine(F.Mean);
            Console.WriteLine(F.StandardDeviation);
            Console.WriteLine(F.Skewness);

            Console.WriteLine(F.LeftProbability(0.1));
            Console.WriteLine(F.LeftProbability(1.0));
            Console.WriteLine(F.LeftProbability(1.0E9));

            BetaDistribution B = new BetaDistribution(F.NumeratorDegreesOfFreedom / 2.0, F.DenominatorDegreesOfFreedom / 2.0);
            double x = 1.3;
            double y = nu1 * x / (nu2 + nu1 * x);
            double u = nu1 * nu2 / MoreMath.Pow(nu2 + nu1 * x, 2);
            Console.WriteLine("{0} v. {1}", F.ProbabilityDensity(x), u * B.ProbabilityDensity(y));
            */

            double x1 = F.InverseLeftProbability(0.6);
            Console.WriteLine(x1);
            double P = F.LeftProbability(x1);
            Console.WriteLine(P);

        }

        public static double InverseFisherDistribution(FisherDistribution F, double P) {
            if (F == null) throw new ArgumentNullException("F");
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");

            double x0;
            if (F.DenominatorDegreesOfFreedom < 3.0) {
                x0 = 3.0;
            } else {
                x0 = F.Mean;
            }

            double x1 = FunctionMath.FindZero(delegate(double x) { return (F.LeftProbability(x) - P); }, x0);

            return (x1);
        }

    }

}