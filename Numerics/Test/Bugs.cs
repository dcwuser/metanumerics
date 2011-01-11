using System;
using System.Collections.Generic;
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

            //ChiSquaredDistribution d = new ChiSquaredDistribution(4);
            //double x = d.InverseLeftProbability(1.0E-10);
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
            double x0 = FunctionMath.FindZero(delegate(double x) { return (x * x * x); }, 1.0);
            Assert.IsTrue(Math.Abs(x0) < TestUtilities.TargetPrecision);
        }

    }

}