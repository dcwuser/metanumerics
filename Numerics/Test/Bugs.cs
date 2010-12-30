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
        public void BugBug () {

            int d = 8;
            SquareMatrix C = CyclicMatrix(d);
            Complex[] lambdas = C.Eigenvalues();
            foreach (Complex lambda in lambdas) {
                Console.WriteLine("{0} ({1} {2})", lambda, ComplexMath.Abs(lambda), ComplexMath.Arg(lambda));
            }

        }

        [TestMethod]
        public void ImproveInverseErf () {

            //foreach (double y in new double[] { 1.0E-32, 1.0E-12, 1.0E-4, 0.1, 1.0 / 3.0, 1.0 / 2.0, 2.0 / 3.0, 0.9 , 1.0 - 1.0E-4, 1.0 - 1.0E-12, 1.0 - 1.0E-32 }) {

                //Console.WriteLine("y={0}", y);
                GammaDistribution d = new GammaDistribution(5.0);
                double x = d.InverseLeftProbability(0.01);
                //double x = AdvancedMath.InverseErfc(y);
                Console.WriteLine("x={0}", x);

            //}

        }

        [TestMethod]
        public void GamAsym () {

            double a = 0.1;
            double Q = 0.01;
            Console.WriteLine("a={0} Q={1}", a, Q);

            double z0 = -Math.Log(Q * AdvancedMath.Gamma(a));
            double lz0 = Math.Log(z0);
            double z1 = z0 - (1.0 - a) * lz0;
            double z2 = z0 - (1.0 - a) * lz0 + Math.Pow(1.0 - a, 2) * (1.0 - lz0) / z0;

            GammaDistribution d = new GammaDistribution(a);
            Console.WriteLine("z0={0} Q(z0)={1}", z0, d.RightProbability(z0));
            Console.WriteLine("z1={0} Q(z1)={1}", z1, d.RightProbability(z1));
            Console.WriteLine("z2={0} Q(z2)={1}", z2, d.RightProbability(z2));

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
        public void TestArray () {

            double[] array = new double[] { 0.0, 1.0, 2.0, 3.0 };
            TestList(array);

        }

        public void TestList (IList<double> list) {

            Console.WriteLine(list.Count);
            Console.WriteLine(list.Contains(0.0));
            Console.WriteLine(list.Contains(5.0));
            Console.WriteLine(list.IndexOf(2.0));
            Console.WriteLine(list.Count);
            list[2] = 3.0;

        }

    }

}