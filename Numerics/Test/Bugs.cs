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
        public void Bug7686 () {
            // This SVD failed with a IndexOutOfBoundsException because we didn't handle SVD for cols > rows and didn't check for this condition on entry.
            RectangularMatrix A = new RectangularMatrix(new double[,] {
                {-418.746, 310.726, 313.969, 1},
                {-418.746, 337.451, 229.786, 1},
                {-305.253, 321.895, 304.895, 1}
            });
            var SVD = A.SingularValueDecomposition();
        }

        [TestMethod]
        public void Bug7685 () {
            // This SVD failed with a NonconvergenceException.
            // Probably this was due to a zero appearing in the superdiagonal band in an intermediate position, since the code change to handle that eliminated the probvlem.
            SquareMatrix A = new SquareMatrix(new double[,] {
                {1.900019E+01, 0.000000E+00, 0.000000E+00, -1.385471E+02, 1.027977E+05, 1.520100E+04},
                {0.000000E+00, 1.900019E+01, 0.000000E+00, -1.026921E+05, -1.499209E+02, 1.410499E+05},
                {0.000000E+00, 0.000000E+00, 1.900019E+01, -1.499527E+04, -1.410719E+05, 0.000000E+00},
                {-1.385471E+02, -1.026921E+05, -1.499527E+04, 5.669219E+08, 1.101144E+08, -7.618976E+08},
                {1.027977E+05, -1.499209E+02, -1.410719E+05, 1.101144E+08, 1.670648E+09, 8.113594E+07},
                {1.520100E+04, 1.410499E+05, 0.000000E+00, -7.618976E+08, 8.113594E+07, 1.126334E+09}
            });
            SingularValueDecomposition SVD = A.SingularValueDecomposition();
            for (int i = 0; i < SVD.Dimension; i++) {
                Console.WriteLine(SVD.SingularValue(i));
            }
        }

        [TestMethod]
        public void Bug7684 () {

            // These values caused incomplete Beta to be called with argument outside the interval [0,1].

            double I1 = 0.9999902;
            double a1 = 0.0000434313636267175;
            double b1 = 18474.36071078790000;
            BetaDistribution D1 = new BetaDistribution(a1, b1);
            double x1 = D1.InverseLeftProbability(I1);
            Console.WriteLine("{0} {1} {2}", x1, D1.LeftProbability(x1), I1);

            double I2 = 0.9998063099306;
            double a2 = 0.00034509911609819255;
            double b2 = 6.8453983996634218;
            BetaDistribution D2 = new BetaDistribution(a2, b2);
            double x2 = D2.InverseLeftProbability(I2);
            Console.WriteLine("{0} {1} {2}", x2, D2.LeftProbability(x2), I2);

        }

        [TestMethod]
        public void Bug7213 () {

            Sample s = new Sample();
            s.Add(0.00590056, 0.00654598, 0.0066506, 0.00679065, 0.008826);
            FitResult r = WeibullDistribution.FitToSample(s);

        }

        [TestMethod]
        public void Bug7208 () {
            // this matrix has two eigenpairs with the same eigenvalue but distinct eigenvectors
            // we would report the same eigenvector for both, making the matrix of eigenvectors
            // non-invertible
            int n = 4;
            SquareMatrix matrix = new SquareMatrix(n + 1);
            double d = (3 + 2 * Math.Cos(2 * Math.PI / n));
            double w = (64 * n) / (40 - d * d) - n;
            double w1 = w / (w + n);
            double w2 = 1 / (w + n);
            matrix[0, 0] = w1;
            for (int i = 1; i < n; i++) {
                matrix[0, i + 1] = w2;
                matrix[i + 1, 0] = 3.0 / 8;
                matrix[i + 1, i + 1] = 3.0 / 8;
                matrix[i, i + 1] = 1.0 / 8;
                matrix[i + 1, i] = 1.0 / 8;
            }
            matrix[0, 1] = w2;
            matrix[1, 0] = 3.0 / 8;
            matrix[1, 1] = 3.0 / 8;
            matrix[1, n] = 1.0 / 8;
            matrix[n, 1] = 1.0 / 8;

            ComplexEigensystem ces = matrix.Eigensystem();

            Console.WriteLine("output");
            SquareMatrix V = new SquareMatrix(n + 1);
            for (int i = 0; i < n + 1; i++) {
                Console.WriteLine("#{0}: {1}", i, ces.Eigenvalue(i));
                for (int j = 0; j < n + 1; j++) {
                    V[i, j] = ces.Eigenvector(i)[j].Re;
                    Console.WriteLine(V[i, j]);
                }
            }

            SquareMatrix inv = V.Inverse();

        }

        [TestMethod]
        public void Bug6162 () {

            // When UncertianMeasurementSample.FitToPolynomial used Cholesky inversion of (A^T A), this inversion
            // would fail when roundoff errors would made the matrix non-positive-definite. We have now changed
            // to QR decomposition, which is more robust.

            //real data
            double[] X_axis = new double[] { 40270.65625, 40270.6569444444, 40270.6576388888, 40270.6583333332, 40270.6590277776,
                40270.659722222, 40270.6604166669, 40270.6611111113, 40270.6618055557, 40270.6625000001 };

            double[] Y_axis = new double[] { 246.824996948242, 246.850006103516, 245.875, 246.225006103516, 246.975006103516,
                247.024993896484, 246.949996948242, 246.875, 247.5, 247.100006103516 };

            UncertainMeasurementSample DataSet = new UncertainMeasurementSample();
            
            for (int i = 0; i < 10; i++) DataSet.Add(X_axis[i], Y_axis[i], 1);
            //for (int i = 0; i < 10; i++) DataSet.Add(X_axis[i] - 40270.0, Y_axis[i] - 247.0, 1);

            FitResult DataFit = DataSet.FitToPolynomial(3);
 
            for (int i = 0; i < DataFit.Dimension; i++)
                Console.WriteLine("a" + i.ToString() + " = " + DataFit.Parameter(i).Value);

            BivariateSample bs = new BivariateSample();
            for (int i = 0; i < 10; i++) bs.Add(X_axis[i], Y_axis[i]);
            FitResult bsFit = bs.PolynomialRegression(3);
            for (int i = 0; i < bsFit.Dimension; i++) Console.WriteLine(bsFit.Parameter(i));

        }

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

        [TestMethod]
        public void Bug5886 () {
            // the inverse CDF of the F-distribution would fail for d2 <= 2
            double d1 = 1.0;
            double d2 = 0.1;
            FisherDistribution F = new FisherDistribution(d1, d2);

            double x1 = F.InverseLeftProbability(0.6);
            Console.WriteLine(x1);
            double P = F.LeftProbability(x1);
            Console.WriteLine(P);
        }

        // Not fixing this bug; use Polynomial interpolation for this scenario instead
        //[TestMethod]
        public void Bug6392 () {
            // bug requests that we support regression with number of points equal to number
            // of fit parameters, i.e. polynomial fit
            var biSample = new BivariateSample();
            biSample.Add(0, 1);
            biSample.Add(1, -1);
            var fitResult = biSample.LinearRegression();
        }


        [TestMethod]
        public void Bug6391 () {
            // this simple PCA caused a NonConvergenceException
            var mvSample = new MultivariateSample(2);
            mvSample.Add(0, 1);
            mvSample.Add(0, -1);
            var pca = mvSample.PrincipalComponentAnalysis();
        }

 
        [TestMethod]
        public void Bug6988 () {
            // due to writing i / n instead of (double) i / n, Sample.LeftProbability was reported as 0 except for the last value
            Sample s = new Sample(0.0, 1.0, 3.0, 4.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(s.LeftProbability(2.0), 0.5));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(s.InverseLeftProbability(0.5), 2.0));
        }

    }

}