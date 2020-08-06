using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public class BugTests {

        [TestMethod]
        public void Bug56 () {

            RectangularMatrix M = new RectangularMatrix(new double[,] {
                { 44.6667,  -392.0000, -66.0000 },
                { -392.0000, 3488.0000, 504.0001 },
                { -66.0000, 504.0001, 216.0001 }
            });
            //M = M.Transpose;
            SingularValueDecomposition S = M.SingularValueDecomposition();

        }

        //[TestMethod]
        [DeploymentItem("Bug53.txt")]
        public void Bug53 () {
            // User sent data for which non-linear fit failed with Nonconvergence exception.
            List<double> x = new List<double>();
            List<double> y = new List<double>();
            using (System.IO.StreamReader reader = System.IO.File.OpenText(@"Bug53.txt")) {
                reader.ReadLine();
                while (!reader.EndOfStream) {
                    string line = reader.ReadLine();
                    string[] values = line.Split(';');
                    y.Add(Double.Parse(values[0]));
                    x.Add(Double.Parse(values[1]));
                }
            }
            NonlinearRegressionResult r = y.NonlinearRegression(
                x,
                (IReadOnlyDictionary<string, double> p, double u) => {
                    return p["0"] * (1.0 - Math.Exp(-p["1"] * (u - p["2"])));
                },
                //new Dictionary<string, double>() { { "0", 100.0 }, { "1", 1.0 }, { "2", 1.0 } }
                new Dictionary<string, double>() { {"0", 1.0 }, {"1", 1.0}, {"2", 0.0} }
            );

        }

        [TestMethod]
        [DeploymentItem("Bug8021.csv")]
        public void Bug8021 () {

            // User presented a matrix with a large number (138) of zero eigenvalues.
            // QR algorithm got tripped up on high degeneracy, but eigenvalues could
            // be revealed before QR by simple index permutation. Added code to isolate
            // these "cheap" eigenvalues before starting QR algorithm.

            int n = 276;
            SquareMatrix A = new SquareMatrix(n);
            using (System.IO.StreamReader reader = System.IO.File.OpenText(@"Bug8021.csv")) {
                int r = 0;
                while (!reader.EndOfStream) {
                    string line = reader.ReadLine();
                    string[] cells = line.Split(',');
                    for (int c = 0; c < cells.Length; c++) {
                        string cell = cells[c];
                        double value = Double.Parse(cell);
                        A[r, c] = value;
                    }
                    r++;
                }
            }

            ComplexEigendecomposition S = A.Eigendecomposition();
            foreach (ComplexEigenpair pair in S.Eigenpairs)  {
                TestUtilities.IsNearlyEigenpair(A, pair.Eigenvector, pair.Eigenvalue);
            }

            Complex[] eigenvalues = A.Eigenvalues();
            double trace = A.Trace();
            Complex sum = 0.0;
            for (int i = 0; i < eigenvalues.Length; i++) sum += eigenvalues[i];
            TestUtilities.IsNearlyEqual(trace, sum);

        }

        [TestMethod]
        public void Bug7953 () {

            // Fitting this sample to a Weibull caused a NonconvergenceException in the root finder that was used inside the fit method.
            // The underlying problem was that our equation to solve involved x^k and k ~ 2000 and (~12)^(~2000) overflows double
            // so all the quantities became Infinity and the root-finder never converged. We changed the algorithm to operate on
            // w = log x - <log x> which keeps quantities much smaller.

            Sample sample = new Sample(
                12.824, 12.855, 12.861, 12.862, 12.863,
                12.864, 12.865, 12.866, 12.866, 12.866,
                12.867, 12.867, 12.868, 12.868, 12.870,
                12.871, 12.871, 12.871, 12.871, 12.872,
                12.876, 12.878, 12.879, 12.879, 12.881
            );

            WeibullFitResult result = WeibullDistribution.FitToSample(sample);
            Console.WriteLine("{0} {1}", result.Scale, result.Shape);
            Console.WriteLine(result.GoodnessOfFit.Probability);
        }


        [TestMethod]
        public void Bug7887 () {
            // A user reported a NoncovergenceException for -0.213170584. I was able to reproduce with -0.213170585, and ultimately other values too.
            // These came from infinitely repeating Halley iterations, so changed I changed the termination criterion to not be strict equality.
            double x = -0.213170585;
            double W = AdvancedMath.LambertW(x);
            Console.WriteLine("{0} -> {1}", x, W);
        }

        [TestMethod]
        public void Bug7788 () {

            // Beta with high parameters used to have incorrect inverse probabilty,
            // which in this case resulted in zero median and consistently zero random values.

            ContinuousDistribution beta = new BetaDistribution(1713.0, 58743.0);
            Console.WriteLine(beta.Mean);
            Console.WriteLine(beta.StandardDeviation);

            Assert.IsTrue(beta.Median != 0.0);
            Assert.IsTrue(Math.Abs(beta.Median - beta.Median) <= beta.StandardDeviation);
            Random rng = new Random(1);
            for (int i = 0; i < 10; i++) {
                Assert.IsTrue(beta.GetRandomValue(rng) != 0.0);
            }

        }

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
            // Probably this was due to a zero appearing in the super-diagonal band in an intermediate position, since the code change to handle that eliminated the problem.
            SquareMatrix A = new SquareMatrix(new double[,] {
                {1.900019E+01, 0.000000E+00, 0.000000E+00, -1.385471E+02, 1.027977E+05, 1.520100E+04},
                {0.000000E+00, 1.900019E+01, 0.000000E+00, -1.026921E+05, -1.499209E+02, 1.410499E+05},
                {0.000000E+00, 0.000000E+00, 1.900019E+01, -1.499527E+04, -1.410719E+05, 0.000000E+00},
                {-1.385471E+02, -1.026921E+05, -1.499527E+04, 5.669219E+08, 1.101144E+08, -7.618976E+08},
                {1.027977E+05, -1.499209E+02, -1.410719E+05, 1.101144E+08, 1.670648E+09, 8.113594E+07},
                {1.520100E+04, 1.410499E+05, 0.000000E+00, -7.618976E+08, 8.113594E+07, 1.126334E+09}
            });
            SingularValueDecomposition SVD = A.SingularValueDecomposition();
            foreach (SingularValueContributor contributor in SVD.Contributors) {
                Console.WriteLine(contributor.SingularValue);
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
            WeibullFitResult r = WeibullDistribution.FitToSample(s);

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

            ComplexEigendecomposition ces = matrix.Eigendecomposition();

            SquareMatrix V = new SquareMatrix(n + 1);
            for (int i = 0; i < n + 1; i++) {
                for (int j = 0; j < n + 1; j++) {
                    V[i, j] = ces.Eigenpairs[i].Eigenvector[j].Re;
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

            UncertainMeasurementFitResult DataFit = DataSet.FitToPolynomial(3);

            BivariateSample bs = new BivariateSample();
            for (int i = 0; i < 10; i++) bs.Add(X_axis[i], Y_axis[i]);
            PolynomialRegressionResult bsFit = bs.PolynomialRegression(3);
            foreach (Parameter p in bsFit.Parameters) Console.WriteLine(p);

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


        [TestMethod]
        public void Bug10 () {

            // Fisher exact test didn't give same probability when rows were permuted.
            // To compute the Fisher exact probability, we iterate over all contingency tables
            // with the same marginal totals and count their probability if it is less than or
            // equal to the probability of the observed matrix.
            // In particular for the symmetric case, the "opposite" table has the exact
            // same probability. But there is floating point noise, so sometimes it's
            // calculated probability is infinitesimally larger and isn't counted.
            // To fix this, we special-case the symmetric case. 

            ContingencyTable t1 = new ContingencyTable(new int[,] {
                { 18, 16 }, { 12, 14 }
            });

            ContingencyTable t2 = new ContingencyTable(new int[,] {
                { 12, 14 }, { 18, 16 }
            });

            Assert.IsTrue(TestUtilities.IsNearlyEqual(t1.Binary.FisherExactTest().Probability, t2.Binary.FisherExactTest().Probability));

        }

        [TestMethod]
        public void Bug2 () {

            // Integrate failed with a NullReference exception when all NaNs were returned.
            // Changed logic to ignore NaNs.

            Func<double, double> integrand = x => Double.NaN;

            IntegrationResult result = FunctionMath.Integrate(integrand, 0.0, 1.0);
            Assert.IsTrue(result.Estimate.Value == 0.0);

        } 
    }

}