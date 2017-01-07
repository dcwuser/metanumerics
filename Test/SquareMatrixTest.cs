using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Test {    
    
    [TestClass]
    public class SquareMatrixTest {

        public static void PrintMatrix (AnyRectangularMatrix M) {
            for (int r = 0; r < M.RowCount; r++) {
                for (int c = 0; c < M.ColumnCount; c++) {
                    Console.Write("{0,12:g8} ", M[r, c]);
                }
                Console.WriteLine();
            }
            Console.WriteLine("--");
        }

        private static SquareMatrix CreateSquareRandomMatrix (int n) {
            return (CreateSquareRandomMatrix(n, 1));
        }

        private static SquareMatrix CreateSquareRandomMatrix (int n, int seed) {
            SquareMatrix M = new SquareMatrix(n);
            Random rng = new Random(seed);
            for (int r = 0; r < n; r++) {
                for (int c = 0; c < n; c++) {
                    M[r, c] = 2.0 * rng.NextDouble() - 1.0;
                }
            }
            return (M);
        }

        private static SquareMatrix CreateVandermondeMatrix (int n) {
            double[] x = new double[n];
            for (int i = 0; i < n; i++) {
                x[i] = i;
            }
            return (CreateVandermondeMatrix(x));
        }

        private static SquareMatrix CreateVandermondeMatrix (double[] x) {
            int d = x.Length;
            SquareMatrix V = new SquareMatrix(d);
            for (int r = 0; r < d; r++) {
                double v = 1.0;
                for (int c = 0; c < d; c++) {
                    V[r, c] = v;
                    v = v * x[r];
                }
            }
            return (V);
        }

        private static SquareMatrix CreateCirculantMatrix (double[] x) {
            int d = x.Length;
            SquareMatrix A = new SquareMatrix(d);
            for (int c = 0; c < d; c++) {
                for (int i = 0; i < d; i++) {
                    int r = (c + i) % d;
                    A[r, c] = x[i];
                }
            }
            return (A);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void SquareMatrixInvalidDimensionTest () {
            SquareMatrix M = new SquareMatrix(0);
        }

        [TestMethod]
        public void SquareMatrixAccess () {
            double[] x = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
            SquareMatrix V = CreateVandermondeMatrix(x);

            // check nullity
            Assert.IsTrue(V != null);

            // check dimension
            Assert.IsTrue(V.Dimension == x.Length);

            // check column
            ColumnVector c = V.Column(1);
            Assert.IsTrue(c.Dimension == x.Length);
            for (int i = 0; i < c.Dimension; i++) {
                Assert.IsTrue(c[i] == V[i, 1]);
                Assert.IsTrue(c[i] == x[i]);
            }

            // check row
            RowVector r = V.Row(0);
            Assert.IsTrue(r.Dimension == x.Length);
            for (int i = 0; i < r.Dimension; i++) {
                Assert.IsTrue(r[i] == V[0, i]);
                Assert.IsTrue(r[i] == 1.0);
            }

            // check clone
            SquareMatrix VC = V.Copy();
            Assert.IsTrue(VC == V);
            Assert.IsFalse(VC != V);

            // check independence of clone
            VC[0, 0] += 1.0;
            Assert.IsFalse(VC == V);
            Assert.IsTrue(VC != V);

        }

        [TestMethod]
        public void SquareMatrixArithmetic () {

            SquareMatrix M = CreateSquareRandomMatrix(5);

            // Addition is same as multiplication by two
            SquareMatrix MA = M + M;
            SquareMatrix M2 = 2.0 * M;
            Assert.IsTrue(MA == M2);

            // Division by two returns us to original
            SquareMatrix MB = MA / 2.0;
            Assert.IsTrue(MB == M);

            // Subraction of self same as multiplication by zero
            SquareMatrix MS = M - M;
            SquareMatrix M0 = 0.0 * M;
            Assert.IsTrue(MS == M0);

            // Negation is same as multiplication by negative one
            SquareMatrix MN = -M;
            SquareMatrix MM = -1.0 * M;
            Assert.IsTrue(MN == MM);

            // check transpose
            SquareMatrix MT = M.Transpose();
            Assert.IsTrue(MT != M);

            // matrix multiplication is not ableian
            SquareMatrix MMT = M * MT;
            SquareMatrix MTM = MT * M;
            Assert.IsFalse(MMT == MTM);

            // check that transpose of transpose is original
            SquareMatrix MTT = MT.Transpose();
            Assert.IsTrue(MTT == M);

        }

        [TestMethod]
        public void SquareVandermondeMatrixInverse () {
            for (int d = 1; d <= 8; d++) {
                SquareMatrix I = TestUtilities.CreateSquareUnitMatrix(d);
                SquareMatrix H = CreateVandermondeMatrix(d);
                DateTime start = DateTime.Now;
                SquareMatrix HI = H.Inverse();
                DateTime finish = DateTime.Now;
                Console.WriteLine("d={0} t={1} ms", d, (finish - start).Milliseconds);
                PrintMatrix(HI);
                PrintMatrix(H * HI);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(H * HI, I));
            }
        }

        [TestMethod]
        public void SquareRandomMatrixInverse () {
            for (int d = 1; d <= 100; d=d+11) {
                SquareMatrix I = TestUtilities.CreateSquareUnitMatrix(d);
                SquareMatrix M = CreateSquareRandomMatrix(d, 1);
                SquareMatrix MI = M.Inverse();
                Assert.IsTrue(TestUtilities.IsNearlyEqual(M * MI, I));
            }
        }


        [TestMethod]
        public void SquareRandomMatrixQRDecomposition () {
            for (int d = 1; d <= 100; d += 11) {

                Console.WriteLine("d={0}", d);

                SquareMatrix M = CreateSquareRandomMatrix(d);

                // QR decompose the matrix
                SquareQRDecomposition QRD = M.QRDecomposition();

                // the dimension should be right
                Assert.IsTrue(QRD.Dimension == M.Dimension);

                // test that the decomposition works
                SquareMatrix Q = QRD.QMatrix;
                SquareMatrix R = QRD.RMatrix;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(Q * R, M));

                // check that the inverse works
                SquareMatrix MI = QRD.Inverse();
                SquareMatrix I = TestUtilities.CreateSquareUnitMatrix(d);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(M * MI, I));

                // test that a solution works
                ColumnVector t = new ColumnVector(d);
                for (int i = 0; i < d; i++) t[i] = i;
                ColumnVector s = QRD.Solve(t);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(M * s, t));

            }
        }

        [TestMethod]
        public void SquareRandomMatrixLUDecomposition () {
            for (int d = 1; d <= 256; d += 11) {

                Console.WriteLine("d={0}", d);

                SquareMatrix M = CreateSquareRandomMatrix(d);

                // LU decompose the matrix
                //Stopwatch sw = Stopwatch.StartNew();
                LUDecomposition LU = M.LUDecomposition();
                //sw.Stop();
                //Console.WriteLine(sw.ElapsedMilliseconds);

                Assert.IsTrue(LU.Dimension == d);

                // test that the decomposition works
                SquareMatrix P = LU.PMatrix();
                SquareMatrix L = LU.LMatrix();
                SquareMatrix U = LU.UMatrix();
                Assert.IsTrue(TestUtilities.IsNearlyEqual(P * M, L * U));

                // check that the inverse works
                SquareMatrix MI = LU.Inverse();
                SquareMatrix I = TestUtilities.CreateSquareUnitMatrix(d);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(M * MI, I));

                // test that a solution works
                ColumnVector t = new ColumnVector(d);
                for (int i = 0; i < d; i++) t[i] = i;
                ColumnVector s = LU.Solve(t);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(M * s, t));

            }
        }

        [TestMethod]
        public void SquareVandermondeMatrixLUDecomposition () {
            // fails now for d = 8 because determinant slightly off
            for (int d = 1; d <= 7; d++) {

                Console.WriteLine("d={0}", d);

                double[] x = new double[d];
                for (int i = 0; i < d; i++) {
                    x[i] = i;
                }
                double det = 1.0;
                for (int i = 0; i < d; i++) {
                    for (int j = 0; j < i; j++) {
                        det = det * (x[i] - x[j]);
                    }
                }

                // LU decompose the matrix
                SquareMatrix V = CreateVandermondeMatrix(d);
                LUDecomposition LU = V.LUDecomposition();

                // test that the decomposition works
                SquareMatrix P = LU.PMatrix();
                SquareMatrix L = LU.LMatrix();
                SquareMatrix U = LU.UMatrix();
                Assert.IsTrue(TestUtilities.IsNearlyEqual(P * V, L * U));

                // check that the determinant agrees with the analytic expression
                Console.WriteLine("det {0} {1}", LU.Determinant(), det);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(LU.Determinant(), det));

                // check that the inverse works
                SquareMatrix VI = LU.Inverse();
                //PrintMatrix(VI);
                //PrintMatrix(V * VI);
                SquareMatrix I = TestUtilities.CreateSquareUnitMatrix(d);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(V * VI, I));

                // test that a solution works
                ColumnVector t = new ColumnVector(d);
                for (int i = 0; i < d; i++) t[i] = 1.0;
                ColumnVector s = LU.Solve(t);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(V * s, t));
                
            }
        }

        // eigensystem of Vandermonde matrices

        [TestMethod]
        public void SquareVandermondeMatrixEigenvalues () {
            for (int d = 1; d <= 8; d++) {
                SquareMatrix H = CreateVandermondeMatrix(d);
                double tr = H.Trace();
                ComplexEigensystem E = H.Eigensystem();
                double sum = 0.0;
                for (int i = 0; i < d; i++) {
                    Console.WriteLine(E.Eigenvalue(i));
                    double e = E.Eigenvalue(i).Re;
                    sum += e;
                    Complex[] vc = E.Eigenvector(i);
                    ColumnVector v = new ColumnVector(d);
                    for (int j = 0; j < d; j++) {
                        //Console.WriteLine("  {0}", vc[j]);
                        v[j] = vc[j].Re;
                    }
                    ColumnVector Hv = H * v;
                    ColumnVector ev = e * v;
                    Assert.IsTrue(TestUtilities.IsNearlyEigenpair(H, v, e));
                }
                Assert.IsTrue(TestUtilities.IsNearlyEqual(tr, sum));
            }
        }


        [TestMethod]
        public void SquareRandomMatrixEigenvalues () {
            for (int d = 1; d <= 67; d = d + 11) {
                SquareMatrix I = TestUtilities.CreateSquareUnitMatrix(d);
                SquareMatrix M = CreateSquareRandomMatrix(d, d);
                double tr = M.Trace();
                DateTime start = DateTime.Now;
                ComplexEigensystem E = M.Eigensystem();
                DateTime finish = DateTime.Now;
                Console.WriteLine("d={0} t={1} ms", d, (finish - start).Milliseconds);
                Assert.IsTrue(E.Dimension == d);
                Complex[] es = new Complex[d];
                for (int i = 0; i < d; i++) {
                    es[i] = E.Eigenvalue(i);
                    Complex[] v = E.Eigenvector(i);
                    Assert.IsTrue(TestUtilities.IsNearlyEigenpair(M, v, es[i]));
                }
                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(es, tr));
            }
        }

        [TestMethod]
        public void SquareUnitMatrixLUDecomposition () {
            for (int d = 1; d <= 10; d++) {
                SquareMatrix I = TestUtilities.CreateSquareUnitMatrix(d);
                Assert.IsTrue(I.Trace() == d);
                LUDecomposition LU = I.LUDecomposition();
                Assert.IsTrue(LU.Determinant() == 1.0);
                SquareMatrix II = LU.Inverse();
                Assert.IsTrue(TestUtilities.IsNearlyEqual(II, I));
            }
        }

        [TestMethod]
        public void SquareUnitMatrixEigensystem () {
            int d = 3;
            SquareMatrix I = TestUtilities.CreateSquareUnitMatrix(d);
            ComplexEigensystem E = I.Eigensystem();
            Assert.IsTrue(E.Dimension == d);
            for (int i = 0; i < d; i++) {
                Complex val = E.Eigenvalue(i);
                Assert.IsTrue(val == 1.0);
                Complex[] vec = E.Eigenvector(i);
                for (int j = 0; j < d; j++) {
                    if (i == j) {
                        Assert.IsTrue(vec[j] == 1.0);
                    } else {
                        Assert.IsTrue(vec[j] == 0.0);
                    }
                }
            }
        }

        [TestMethod]
        public void SquareMatrixDifficultEigensystem () {
            // this requires > 30 iterations to converge, looses more accuracy than it should
            // this example throws a NonConvergenceException without the ad hoc shift but works once the ad hoc shift is included
            SquareMatrix M = new SquareMatrix(4);
            M[0,1] = 2;
            M[0,3] = -1;
            M[1,0] = 1;
            M[2,1] = 1;
            M[3,2] = 1;
            // The eigenvalues of this matrix are -1, -1, 1, 1.
            // There are only two distinct eigenvectors: (1, -1, 1, -1) with eigenvalue -1, and (1, 1, 1, 1) with eigenvalue 1.
            ComplexEigensystem E = M.Eigensystem();
            for (int i = 0; i < E.Dimension; i++) {
                Console.WriteLine(E.Eigenvalue(i));
                for (int j = 0; j < E.Dimension; j++) {
                    Console.Write("{0} ", E.Eigenvector(i)[j]);
                }
                Console.WriteLine();
                Assert.IsTrue(TestUtilities.IsNearlyEigenpair(M, E.Eigenvector(i), E.Eigenvalue(i)));
            }
        }


        [TestMethod]
        public void KnownEigenvalues () {

            // This example is presented in http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf
            // It has eigenvalues 1 \pm 2i, 3, 4, 5 \pm 6i

            double[,] souce = new double[,] {
                { 7, 3, 4, -11, -9, -2},
                {-6, 4, -5, 7, 1, 12},
                {-1, -9, 2, 2, 9, 1},
                {-8, 0, -1, 5, 0, 8},
                {-4, 3, -5, 7, 2, 10},
                {6, 1, 4, -11, -7, -1}
            };

            SquareMatrix A = new SquareMatrix(6);
            for (int r = 0; r < souce.GetLength(0); r++) {
                for (int c = 0; c < souce.GetLength(1); c++) {
                    A[r, c] = souce[r, c];
                }
            }

            Complex[] eigenvalues = A.Eigenvalues();

            foreach (Complex eigenvalue in eigenvalues) {
                Console.WriteLine(eigenvalue);
            }
        }

        [TestMethod]
        public void DegenerateEigenvalues () {

            double[,] souce = new double[,] {
                {1, 0, 1},
                {0, 2, 0},
                {1, 0, 1}
            };
            SquareMatrix A = new SquareMatrix(3);
            for (int r = 0; r < souce.GetLength(0); r++) {
                for (int c = 0; c < souce.GetLength(1); c++) {
                    A[r, c] = souce[r, c];
                }
            }

            ComplexEigensystem e = A.Eigensystem();


        }


        [TestMethod]
        public void DifficultEigenvalue () {

            // This is from a paper describing difficult eigenvalue problems.
            // https://www.mathworks.com/company/newsletters/news_notes/pdf/sum95cleve.pdf

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
        public void SquareMatrixNotInvertable () {

            SquareMatrix A = new SquareMatrix(2);
            A[1, 1] = 1.0;

            // Inverting should throw
            try {
                A.Inverse();
                Assert.IsTrue(false);
            } catch (DivideByZeroException) {

            }

            // LU Decomposing should throw
            try {
                A.LUDecomposition();
                Assert.IsTrue(false);
            } catch (DivideByZeroException) {

            }

            // SVD should succeed, and give infinite condition number
            SingularValueDecomposition SVD = A.SingularValueDecomposition();
            Assert.IsTrue(Double.IsInfinity(SVD.ConditionNumber));

        }

        [TestMethod]
        public void SquareMatrixStochasticEigensystem () {

            // this is a simplifed form of a Markov matrix that arose in the Monopoly problem
            // and failed to converge

            int n = 12;
            // originally failed for n=12 due to lack of 2x2 detection at top
            // now tested for n up to 40, but n=14 appears to still be a problem,
            // probably due to a six-times degenerate eigenvalue

            SquareMatrix R = new SquareMatrix(n);
            for (int c = 0; c < n; c++) {
                R[(c + 2) % n, c] = 1.0 / 36.0;
                R[(c + 3) % n, c] = 2.0 / 36.0;
                R[(c + 4) % n, c] = 3.0 / 36.0;
                R[(c + 5) % n, c] = 4.0 / 36.0;
                R[(c + 6) % n, c] = 5.0 / 36.0;
                R[(c + 7) % n, c] = 6.0 / 36.0;
                R[(c + 8) % n, c] = 5.0 / 36.0;
                R[(c + 9) % n, c] = 4.0 / 36.0;
                R[(c + 10) % n, c] = 3.0 / 36.0;
                R[(c + 11) % n, c] = 2.0 / 36.0;
                R[(c + 12) % n, c] = 1.0 / 36.0;
            }

            //Complex[] v = R.Eigenvalues();
            ComplexEigensystem E = R.Eigensystem();

            for (int i = 0; i < E.Dimension; i++) {
                Console.WriteLine(E.Eigenvalue(i));
                Assert.IsTrue(TestUtilities.IsNearlyEigenpair(R, E.Eigenvector(i), E.Eigenvalue(i)));
            }

        }

        [TestMethod]
        public void CirculantEigenvalues () {

            int n = 12;

            double[] x = new double[n];
            for (int i = 0; i < x.Length; i++) x[i] = 1.0 / (i + 1);
            SquareMatrix A = CreateCirculantMatrix(x);

            Complex[] u = RootsOfUnity(n);
            Complex[] v = new Complex[n];
            for (int i = 0; i < v.Length; i++) {
                for (int j = 0; j < n; j++) {
                    v[i] += x[j] * u[(i * (n - j)) % n];
                }
            }

            Complex[][] w = new Complex[n][];
            for (int i = 0; i < n; i++) {
                w[i] = new Complex[n];
                for (int j = 0; j < n; j++) {
                    w[i][j] = u[i * j % n];
                }
            }

            Complex[] eValues = A.Eigenvalues();
            Complex eProduct = 1.0;
            foreach (Complex eValue in eValues) {
                eProduct *= eValue;
            }

            // v and eValues should be equal. By inspection they are,
            // but how to verify this given floating point jitter?

            // Verify that eigenvalue product equals determinant.
            SquareQRDecomposition QR = A.QRDecomposition();
            double det = QR.Determinant();

            Assert.IsTrue(TestUtilities.IsNearlyEqual(eProduct, det));
        }

        public Complex[] RootsOfUnity (int n) {
            Complex[] u = new Complex[n];
            for (int i = 0; i < n; i++) {
                u[i] = new Complex(MoreMath.Cos(2.0 * Math.PI * i / n), MoreMath.Sin(2.0 * Math.PI * i / n));
            }
            return (u);
        }

        public static void PrintMatrix (double[,] A) {
            for (int r = 0; r < A.GetLength(0); r++) {
                for (int c = 0; c < A.GetLength(1); c++) {
                    Console.Write("  {0}", A[r, c]);
                }
                Console.WriteLine();
            }
            Console.WriteLine();
        }

        public void TimedEigenvalues () {

            int n = 200;
            Random rng = new Random(1);

            SquareMatrix A = new SquareMatrix(n);
            for (int r = 0; r < n; r++) {
                for (int c = 0; c < n; c++) {
                    A[r, c] = rng.NextDouble();
                }
            }

            Stopwatch s = Stopwatch.StartNew();
            Complex[] eigenvalues = A.Eigenvalues();
            s.Stop();
            Console.WriteLine(s.ElapsedMilliseconds);

        }

        [TestMethod]
        public void SquareMatrixNorms () {

            SquareMatrix Z = new SquareMatrix(3);
            Assert.IsTrue(Z.OneNorm() == 0.0);
            Assert.IsTrue(Z.InfinityNorm() == 0.0);
            Assert.IsTrue(Z.FrobeniusNorm() == 0.0);

            SquareMatrix A = CreateSquareRandomMatrix(4);
            Assert.IsTrue(A.OneNorm() > 0.0);
            Assert.IsTrue(A.InfinityNorm() > 0.0);
            Assert.IsTrue(A.FrobeniusNorm() > 0.0);

            SquareMatrix B = CreateVandermondeMatrix(4);
            Assert.IsTrue(B.OneNorm() > 0.0);
            Assert.IsTrue(B.InfinityNorm() > 0.0);
            Assert.IsTrue(B.FrobeniusNorm() > 0.0);

            SquareMatrix S = A + B;
            Assert.IsTrue(S.OneNorm() <= A.OneNorm() + B.OneNorm());
            Assert.IsTrue(S.InfinityNorm() <= A.InfinityNorm() + B.InfinityNorm());
            Assert.IsTrue(S.FrobeniusNorm() <= A.FrobeniusNorm() + B.FrobeniusNorm());

            SquareMatrix P = A * B;
            Assert.IsTrue(P.FrobeniusNorm() <= A.FrobeniusNorm() * B.FrobeniusNorm());
            // Frobenium norm is sub-multiplicative

        }

        [TestMethod]
        public void HilbertMatrixSVD () {

            int n = 4;
            SquareMatrix H = new SquareMatrix(n);
            for (int r = 0; r < n; r++) {
                for (int c = 0; c < n; c++) {
                    H[r, c] = 1.0 / (r + c + 1);
                }
            }

            SingularValueDecomposition SVD = H.SingularValueDecomposition();
            for (int i = 0; i < n; i++) {
                Console.WriteLine(SVD.SingularValue(i));
            }

            SquareMatrix S = SVD.LeftTransformMatrix().Transpose() * H * SVD.RightTransformMatrix();
            for (int i = 0; i < SVD.Dimension; i++) {
                Console.WriteLine(S[i, i]);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(S[i, i], SVD.SingularValue(i)));
            }

        }

        [TestMethod]
        public void CompanionMatrixEigenvalues () {


            SquareMatrix CM = new SquareMatrix(3);
            CM[0, 2] = -1.0;
            CM[1, 0] = 1.0;
            CM[1, 2] = -3.0;
            CM[2, 1] = 1.0;
            CM[2, 2] = -3.0;
            ComplexEigensystem EM = CM.Eigensystem();
            for (int i = 0; i < EM.Dimension; i++) {
                Console.WriteLine("! {0}", EM.Eigenvalue(i));
            }


            for (int d = 2; d <= 8; d++) {

                Console.WriteLine(d);

                SquareMatrix C = new SquareMatrix(d);
                for (int r = 1; r < d; r++) {
                    C[r, r - 1] = 1.0;
                }
                for (int r = 0; r < d; r++) {
                    C[r, d - 1] = -AdvancedIntegerMath.BinomialCoefficient(d, r);
                }

                ComplexEigensystem e = C.Eigensystem();
                for (int i = 0; i < e.Dimension; i++) {
                    Console.WriteLine("!! {0}", e.Eigenvalue(i));
                }


            }

        }


        [TestMethod]
        public void MatrixPowerSpecialCases () {

            SquareMatrix A = CreateSquareRandomMatrix(3);

            // Zero power is identity
            SquareMatrix A0 = A.Power(0);

            // One power is self
            SquareMatrix A1 = A.Power(1);
            Assert.IsTrue(A == A1);

            // Two powers same as multiplying by self
            SquareMatrix A2 = A.Power(2);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(A2, A * A));

        }

        [TestMethod]
        public void MatrixPeriodTest () {

            // Period 1 (idempotent) matrix
            SquareMatrix A = new SquareMatrix(new double[,] {
                { 1.0, 1.0 },
                { 0.0, 0.0 }
            });

            SquareMatrix A2 = A.Power(2);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(A2, A));

            // Period 2 (self-inverse) matrix
            SquareMatrix B = new SquareMatrix(2);
            B[0, 0] = 0.0; B[0, 1] = 1.0;
            B[1, 0] = 1.0; B[1, 1] = 0.0;

            SquareMatrix B3 = B.Power(3);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(B3, B));

            // Period 3 matrix
            // Use the cyclic permutation matrix of order 3
            SquareMatrix C = new SquareMatrix(3);
            C[0, 1] = 1.0; C[1, 2] = 1.0; C[2, 0] = 1.0;

            SquareMatrix C4 = C.Power(4);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(C4, C));

        }

        [TestMethod]
        public void SvdOfRankOneMatrix () {

            // Create a rank-1 matrix
            ColumnVector v = new ColumnVector(1.0, 2.0, 3.0, 4.0);
            RectangularMatrix A = v * v.Transpose();

            SingularValueDecomposition SVD = A.SingularValueDecomposition();

            // The rank should be 1
            Assert.IsTrue(SVD.Rank == 1);

            // Only the first singular value should be non-zero
            Console.WriteLine(SVD.SingularValue(0));
            for (int i = 1; i < SVD.Dimension; i++) {
                Console.WriteLine(SVD.SingularValue(i));
                Assert.IsTrue(SVD.SingularValue(i) < TestUtilities.TargetPrecision * SVD.SingularValue(0));
            }

        }

        [TestMethod]
        public void SquareMatrixSVD () {

            for (int d = 4; d < 50; d += 7) {

                SquareMatrix A = CreateSquareRandomMatrix(d, d);

                SingularValueDecomposition SVD = A.SingularValueDecomposition();

                Assert.IsTrue(SVD.Dimension == A.Dimension);

                ColumnVector x = new ColumnVector(d);
                for (int i = 0; i < d; i++) x[i] = i;
                ColumnVector b = A * x;

                ColumnVector b1 = SVD.LeftTransformMatrix().Transpose() * b;
                for (int i = 0; i < SVD.Rank; i++) b1[i] = b1[i] / SVD.SingularValue(i);
                ColumnVector b2 = SVD.RightTransformMatrix() * b1;

                ColumnVector y = SVD.Solve(b);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(x, y));

            }

        }

    }

}
