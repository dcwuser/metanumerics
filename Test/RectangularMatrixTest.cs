using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Matrices;

using Meta.Numerics.Statistics;

namespace Test {


    [TestClass]
    public class RectangularMatrixTest {

        private TestContext testContextInstance;

        public TestContext TestContext {
            get {
                return testContextInstance;
            }
            set {
                testContextInstance = value;
            }
        }

        #region Additional test attributes
        // 
        //You can use the following additional attributes as you write your tests:
        //
        //Use ClassInitialize to run code before running the first test in the class
        //[ClassInitialize()]
        //public static void MyClassInitialize(TestContext testContext)
        //{
        //}
        //
        //Use ClassCleanup to run code after all tests in a class have run
        //[ClassCleanup()]
        //public static void MyClassCleanup()
        //{
        //}
        //
        //Use TestInitialize to run code before running each test
        //[TestInitialize()]
        //public void MyTestInitialize()
        //{
        //}
        //
        //Use TestCleanup to run code after each test has run
        //[TestCleanup()]
        //public void MyTestCleanup()
        //{
        //}
        //
        #endregion

        private RectangularMatrix GenerateRandomMatrix (int rd, int cd) {
            RectangularMatrix M = new RectangularMatrix(rd, cd);
            Random rng = new Random(1);
            for (int r = 0; r < rd; r++) {
                for (int c = 0; c < cd; c++) {
                    M[r, c] = 2.0 * rng.NextDouble() - 1.0;
                }
            }
            return (M);
        }

        [TestMethod]
        public void RectangularMatrixAccess () {

            // create a matrix via outer product
            ColumnVector cv = new ColumnVector(new double[] { 1, 2 });
            RowVector rv = new RowVector(new double[] { 3, 4, 5 });
            RectangularMatrix M = cv * rv;

            // check dimensions
            Assert.IsTrue(M.RowCount == cv.Dimension);
            Assert.IsTrue(M.ColumnCount == rv.Dimension);

            // check values
            for (int r = 0; r < M.RowCount; r++) {
                for (int c = 0; c < M.ColumnCount; c++) {
                    Assert.IsTrue(M[r, c] == cv[r] * rv[c]);
                }
            }

            // extract a column
            ColumnVector mc = M.Column(1);
            Assert.IsTrue(mc.Dimension == M.RowCount);
            for (int i = 0; i < mc.Dimension; i++) {
                Assert.IsTrue(mc[i] == M[i, 1]);
            }

            // extract a row
            RowVector mr = M.Row(1);
            Assert.IsTrue(mr.Dimension == M.ColumnCount);
            for (int i = 0; i < mr.Dimension; i++) {
                Assert.IsTrue(mr[i] == M[1, i]);
            }

            // test clone
            RectangularMatrix MC = M.Copy();
            Assert.IsTrue(MC.RowCount == M.RowCount);
            Assert.IsTrue(MC.ColumnCount == M.ColumnCount);

            // test equality of clone
            Assert.IsTrue(MC == M);
            Assert.IsFalse(MC != M);

            // test independence of clone
            MC[0, 0] += 1.0;
            Assert.IsFalse(MC == M);
            Assert.IsTrue(MC != M);

        }

        [TestMethod]
        public void RectangularMatrixEquality () {

            RectangularMatrix A = new RectangularMatrix(new double[,] {
                { 1.0, 2.0, 3.0 },
                { 4.0, 5.0, 6.0 }
            });
            RectangularMatrix AC = A.Copy();
            RectangularMatrix B = 2.0 * A;

            // Equality operator
            Assert.IsTrue(A == AC);
            Assert.IsTrue(AC == A);
            Assert.IsFalse(A == B);
            Assert.IsFalse(B == A);
            Assert.IsFalse(A == null);
            Assert.IsFalse(null == A);
            Assert.IsTrue((RectangularMatrix) null == (RectangularMatrix) null);

            // Inequality operator
            Assert.IsFalse(A != AC);
            Assert.IsFalse(AC != A);
            Assert.IsTrue(A != B);
            Assert.IsTrue(B != A);
            Assert.IsTrue(A != null);
            Assert.IsTrue(null != A);
            Assert.IsFalse((RectangularMatrix) null != (RectangularMatrix) null);

            // Equals method
            Assert.IsTrue(A.Equals(AC));
            Assert.IsFalse(A.Equals(B));
            Assert.IsFalse(A.Equals(new object()));
            Assert.IsFalse(A.Equals(null));

        }

        [TestMethod]
        public void RectangularMatrixArithmetic () {

            RectangularMatrix M = GenerateRandomMatrix(2, 3);

            RectangularMatrix MA = M + M;
            RectangularMatrix M2 = 2.0 * M;
            Assert.IsTrue(MA == M2);

            RectangularMatrix MB = MA / 2.0;
            Assert.IsTrue(MB == M);

            RectangularMatrix MS = M - M;
            RectangularMatrix M0 = 0.0 * M;
            Assert.IsTrue(MS == M0);

            RectangularMatrix MN = -M;
            RectangularMatrix MM = -1.0 * M;
            Assert.IsTrue(MN == MM);

            RectangularMatrix MT = M.Transpose();
            Assert.IsTrue(MT.RowCount == M.ColumnCount);
            Assert.IsTrue(MT.ColumnCount == M.RowCount);

            RectangularMatrix MMT = M * MT;
            Assert.IsTrue(MMT.RowCount == M.RowCount);
            Assert.IsTrue(MMT.ColumnCount == MT.ColumnCount);

        }

        [TestMethod]
        public void RectangularQRDecomposition () {

            RectangularMatrix M = GenerateRandomMatrix(30, 10);

            QRDecomposition QRD = M.QRDecomposition();
            Assert.IsTrue(QRD.RowCount == M.RowCount);
            Assert.IsTrue(QRD.ColumnCount == M.ColumnCount);

            SquareMatrix Q = QRD.QMatrix;
            Assert.IsTrue(Q.Dimension == M.RowCount);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(Q * Q.Transpose(), TestUtilities.CreateSquareUnitMatrix(Q.Dimension)));

            RectangularMatrix R = QRD.RMatrix;
            Assert.IsTrue(R.RowCount == M.RowCount);
            Assert.IsTrue(R.ColumnCount == M.ColumnCount);

            RectangularMatrix QR = Q * R;
            Assert.IsTrue(TestUtilities.IsNearlyEqual(QR, M));

        }

        private void WriteMatrix (AnyRectangularMatrix A) {
            Console.WriteLine("--");
            for (int r = 0; r < A.RowCount; r++) {
                for (int c = 0; c < A.ColumnCount; c++) {
                    Console.Write("{0}  ", A[r, c]);
                }
                Console.WriteLine();
            }

        }

        [TestMethod]
        public void RectangularMatrixNorms () {

            RectangularMatrix Z = new RectangularMatrix(3, 4);
            Assert.IsTrue(Z.OneNorm() == 0.0);
            Assert.IsTrue(Z.InfinityNorm() == 0.0);
            Assert.IsTrue(Z.FrobeniusNorm() == 0.0);

            RectangularMatrix A = GenerateRandomMatrix(4, 5);
            Assert.IsTrue(A.OneNorm() > 0.0);
            Assert.IsTrue(A.InfinityNorm() > 0.0);
            Assert.IsTrue(A.FrobeniusNorm() > 0.0);

            RectangularMatrix B = GenerateRandomMatrix(5, 6);
            Assert.IsTrue(B.OneNorm() > 0.0);
            Assert.IsTrue(B.InfinityNorm() > 0.0);
            Assert.IsTrue(B.FrobeniusNorm() > 0.0);

            // Frobenius norm is sub-multiplicative
            RectangularMatrix P = A * B;
            Assert.IsTrue(P.FrobeniusNorm() <= A.FrobeniusNorm() * B.FrobeniusNorm());

        }

        [TestMethod]
        public void BigSVD () {

            RectangularMatrix R = GenerateRandomMatrix(500, 100);

            Stopwatch s = Stopwatch.StartNew();
            SingularValueDecomposition SVD = R.SingularValueDecomposition();
            s.Stop();

            Console.WriteLine(s.Elapsed);
            Console.WriteLine(SVD.ConditionNumber);
        }

        [TestMethod]
        public void RandomRectangularSVD () {

            for (int c = 1; c < 64; c += 11) {
                Console.WriteLine(c);

                RectangularMatrix R = GenerateRandomMatrix(64, c);

                SingularValueDecomposition SVD = R.SingularValueDecomposition();

                Assert.IsTrue(SVD.RowCount == R.RowCount);
                Assert.IsTrue(SVD.ColumnCount == SVD.ColumnCount);
                Assert.IsTrue(SVD.Dimension == SVD.ColumnCount);

                SquareMatrix U = SVD.LeftTransformMatrix;
                Assert.IsTrue(U.Dimension == R.RowCount);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(U.Transpose() * U, TestUtilities.CreateSquareUnitMatrix(U.Dimension)));

                SquareMatrix V = SVD.RightTransformMatrix;
                Assert.IsTrue(V.Dimension == R.ColumnCount);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(V.Transpose() * V, TestUtilities.CreateSquareUnitMatrix(V.Dimension)));

                RectangularMatrix S = U.Transpose() * R * V;
                for (int i = 0; i < SVD.Dimension; i++) {
                    double w = SVD.SingularValue(i);
                    Assert.IsTrue(w >= 0.0);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(S[i, i], w));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(R * SVD.RightSingularVector(i), w * SVD.LeftSingularVector(i)));
                }

            }

        }

        [TestMethod]
        public void PC () {

            Random rng = new Random(1);
            double s = 1.0 / Math.Sqrt(2.0);

            MultivariateSample MS = new MultivariateSample(2);
            RectangularMatrix R = new RectangularMatrix(1000, 2);
            for (int i = 0; i < 1000; i++) {
                double r1 = 2.0 * rng.NextDouble() - 1.0;
                double r2 = 2.0 * rng.NextDouble() - 1.0;
                double x = r1 * 4.0 * s - r2 * 9.0 * s;
                double y = r1 * 4.0 * s + r2 * 9.0 * s;
                R[i, 0] = x; R[i, 1] = y;
                MS.Add(x, y);
            }

            Console.WriteLine("x {0} {1}", MS.Column(0).Mean, MS.Column(0).Variance);
            Console.WriteLine("y {0} {1}", MS.Column(1).Mean, MS.Column(1).Variance);

            Console.WriteLine("SVD");

            SingularValueDecomposition SVD = R.SingularValueDecomposition();
            for (int i = 0; i < SVD.Dimension; i++) {
                Console.WriteLine("{0} {1}", i, SVD.SingularValue(i));
                ColumnVector v = SVD.RightSingularVector(i);
                Console.WriteLine("  {0} {1}", v[0], v[1]);
            }

            Console.WriteLine("PCA");

            PrincipalComponentAnalysis PCA = MS.PrincipalComponentAnalysis();
            for (int i = 0; i < PCA.Dimension; i++) {
                PrincipalComponent PC = PCA.Components[i];
                Console.WriteLine("  {0} {1} {2} {3}", PC.Index, PC.Weight, PC.VarianceFraction, PC.CumulativeVarianceFraction);
                RowVector v = PC.NormalizedVector;
                Console.WriteLine("  {0} {1}", v[0], v[1]);
            }

            // reconstruct
            SquareMatrix U = SVD.LeftTransformMatrix;
            SquareMatrix V = SVD.RightTransformMatrix;
            double x1 = U[0, 0] * SVD.SingularValue(0) * V[0, 0] + U[0, 1] * SVD.SingularValue(1) * V[0, 1];
            Console.WriteLine("x1 = {0} {1}", x1, R[0, 0]);
            double y1 = U[0, 0] * SVD.SingularValue(0) * V[1, 0] + U[0, 1] * SVD.SingularValue(1) * V[1, 1];
            Console.WriteLine("y1 = {0} {1}", y1, R[0, 1]);
            double x100 = U[100, 0] * SVD.SingularValue(0) * V[0, 0] + U[100, 1] * SVD.SingularValue(1) * V[0, 1];
            Console.WriteLine("x100 = {0} {1}", x100, R[100, 0]);
            double y100 = U[100, 0] * SVD.SingularValue(0) * V[1, 0] + U[100, 1] * SVD.SingularValue(1) * V[1, 1];
            Console.WriteLine("y100 = {0} {1}", y100, R[100, 1]);

            ColumnVector d1 = U[0, 0] * SVD.SingularValue(0) * SVD.RightSingularVector(0) +
                U[0, 1] * SVD.SingularValue(1) * SVD.RightSingularVector(1);
            Console.WriteLine("d1 = ({0} {1})", d1[0], d1[1]);
            ColumnVector d100 = U[100, 0] * SVD.SingularValue(0) * SVD.RightSingularVector(0) +
                U[100, 1] * SVD.SingularValue(1) * SVD.RightSingularVector(1);
            Console.WriteLine("d100 = ({0} {1})", d100[0], d100[1]);

            Console.WriteLine("compare");
            MultivariateSample RS = PCA.TransformedSample();
            IEnumerator<double[]> RSE = RS.GetEnumerator();
            RSE.MoveNext();
            double[] dv1 = RSE.Current;
            Console.WriteLine("{0} {1}", dv1[0], dv1[1]);
            Console.WriteLine("{0} {1}", U[0, 0], U[0, 1]);
            RSE.Dispose();

        }

        [TestMethod]
        public void SmallSVD () {

            SquareMatrix A0 = new SquareMatrix(1);
            A0[0, 0] = 0.0;
            SingularValueDecomposition SVD0 = A0.SingularValueDecomposition();
            Console.WriteLine(SVD0.SingularValue(0));
            Assert.IsTrue(SVD0.SingularValue(0) == 0.0);

            SquareMatrix A1 = new SquareMatrix(1);
            A1[0, 0] = 1.0;
            SingularValueDecomposition SVD1 = A1.SingularValueDecomposition();
            Console.WriteLine(SVD1.SingularValue(0));
            //Assert.IsTrue(SVD1.SingularValue(0) == 1.0);

            SquareMatrix A2 = new SquareMatrix(2);
            A2[0, 0] = 0.0; A2[0, 1] = 1.0;
            A2[1, 0] = 0.0; A2[1, 1] = 1.0;
            // Singular values Sqrt(2), 0
            SingularValueDecomposition SVD2 = A2.SingularValueDecomposition();
            SquareMatrix S2 = SVD2.LeftTransformMatrix.Transpose() * A2 * SVD2.RightTransformMatrix;
            for (int i = 0; i < SVD2.Dimension; i++) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(S2[i, i], SVD2.SingularValue(i)));
            }

        }

        [TestMethod]
        public void MatrixSelfMultiplication () {

            RectangularMatrix A = GenerateRandomMatrix(3, 4);

            SymmetricMatrix AAT = A.MultiplySelfByTranspose();
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AAT, A * A.Transpose()));

            SymmetricMatrix ATA = A.MultiplyTransposeBySelf();
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ATA, A.Transpose() * A));

        }


        [TestMethod]
        public void MatrixArrayConversion () {

            // start with a .NET Array
            double[,] A = new double[,] {
                { 0, 1, 2 },
                { 3, 4, 5 }
            };

            // Convert it to a Matrix
            RectangularMatrix B = new RectangularMatrix(A);

            Assert.IsTrue(B.RowCount == A.GetLength(0));
            Assert.IsTrue(B.ColumnCount == A.GetLength(1));

            for (int r = 0; r < B.RowCount; r++) {
                for (int c = 0; c < B.ColumnCount; c++) {
                    Assert.IsTrue(B[r, c] == A[r, c]);
                }
            }

            // Convert that Matrix back to an Array
            double[,] C = B.ToArray();

            Assert.IsTrue(C.Rank == 2);
            Assert.IsTrue(C.GetLength(0) == B.RowCount);
            Assert.IsTrue(C.GetLength(1) == B.ColumnCount);

            for (int r = 0; r < B.RowCount; r++) {
                for (int c = 0; c < B.ColumnCount; c++) {
                    Assert.IsTrue(C[r, c] == A[r, c]);
                }
            }

        }

        [TestMethod]
        public void CastRectangularToSquare () {

            RectangularMatrix A = GenerateRandomMatrix(2, 3);

            RectangularMatrix AT = A.Transpose();

            RectangularMatrix ATA1 = AT * A;

            SquareMatrix ATA2 = (SquareMatrix) ATA1;
            Assert.IsTrue(ATA2 == ATA1);
       
            SymmetricMatrix ATA3 = A.MultiplyTransposeBySelf();
            Assert.IsTrue(ATA2 == ATA3);
        
            RectangularMatrix AAT1 = A * AT;

            SquareMatrix AAT2 = (SquareMatrix) AAT1;
            Assert.IsTrue(AAT2 == AAT1);

            SymmetricMatrix AAT3 = A.MultiplySelfByTranspose();
            Assert.IsTrue(AAT2 == AAT3);

        }
    }

}
