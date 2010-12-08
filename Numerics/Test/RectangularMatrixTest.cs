using Meta.Numerics.Matrices;

using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Test {


    [TestClass()]
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
                    M[r,c] = 2.0*rng.NextDouble() - 1.0;
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
            RectangularMatrix MC = M.Clone();
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
        public void RectangularMatrixArithmetic () {

            RectangularMatrix M = GenerateRandomMatrix(2, 3);

            RectangularMatrix MA = M + M;
            RectangularMatrix M2 = 2.0 * M;
            Assert.IsTrue(MA == M2);

            RectangularMatrix MS = M - M;
            RectangularMatrix M0 = 0.0 * M;
            Assert.IsTrue(MS == M0);

            RectangularMatrix MT = M.Transpose();
            Assert.IsTrue(MT.RowCount == M.ColumnCount);
            Assert.IsTrue(MT.ColumnCount == M.RowCount);

            RectangularMatrix MM = M * MT;
            Assert.IsTrue(MM.RowCount == M.RowCount);
            Assert.IsTrue(MM.ColumnCount == MT.ColumnCount);

        }

        [TestMethod]
        public void RectangularQRDecomposition () {

            RectangularMatrix M = GenerateRandomMatrix(30, 10);

            QRDecomposition QRD = M.QRDecompose();
            Assert.IsTrue(QRD.RowCount == M.RowCount);
            Assert.IsTrue(QRD.ColumnCount == M.ColumnCount);

            SquareMatrix Q = QRD.QMatrix();
            Assert.IsTrue(Q.Dimension == M.RowCount);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(Q * Q.Transpose(), TestUtilities.CreateSquareUnitMatrix(Q.Dimension)));

            RectangularMatrix R = QRD.RMatrix();
            Assert.IsTrue(R.RowCount == M.RowCount);
            Assert.IsTrue(R.ColumnCount == M.ColumnCount);

            RectangularMatrix QR = Q * R;
            Assert.IsTrue(TestUtilities.IsNearlyEqual(QR, M));

        }

        private void WriteMatrix (RectangularMatrixBase A) {
            Console.WriteLine("--");
            for (int r = 0; r < A.RowCount; r++) {
                for (int c = 0; c < A.ColumnCount; c++) {
                    Console.Write("{0}  ", A[r, c]);
                }
                Console.WriteLine();
            }

        }
    }
}
