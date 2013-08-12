using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Matrices;

namespace Test {

    [TestClass()]
    public class TridiagonalMatrixTest {


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

        public TridiagonalMatrix CreateRandomTridiagonalMatrix (int n) {
            return (CreateRandomTridiagonalMatrix(n, new Random(1)));
        }

        public TridiagonalMatrix CreateRandomTridiagonalMatrix (int n, Random rng) {
            TridiagonalMatrix T = new TridiagonalMatrix(n);
            for (int i = 0; i < (n-1); i++) {
                T[i, i] = 2.0 * rng.NextDouble() - 1.0;
                T[i + 1, i] = 2.0 * rng.NextDouble() - 1.0;
                T[i, i + 1] = 2.0 * rng.NextDouble() - 1.0;
            }
            T[n - 1, n - 1] = 2.0 * rng.NextDouble() - 1.0;
            return (T);
        }

        [TestMethod]
        public void TridiagonalMatrixAccessTest () {

            TridiagonalMatrix T = CreateRandomTridiagonalMatrix(4);

            // check nullity
            Assert.IsTrue(T != null);

            // check dimension
            Assert.IsTrue(T.Dimension == 4);

            // check off-diagonal element
            Assert.IsTrue(T[0, 3] == 0.0);

        }

        [TestMethod]
        [ExpectedException(typeof(InvalidOperationException))]
        public void TridiagonalMatrixInvalidSetTest () {
            TridiagonalMatrix T = new TridiagonalMatrix(4);
            T[0, 3] = 1.0;
        }

        [TestMethod]
        public void TridiagonalMatrixArithmeticTest () {

            TridiagonalMatrix T = CreateRandomTridiagonalMatrix(4);

            //SquareMatrixTest.PrintMatrix(T);
            TridiagonalMatrix TA = T + T;
            //SquareMatrixTest.PrintMatrix(TA);
            TridiagonalMatrix T2 = 2.0 * T;
            //SquareMatrixTest.PrintMatrix(T2);
            Assert.IsTrue(TA == T2);

            TridiagonalMatrix TM = T - T;
            TridiagonalMatrix T0 = 0.0 * T;
            Assert.IsTrue(TM == T0);

        }

        [TestMethod]
        public void TridiagonalMatrixDuplicationTest () {

            TridiagonalMatrix T = CreateRandomTridiagonalMatrix(4);

            // check that clone is clone
            TridiagonalMatrix TC = T.Copy();
            Assert.IsTrue(TC == T);

            // check that clone is independend
            TC[0, 0] += 1.0;
            Assert.IsTrue(TC != T);

            // check that transpose of transpose is original
            TridiagonalMatrix TT = T.Transpose();
            Assert.IsTrue(TT != T);
            TridiagonalMatrix TTT = TT.Transpose();
            Assert.IsTrue(TTT == T);

            // check that transpose is independent
            TTT[0, 0] += 1.0;
            Assert.IsTrue(TTT != T);

        }

        [TestMethod]
        public void TridiagonalMatrixLUDecompositionTest () {

            for (int d = 3; d < 100; d = 2 * d) {

                Console.WriteLine("d={0}", d);

                TridiagonalMatrix T = CreateRandomTridiagonalMatrix(d);
                Assert.IsTrue(T.Dimension == d);


                TridiagonalLUDecomposition LU = T.LUDecomposition();
                Assert.IsTrue(LU.Dimension == d);

                // check determinant
                Assert.IsTrue(TestUtilities.IsNearlyEqual(LU.Determinant(), T.Determinant()));

                //SquareMatrix P = LU.PMatrix();
                //SquareMatrix L = LU.LMatrix();
                //SquareMatrix U = LU.UMatrix();

                //SquareMatrixTest.PrintMatrix(T);
                //SquareMatrixTest.PrintMatrix(L);
                //SquareMatrixTest.PrintMatrix(U);
                //SquareMatrixTest.PrintMatrix(L * U);


                // check solution to decomposition
                ColumnVector b = new ColumnVector(d);
                for (int i = 0; i < d; i++) {
                    b[i] = i;
                }
                ColumnVector x = LU.Solve(b);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(T * x, b));

                // test inverse
                SquareMatrix TI = LU.Inverse();
                SquareMatrix I = new SquareMatrix(d);
                for (int i = 0; i < d; i++) I[i,i] = 1.0;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(TI * T, I));

            }


        }

    }


}
