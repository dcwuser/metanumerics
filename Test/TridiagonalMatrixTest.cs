using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;


using Meta.Numerics.Matrices;
namespace Test
{

#if FUTURE

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

            Assert.IsTrue(T.Dimension == 4);

            TridiagonalMatrix TC = T.Clone();

            Assert.IsTrue(T[0, 0] == TC[0, 0]);

        }

        [TestMethod]
        public void TridiagonalMatrixLUDecompositionTest () {

            int d = 4;

            TridiagonalMatrix T = CreateRandomTridiagonalMatrix(d);
            TridiagonalLUDecomposition LU = T.LUDecompose();

            Assert.IsTrue(LU.Dimension == T.Dimension);

            //Assert.IsTrue(TestUtilities.IsNearlyEqual(LU.Determinant(), T.Determinant()));

            SquareMatrix P = LU.PMatrix();
            SquareMatrix L = LU.LMatrix();
            SquareMatrix U = LU.UMatrix();

            SquareMatrixTest.PrintMatrix(T);
            SquareMatrixTest.PrintMatrix(L);
            SquareMatrixTest.PrintMatrix(U);
            SquareMatrixTest.PrintMatrix(L * U);


            //ColumnVector b = new ColumnVector(TestUtilities.GenerateRealValues(0.1, 10.0, d));
            ColumnVector b = new ColumnVector( new double[] { 1.0, 2.0, 3.0, 4.0 } );
            ColumnVector x = LU.Solve(b);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(T * x, b));

            ColumnVector e0 = new ColumnVector( new double[] { 1.0, 0.0, 0.0, 0.0 } );
            ColumnVector v0 = LU.Solve(e0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(T * v0, e0));

            //LU.Inverse();


        }

    }

#endif

}
