using Meta.Numerics.Matrices;
using Microsoft.VisualStudio.TestTools.UnitTesting;
namespace Test
{


    [TestClass()]
    public class MixedMatrixTest {


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

        private static Matrix R;

        private static SquareMatrix M;

        private static SymmetricMatrix S;


        static MixedMatrixTest () {

            R = new Matrix(2, 3);
            R[0, 0] = 0;
            R[0, 1] = 1;
            R[0, 2] = 2;
            R[1, 0] = 3;
            R[1, 1] = 4;
            R[1, 2] = 5;

            M = new SquareMatrix(3);
            M[0, 0] = 0;
            M[0, 1] = 1;
            M[0, 2] = 2;
            M[1, 0] = 3;
            M[1, 1] = 4;
            M[1, 2] = 5;
            M[2, 0] = 6;
            M[2, 1] = 7;
            M[2, 2] = 8;

            S = new SymmetricMatrix(3);
            M[0, 0] = 0;
            M[0, 1] = 1;
            M[0, 2] = 2;
            M[1, 1] = 3;
            M[1, 2] = 4;
            M[2, 2] = 5;

        }

        [TestMethod]
        public void MixedMatrixAddition () {

            SquareMatrix MS = M + S;
            SquareMatrix SM = S + M;

        }

        [TestMethod]
        public void MixedMatrixSubtraction () {

            SquareMatrix MS = M - S;
            SquareMatrix SM = S - M;

        }

        [TestMethod]
        public void MixedMatrixMultiplication () {

            Matrix RM = R * M;
            Matrix RS = R * S;
            SquareMatrix MS = M * S;
            SquareMatrix SM = S * M;

        }

    }

}
