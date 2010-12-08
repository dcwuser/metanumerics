using Meta.Numerics.Matrices;


using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Test {


    [TestClass()]
    public class VectorTest {


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


        private RowVector R = new RowVector(new double[] { 1, 1 });
        private ColumnVector C = new ColumnVector(new double[] { 1, -1, 1 });
        private RectangularMatrix M = new RectangularMatrix(new double[,] { { 1, 2, 3 }, { 4, 5, 6 } });


        [TestMethod]
        public void ColumnVectorAccessTest () {

            ColumnVector v = new ColumnVector(3);

            // test cloning and equality/inequality testing
            ColumnVector vc = v.Clone();
            Assert.IsTrue(vc == v);
            Assert.IsFalse(vc != v);

            // test independence clone and equality/inequality testing
            vc[0] += 1.0;
            Assert.IsFalse(vc == v);
            Assert.IsTrue(vc != v);


        }

        [TestMethod]
        public void ColumnVectorArithmeticTest () {

            // addition and multiplication by a double
            ColumnVector CA = C + C;
            ColumnVector C2 = 2.0 * C;
            Assert.IsTrue(CA == C2);

            // subtraction and multiplication by a double
            ColumnVector CS = C - C;
            ColumnVector C0 = 0.0 * C;
            Assert.IsTrue(CS == C0);

            // multiply a column by a matrix
            ColumnVector MC = M * C;
            Assert.IsTrue(MC.Dimension == M.RowCount);

            // dot multiply
            RowVector CT = C.Transpose();
            double n = CT * C;
            Assert.IsTrue(n > 0);

        }

        [TestMethod]
        public void RowVectorArithmeticTest () {

            // addition and multiplication by a double
            RowVector RA = R + R;
            RowVector R2 = 2.0 * R;
            Assert.IsTrue(RA == R2);

            // subtraction and multiplication by a double
            RowVector RS = R - R;
            RowVector R0 = 0.0 * R;
            Assert.IsTrue(RS == R0);

            // multiply a column by a matrix
            RowVector RM = R * M;
            Assert.IsTrue(RM.Dimension == M.ColumnCount);

            // dot multiply
            ColumnVector RT = R.Transpose();
            double n = R * RT;
            Assert.IsTrue(n > 0);


        }

        [TestMethod]
        public void MixedVectorArithmeticTest () {

            // outper product
            RectangularMatrix CR = C * R;
            Assert.IsTrue(CR.RowCount == C.Dimension);
            Assert.IsTrue(CR.ColumnCount == R.Dimension);

            // inner product
            double x = R * M * C;

        }
    }
}
