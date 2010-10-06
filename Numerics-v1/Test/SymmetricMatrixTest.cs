using System;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Matrices;
using Meta.Numerics.Functions;

namespace Test {

    
    [TestClass()]
    public class SymmetricMatrixTest {


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


        private static SymmetricMatrix CreateSymmetricRandomMatrix (int n, int seed) {
            SymmetricMatrix M = new SymmetricMatrix(n);
            Random rng = new Random(seed);
            for (int r = 0; r < n; r++) {
                for (int c = 0; c <= r; c++) {
                    M[r, c] = 2.0 * rng.NextDouble() - 1.0;
                }
            }
            return (M);
        }

        [TestMethod]
        public void SymmetricMatrixAccess () {

            SymmetricMatrix M = CreateSymmetricRandomMatrix(4, 1);

            /*
            // check column
            ColumnVector c = M.Column(1);
            Assert.IsTrue(c.Dimension == M.Dimension);
            for (int i = 0; i < c.Dimension; i++) {
                Assert.IsTrue(c[i] == M[i, 1]);
            }

            // check row
            RowVector r = M.Row(1);
            Assert.IsTrue(r.Dimension == M.Dimension);
            for (int i = 0; i < r.Dimension; i++) {
                Assert.IsTrue(r[i] == c[i]);
            }
            */

            // check clone
            SymmetricMatrix MC = M.Clone();
            Assert.IsTrue(MC == M);
            Assert.IsFalse(MC != M);

            // check independence of clone
            MC[0, 1] += 1.0;
            Assert.IsFalse(MC == M);
            Assert.IsTrue(MC != M);

            // check that update was symmetric
            Assert.IsTrue(M[0, 1] == M[1, 0]);

        }

        [TestMethod]
        public void SymmetricMatrixArithmeticTest () {
            SymmetricMatrix M = CreateSymmetricRandomMatrix(5,1);

            // addition is same a multiplication by two
            SymmetricMatrix MA = M + M;
            SymmetricMatrix M2 = 2.0 * M;
            Assert.IsTrue(MA == M2);

            // subraction of self same as multiplication by zero
            SymmetricMatrix MS = M - M;
            SymmetricMatrix M0 = 0.0 * M;
            Assert.IsTrue(MS == M0);

            // matrix multiplication
            SquareMatrix MM = M * M;

        }

        [TestMethod]
        public void SymmetricHilbertMatrixInverse () {
            for (int d = 1; d <= 4; d++) {
                Console.WriteLine("d={0}", d);
                SquareMatrix I = TestUtilities.CreateSquareUnitMatrix(d);
                SymmetricMatrix H = TestUtilities.CreateSymmetricHilbertMatrix(d);
                SymmetricMatrix HI = H.Inverse();
                Assert.IsTrue(TestUtilities.IsNearlyEqual(H * HI, I));
            }
            // fails for d > 4! look into this
        }

        [TestMethod]
        public void SymmetricRandomMatrixInverse () {
            for (int d = 1; d <= 100; d = d + 11) {
                Console.WriteLine("d={0}", d);
                SquareMatrix I = TestUtilities.CreateSquareUnitMatrix(d);
                SymmetricMatrix M = TestUtilities.CreateSymmetricRandomMatrix(d, 1);
                SymmetricMatrix MI = M.Inverse();
                Assert.IsTrue(TestUtilities.IsNearlyEqual(M * MI, I));
            }

        }

        [TestMethod]
        public void SymmetricRandomMatrixEigenvalue () {
            for (int d = 1; d <= 100; d = d + 11) {
                Console.WriteLine("d={0}", d);
                SymmetricMatrix M = CreateSymmetricRandomMatrix(d, 1);
                double tr = M.Trace();
                DateTime start = DateTime.Now;
                RealEigensystem E = M.Eigensystem();
                DateTime finish = DateTime.Now;
                Console.WriteLine("  {0} ms", (finish - start).Milliseconds);
                Assert.IsTrue(E.Dimension == M.Dimension);
                double[] es = new double[d];
                for (int i = 0; i < d; i++) {
                    double e = E.Eigenvalue(i);
                    ColumnVector v = E.Eigenvector(i);
                    Assert.IsTrue(TestUtilities.IsNearlyEigenpair(M, v, e));
                    es[i] = e;
                }
                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(es, tr));
            }
        }

        [TestMethod]
        public void SymmetricHilbertMatrixEigenvalues () {
            for (int d = 1; d <= 8; d++) {
                Console.WriteLine("d={0}", d);
                SymmetricMatrix H = TestUtilities.CreateSymmetricHilbertMatrix(d);
                double tr = H.Trace();
                RealEigensystem E = H.Eigensystem();
                double[] es = new double[d];
                for (int i = 0; i < d; i++) {
                    double e = E.Eigenvalue(i);
                    ColumnVector v = E.Eigenvector(i);
                    Assert.IsTrue(TestUtilities.IsNearlyEigenpair(H, v, e));
                    es[i] = e;
                }
                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(es, tr));
            }
        }


        [TestMethod]
        public void SymmetricMatrixDecomposition () {
            for (int d = 1; d <= 4; d++) {
                SymmetricMatrix H = TestUtilities.CreateSymmetricHilbertMatrix(d);

                CholeskyDecomposition CD = H.CholeskyDecomposition();
                Assert.IsTrue(CD != null, String.Format("d={0} not positive definite", d));
                Assert.IsTrue(CD.Dimension == d);
                SymmetricMatrix HI = CD.Inverse();
                IMatrix I = TestUtilities.CreateSquareUnitMatrix(d);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(H * HI, I));
            }
        }

        [TestMethod]
        public void SymmetricRandomMatrixCholeskyDecomposition () {

            int d = 100;
            Random rng = new Random(d);
            ColumnVector[] V = new ColumnVector[d];
            for (int i=0; i < d; i++) {
                V[i] = new ColumnVector(d);
                for (int j = 0; j < d; j++) {
                    V[i][j] = rng.NextDouble();
                }
            }

            SymmetricMatrix A = new SymmetricMatrix(d);
            for (int i = 0; i < d; i++) {
                for (int j = 0; j <= i; j++) {
                    A[i, j] = V[i].Transpose() * V[j];
                }
            }

            Stopwatch s = Stopwatch.StartNew();
            CholeskyDecomposition CD = A.CholeskyDecomposition();
            s.Stop();
            Console.WriteLine("{0} {1}", d, s.ElapsedMilliseconds);

            Assert.IsTrue(CD != null);

        }

        [TestMethod]
        public void CatalanHankelMatrixDeterminant () {

            for (int d = 1; d <= 8; d++) {

                SymmetricMatrix S = new SymmetricMatrix(d);
                for (int r = 0; r < d; r++) {
                    for (int c = 0; c <= r; c++) {
                        int n = r + c;
                        S[r, c] = AdvancedIntegerMath.BinomialCoefficient(2*n, n) / (n + 1);
                    }
                }

                CholeskyDecomposition CD = S.CholeskyDecomposition();
                Assert.IsTrue(TestUtilities.IsNearlyEqual(CD.Determinant(), 1.0));

            }

        }

    }
}









