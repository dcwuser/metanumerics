using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Matrices;

namespace Test {
    [TestClass]
    public class SparseSquareMatrixTest {

        public SparseSquareMatrix CreateRandomSparseSquareMatrix (int dim, int fill, Random rng) {

            SparseSquareMatrix S = new SparseSquareMatrix(dim);

            for (int i = 0; i < fill; i++) {
                int r = (int)Math.Floor(rng.NextDouble() * dim);
                int c = (int)Math.Floor(rng.NextDouble() * dim);
                S[r, c] = 2.0 * rng.NextDouble() - 1.0;
            }

            return (S);

        }


        [TestMethod]
        public void SparseSquareMatrixManipulation () {

            // 0 0 0 0 0 0
            // 0 1 0 0 4 0
            // 0 5 0 0 0 0
            // 0 0 0 0 0 0
            // 3 2 0 0 6 0
            // 0 0 0 0 0 0

            SparseSquareMatrix S = new SparseSquareMatrix(6);
            Assert.IsTrue(S.Dimension == 6);
            Assert.IsTrue(S.FillCount == 0);
            Assert.IsTrue(S.FillFraction == 0.0);

            // assign values
            S[1, 1] = 1.0;
            S[4, 1] = 2.0;
            S[4, 0] = 3.0;
            S[1, 4] = 4.0;
            S[2, 1] = 5.0;
            S[4, 4] = 6.0;

            Assert.IsTrue(S[1, 1] == 1.0);
            Assert.IsTrue(S[4, 1] == 2.0);
            Assert.IsTrue(S.FillCount == 6);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(S.FillFraction, 6.0 / 36.0));

            // change a value
            S[4, 4] = 7.0;
            Assert.IsTrue(S[4, 4] == 7.0);
            Assert.IsTrue(S.FillCount == 6);

            // remove a value
            S[4, 1] = 0.0;
            Assert.IsTrue(S[4, 1] == 0.0);
            Assert.IsTrue(S.FillCount == 5);

        }


        [TestMethod]
        public void SparseSquareMatrixArithmetic () {

            SparseSquareMatrix S1 = CreateRandomSparseSquareMatrix(8, 16, new Random(2));

            SparseSquareMatrix S2 = 2.0 * S1;

            Assert.IsTrue(S2 != S1);

            SparseSquareMatrix S3 = 0.5 * S2;

            Assert.IsTrue(S3 == S1);

        }

        [TestMethod]
        public void SparseSquareMatrixRowsAndColumns () {

            int d = 6;

            SparseSquareMatrix S = CreateRandomSparseSquareMatrix(d, 2 * d, new Random(1));

            for (int r = 0; r < S.RowCount; r++) {
                RowVector rv = S.Row(r);
                for (int c = 0; c < S.ColumnCount; c++) {
                    Assert.IsTrue(rv[c] == S[r, c]);
                }
            }

            for (int c = 0; c < S.ColumnCount; c++) {
                ColumnVector cv = S.Column(c);
                for (int r = 0; r < S.RowCount; r++) {
                    Assert.IsTrue(cv[r] == S[r, c]);
                }
            }

        }

        [TestMethod]
        public void SparseSquareMatrixCopy () {

            SparseSquareMatrix A = new SparseSquareMatrix(5);
            A[1, 1] = 1.0;
            A[3, 1] = 2.0;
            A[1, 3] = 3.0;
            A[3, 3] = 4.0;
            A[2, 4] = 5.0;

            // make a copy
            SparseSquareMatrix AC = A.Copy();

            // test that the copy agrees

            Assert.IsTrue(AC.Dimension == A.Dimension);
            Assert.IsTrue(AC.FillCount == A.FillCount);

            Assert.IsTrue(AC[2, 2] == A[2, 2]);
            Assert.IsTrue(AC[3, 3] == A[3, 3]);
            
            // test the copy's independence

            A[3, 3] += 1.0;
            Assert.IsTrue(AC[3,3] != A[3,3]);

            A[1, 3] = 0.0;
            Assert.IsTrue(AC[1, 3] != A[1, 3]);
            Assert.IsTrue(AC.FillCount != A.FillCount);



        }

        [TestMethod]
        public void SparseSquareMatrixAgreement () {

            int d = 6;
            SparseSquareMatrix A = new SparseSquareMatrix(d);
            SquareMatrix B = new SquareMatrix(d);

            Random rng = new Random(1);
            for (int i = 0; i < 2 * d; i++) {
                int r = (int) Math.Floor(rng.NextDouble() * d);
                int c = (int) Math.Floor(rng.NextDouble() * d);
                A[r, c] = 2.0 * rng.NextDouble() - 1.0;
                B[r, c] = A[r, c];
            }

            RowVector u = new RowVector(d);
            ColumnVector v = new ColumnVector(d);
            for (int i = 0; i < d; i++) {
                u[i] = 2.0 * rng.NextDouble() - 1.0;
                v[i] = 2.0 * rng.NextDouble() - 1.0;
            }

            RowVector uA = u * A;
            RowVector uB = u * B;
            Assert.IsTrue(TestUtilities.IsNearlyEqual(uA, uB));

            ColumnVector Av = A * v;
            ColumnVector Bv = B * v;
            Assert.IsTrue(TestUtilities.IsNearlyEqual(Av, Bv));

        }

        [TestMethod]
        public void SparseSquareMatrixSolutionAgreement () {

            Random rng = new Random(2);

            int n = 16;
            //for (int n = 8; n < 100; n+= 11) {

                // create dense and sparse matrices with the same entries
                SquareMatrix M = new SquareMatrix(n);
                SparseSquareMatrix S = new SparseSquareMatrix(n);
                for (int r = 0; r < n; r++) {
                    for (int c = 0; c < n; c++) {
                        if (rng.NextDouble() < 0.5) {
                            M[r, c] = rng.NextDouble();
                            S[r, c] = M[r, c];
                        }
                    }
                }

                // pick a RHS
                ColumnVector b = new ColumnVector(n);
                for (int i = 0; i < n; i++) b[i] = rng.NextDouble();

                // solve each
                ColumnVector Mx = M.LUDecomposition().Solve(b);
                ColumnVector Sx = S.Solve(b);

                // the solutions should be the same
                Assert.IsTrue(TestUtilities.IsNearlyEqual(Mx, Sx, TestUtilities.TargetPrecision * 100.0));

            //}

        }

        [TestMethod]
        public void SparseSquareMatrixPotential () {

            // An 2D electrostatic boundary value problem in cartesian coordinates
            // A square of length n, with specified constant potentials on each wall
            // Discritized Laplace equation is
            //  u_{x,y-1} + u_{x+1,y} + u_{x,y+1} + u_{x-1,y} - 4 u_{x,y} = 0
            // Number interior points sequentially row-wise i.e. i = x + n * y
            // Points nearest the walls will pick up a boundary value, which because it is not a variable we move to the RHS
            // Points closer to the interior pick up no boundary value, so their RHS is zero

            int n = 100;

            double pn = 0.0; double pe = 1.0; double ps = 0.0; double pw = 1.0;

            SparseSquareMatrix A = new SparseSquareMatrix(n * n);
            ColumnVector b = new ColumnVector(n * n);

            // set up A and b
            for (int y = 0; y < n; y++) {
                for (int x = 0; x < n; x++) {
                    int i = x + n * y;
                    // center value
                    A[i, i] = 4.0;
                    // north
                    if (y == 0) {
                        b[i] += pn;
                    } else {
                        int j = x + n * (y - 1);
                        A[i, j] = -1.0;
                    }
                    // east
                    if (x == (n-1)) {
                        b[i] += pe;
                    } else {
                        int j = (x + 1) + n * y;
                        A[i, j] = -1.0;
                    }
                    // south
                    if (y == (n-1)) {
                        b[i] += ps;
                    } else {
                        int j = x + n * (y + 1);
                        A[i, j] = -1.0;
                    }
                    // west
                    if (x == 0) {
                        b[i] += pw;
                    } else {
                        int j = (x - 1) + n * y;
                        A[i, j] = -1.0;
                    }
                }
            }

            ColumnVector u = A.Solve(b);

            for (int y = 0; y < 10; y++) {
                for (int x = 0; x < 10; x++) {
                    int i = x + n * y;
                    Console.Write("{0} ", u[i]);
                }
                Console.WriteLine();
            }

            Console.WriteLine(PotentialSolution(1, 2, n));

        }


        private static double PotentialSolution (int x, int y, int n) {

            double f = 0.0;
            for (int m = 1; m < 128; m += 2) {
                double f_old = f;
                double p = m * Math.PI / (n + 1);
                double q = m * Math.PI / 2.0;
                f += Math.Cosh(p * (x + 1) - q) * Math.Sin(p * (y + 1)) / Math.Cosh(q) / m;
                Console.WriteLine(f);
                //if (f == f_old) return (4 / Math.PI * f);
            }
            return (4 / Math.PI * f);
        }


    }
}
