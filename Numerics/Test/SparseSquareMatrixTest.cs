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

    }
}
