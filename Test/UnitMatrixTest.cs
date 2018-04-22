using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Matrices;

namespace Test {

    [TestClass]
    public class UnitMatrixTest {

        [TestMethod]
        public void UnitMatrixEntries () {

            int n = 3;
            UnitMatrix I = UnitMatrix.OfDimension(n);

            Assert.IsTrue(I.Dimension == n);
            Assert.IsTrue(I.RowCount == n);
            Assert.IsTrue(I.ColumnCount == n);

            for (int r = 0; r < n; r++) {
                for (int c = 0; c < n; c++) {
                    if (r == c) {
                        Assert.IsTrue(I[r, c] == 1.0);
                    } else {
                        Assert.IsTrue(I[r, c] == 0.0);
                    }
                }
            }

        }

        [TestMethod]
        public void UnitMatrixImmutable() {
            UnitMatrix I = UnitMatrix.OfDimension(2);
            try {
                I[0, 0] += 1.0;
                Assert.Fail();
            } catch (InvalidOperationException) { }
        }


        [TestMethod]
        public void UnitMatrixArithmetic () {

            SquareMatrix A = new SquareMatrix(new double[,] {
                { 1.0 , -2.0 },
                { -3.0, 4.0 }
            });
            SquareMatrix AI = A * UnitMatrix.OfDimension(A.Dimension);
            Assert.IsTrue(A == AI);
            SquareMatrix IA = UnitMatrix.OfDimension(A.Dimension) * A;
            Assert.IsTrue(IA == A);

            ColumnVector c = new ColumnVector(1.0, 2.0, 3.0);
            ColumnVector Ic = UnitMatrix.OfDimension(c.Dimension) * c;
            Assert.IsTrue(Ic == c);

            RowVector r = new RowVector(0.0, 1.0);
            RowVector rI = r * UnitMatrix.OfDimension(r.Dimension);
            Assert.IsTrue(rI == r);

            Assert.IsTrue(UnitMatrix.OfDimension(A.Dimension) == UnitMatrix.OfDimension(A.Dimension));

            Assert.IsTrue(0.0 * UnitMatrix.OfDimension(3) == UnitMatrix.OfDimension(3) - UnitMatrix.OfDimension(3));

            Assert.IsTrue(1.0 * UnitMatrix.OfDimension(3) == UnitMatrix.OfDimension(3));

            Assert.IsTrue(2.0 * UnitMatrix.OfDimension(3) == UnitMatrix.OfDimension(3) + UnitMatrix.OfDimension(3));

        }

        [TestMethod]
        public void UnitMatrixNorms () {

            UnitMatrix I = UnitMatrix.OfDimension(4);
            SquareMatrix A = I.ToSquareMatrix();

            Assert.IsTrue(I.OneNorm() == A.OneNorm());
            Assert.IsTrue(I.InfinityNorm() == A.InfinityNorm());
            Assert.IsTrue(I.FrobeniusNorm() == A.FrobeniusNorm());
            Assert.IsTrue(I.MaxNorm() == A.MaxNorm());

        }

        [TestMethod]
        public void UnitMatrixConversions () {

            UnitMatrix I = UnitMatrix.OfDimension(3);

            SquareMatrix A = I.ToSquareMatrix();
            Assert.IsTrue(I == A);
            A[0, 1] += 2.0;
            Assert.IsTrue(I != A);

            SymmetricMatrix B = I.ToSymmetricMatrix();
            Assert.IsTrue(I == B);
            B[0, 1] += 2.0;
            Assert.IsTrue(I != B);

            DiagonalMatrix C = I.ToDiagonalMatrix();
            Assert.IsTrue(I == C);
            C[2, 2] -= 1.0;
            Assert.IsTrue(I != C);

        }
    }
}
