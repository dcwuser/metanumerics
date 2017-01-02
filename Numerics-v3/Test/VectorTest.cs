using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Matrices;


namespace Test {

    [TestClass]
    public class VectorTest {

        private RowVector R = new RowVector(new double[] { 1, 1 });
        private ColumnVector C = new ColumnVector(new double[] { 1, -1, 1 });
        private RectangularMatrix M = new RectangularMatrix(new double[,] { { 1, 2, 3 }, { 4, 5, 6 } });


        [TestMethod]
        public void ColumnVectorCopy () {

            ColumnVector v = new ColumnVector(3);

            // test cloning and equality/inequality testing
            ColumnVector vc = v.Copy();
            Assert.IsTrue(vc == v);
            Assert.IsFalse(vc != v);

            // test independence clone and equality/inequality testing
            vc[0] += 1.0;
            Assert.IsFalse(vc == v);
            Assert.IsTrue(vc != v);


        }

        [TestMethod]
        public void RowVectorCopy () {

            RowVector v = new RowVector(3);

            // test cloning and equality/inequality testing
            RowVector vc = v.Copy();
            Assert.IsTrue(vc == v);
            Assert.IsFalse(vc != v);

            // test independence clone and equality/inequality testing
            vc[0] += 1.0;
            Assert.IsFalse(vc == v);
            Assert.IsTrue(vc != v);

        }

        [TestMethod]
        public void ColumnVectorArithmetic () {

            // addition and multiplication by a double
            ColumnVector CA = C + C;
            ColumnVector C2 = 2.0 * C;
            Assert.IsTrue(CA == C2);

            // subtraction and multiplication by a double
            ColumnVector CS = C - C;
            ColumnVector C0 = 0.0 * C;
            Assert.IsTrue(CS == C0);

            // negation and division by -1
            ColumnVector CN = -C;
            ColumnVector CD = C / (-1.0);
            Assert.IsTrue(CN == CD);

            // multiply a column by a matrix
            ColumnVector MC = M * C;
            Assert.IsTrue(MC.Dimension == M.RowCount);

            // dot multiply
            RowVector CT = C.Transpose();
            double dot = CT * C;
            Assert.IsTrue(dot > 0);

        }

        [TestMethod]
        public void RowVectorArithmetic () {

            // addition and multiplication by a double
            RowVector RA = R + R;
            RowVector R2 = 2.0 * R;
            Assert.IsTrue(RA == R2);

            // subtraction and multiplication by a double
            RowVector RS = R - R;
            RowVector R0 = 0.0 * R;
            Assert.IsTrue(RS == R0);

            // negation and multiplication by -1
            RowVector RN = -R;
            RowVector RD = R / (-1.0);
            Assert.IsTrue(RN == RD);

            // multiply a column by a matrix
            RowVector RM = R * M;
            Assert.IsTrue(RM.Dimension == M.ColumnCount);

            // dot multiply
            ColumnVector RT = R.Transpose();
            double dot = R * RT;
            Assert.IsTrue(dot > 0);

        }

        [TestMethod]
        public void MixedVectorArithmetic () {

            // outer product
            RectangularMatrix CR = C * R;
            Assert.IsTrue(CR.RowCount == C.Dimension);
            Assert.IsTrue(CR.ColumnCount == R.Dimension);

            // inner product
            double x = R * M * C;

        }

        [TestMethod]
        public void VectorNorm () {

            RowVector nR = R / R.Norm();
            Assert.IsTrue(TestUtilities.IsNearlyEqual(nR.Norm(), 1.0));

            ColumnVector nC = C / C.Norm();
            Assert.IsTrue(TestUtilities.IsNearlyEqual(nC.Norm(), 1.0));

        }

        [TestMethod]
        public void VectorAsCollection () {

            ColumnVector v = new ColumnVector(new double[] { 1.0, 2.0, 3.0 });

            IList<double> vl = v as IList<double>;
            for (int i = 0; i < vl.Count; i++) {
                Assert.IsTrue(vl[i] == v[i]);
            }

            ICollection<double> vc = v as ICollection<double>;
            Assert.IsTrue(vc.Count == v.Dimension);
            Assert.IsTrue(vc.Contains(1.0));

            // IEnumerable
            int index = 0;
            foreach (double vi in v) {
                Assert.IsTrue(vi == v[index]);
                index++;
            }


        }

        [TestMethod]
        public void ExtendedVectorNorm () {

            // Check that vector norm is correctly computed even when vector is enormous (square would overflow) or tiny (square woudl underflow)

            double x = Double.MaxValue / 32.0;

            // Use the pythagorian quadrouple (2, 3, 6, 7)
            ColumnVector v = new ColumnVector(2.0, 3.0, 6.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(v.Norm(), 7.0));
            ColumnVector bv = x * v;
            Assert.IsTrue(TestUtilities.IsNearlyEqual(bv.Norm(), 7.0 * x));
            ColumnVector sv = v / x;
            Assert.IsTrue(TestUtilities.IsNearlyEqual(sv.Norm(), 7.0 / x));

            // Use the pythagorian 7-tuple (1, 2, 3, 4, 5, 27, 28)
            RowVector u = new RowVector(1.0, 2.0, 3.0, 4.0, 5.0, 27.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(u.Norm(), 28.0));
            RowVector bu = x * u;
            Assert.IsTrue(TestUtilities.IsNearlyEqual(bu.Norm(), 28.0 * x));
            RowVector su = u / x;
            Assert.IsTrue(TestUtilities.IsNearlyEqual(su.Norm(), 28.0 / x));

        }

    }
}
