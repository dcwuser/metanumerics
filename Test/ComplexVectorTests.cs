using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Matrices;

namespace Test {

    [TestClass]
    public class ComplexVectorTests {

        [TestMethod]
        public void ComplexColumnVectorManipulations () {

            ComplexColumnVector u = new ComplexColumnVector(3);
            u[0] = Complex.I;
            u[2] = Complex.One;
            Assert.IsTrue(u.Dimension == 3);

            ComplexColumnVector v = u.Copy();
            Assert.IsTrue(v == u);
            v[0] += 1.0;
            Assert.IsTrue(v != u);

        }

        [TestMethod]
        public void ComplexColunVectorAsMatrix () {

            ComplexColumnVector u = new ComplexColumnVector(Complex.I, Complex.Zero, Complex.One);

            Assert.IsTrue(u.ColumnCount == 1);
            Assert.IsTrue(u.RowCount == u.Dimension);

            for (int r = 0; r < u.Dimension; r++) {
                Assert.IsTrue(u[r, 0] == u[r]);
            }

            for (int r = 0; r < u.Dimension; r++) {
                u[r, 0] = 0.0;
            }

        }

        [TestMethod]
        public void ComplexColumnVectorArithmetic () {

            ComplexColumnVector u = new ComplexColumnVector(Complex.I, Complex.Zero, Complex.One);
            Assert.IsTrue(u.Dimension == 3);

            Assert.IsTrue(1.0 * u == u);

        }

        [TestMethod]
        public void ComplexVector () {

            ComplexColumnVector v = new ComplexColumnVector(
                Complex.One, Complex.Zero, Complex.I
            );

            Assert.IsTrue(v.Dimension == 3);
            Assert.IsTrue(v[2] == Complex.I);

            ComplexColumnVector v2 = 2.0 * v;
            Assert.IsTrue(v2[1] == 0.0);

            ComplexColumnVector vi = Complex.I * v;
            Assert.IsTrue(vi[2] == -1.0);

            v[1] += 1.0;
            Assert.IsTrue(v[1] == 1.0);

        }

    }
}
