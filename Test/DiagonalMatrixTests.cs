using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Matrices;

namespace Test {

    [TestClass]
    public class DiagonalMatrixTests {

        [TestMethod]
        public void DiagonalMatrixManipulations () {

            DiagonalMatrix D = new DiagonalMatrix(3.0, 2.0, 1.0);
            Assert.IsTrue(D.Dimension == 3);
            Assert.IsTrue(D[0, 0] == 3.0);
            Assert.IsTrue(D[1, 2] == 0.0);
            D[1, 1] -= 1.0;
            Assert.IsTrue(D[1, 1] == 1.0);

        }

    }
}
