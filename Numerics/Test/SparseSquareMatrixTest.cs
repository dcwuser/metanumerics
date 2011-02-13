using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Matrices;

namespace Test {
    [TestClass]
    public class SparseSquareMatrixTest {

        [TestMethod]
        public void TestMethod1 () {

            // 0 0 0 0 0 0
            // 0 1 0 0 4 0
            // 0 5 0 0 0 0
            // 0 0 0 0 0 0
            // 3 2 0 0 6 0
            // 0 0 0 0 0 0

            SparseSquareMatrix S = new SparseSquareMatrix(6);
            Assert.IsTrue(S.FillCount == 0);

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

            // change a value
            S[4, 4] = 7.0;
            Assert.IsTrue(S[4, 4] == 7.0);
            Assert.IsTrue(S.FillCount == 6);

            // remove a value
            S[4, 1] = 0.0;
            Assert.IsTrue(S[4, 1] == 0.0);
            Assert.IsTrue(S.FillCount == 5);

        }

    }
}
