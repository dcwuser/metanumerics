using System;

using TestClassAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestClassAttribute;
using TestMethodAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestMethodAttribute;
using ExpectedExceptionAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.ExpectedExceptionAttribute;
using Assert = Microsoft.VisualStudio.TestTools.UnitTesting.Assert;

using Meta.Numerics;
using Meta.Numerics.Statistics;

namespace Test {
    
    [TestClass]
    public class BinaryContingencyTableTest {

        private static ContingencyTable CreateExperiment (double p, double q0, double q1, int N) {

            ContingencyTable e = new ContingencyTable(2, 2);

            Random rng = new Random(1);
            for (int i = 0; i < N; i++) {

                int r, c;
                if (rng.NextDouble() < p) {
                    r = 0;
                    if (rng.NextDouble() < q0) {
                        c = 0;
                    } else {
                        c = 1;
                    }
                } else {
                    r = 1;
                    if (rng.NextDouble() < q1) {
                        c = 0;
                    } else {
                        c = 1;
                    }
                }
                e[r, c] += 1;

            }


            return (e);

        }

        [TestMethod]
        public void BinaryContingencyNullTest () {

            // Create a table with no significant association and test for it.

            ContingencyTable e0 = CreateExperiment(0.25, 0.33, 0.33, 128);

            Assert.IsTrue(e0.Total == 128);
            Assert.IsTrue(e0.RowTotal(0) + e0.RowTotal(1) == e0.Total);
            Assert.IsTrue(e0.ColumnTotal(0) + e0.ColumnTotal(1) == e0.Total);

            UncertainValue lnr = e0.Binary.LogOddsRatio;
            Assert.IsTrue(lnr.ConfidenceInterval(0.95).ClosedContains(0.0));

            UncertainValue r = e0.Binary.OddsRatio;
            Assert.IsTrue(r.ConfidenceInterval(0.95).ClosedContains(1.0));

            TestResult p = e0.PearsonChiSquaredTest();
            Assert.IsTrue(p.Probability > 0.05);

            TestResult f = e0.Binary.FisherExactTest();
            Assert.IsTrue(f.Probability > 0.05);

        }

        [TestMethod]
        [ExpectedException(typeof(InvalidOperationException))]
        public void BinaryContingencyInvalidConstructionTest () {

            int[,] M = new int[2,3];
            ContingencyTable table = new ContingencyTable(2, 3);
            BinaryContingencyTableOperations binary = table.Binary;
        }

        [TestMethod]
        public void BinaryContingencyTest () {

            // Create a table with significant association and test for it.

            ContingencyTable e1 = CreateExperiment(0.50, 0.50, 0.75, 128);

            Assert.IsTrue(e1.RowTotal(0) + e1.RowTotal(1) == e1.Total);
            Assert.IsTrue(e1.ColumnTotal(0) + e1.ColumnTotal(1) == e1.Total);

            UncertainValue lnr = e1.Binary.LogOddsRatio;
            Assert.IsTrue(!lnr.ConfidenceInterval(0.95).ClosedContains(0.0));

            UncertainValue r = e1.Binary.OddsRatio;
            Assert.IsTrue(!r.ConfidenceInterval(0.95).ClosedContains(1.0));

            TestResult p = e1.PearsonChiSquaredTest();
            Assert.IsTrue(p.Probability < 0.05);

            TestResult f = e1.Binary.FisherExactTest();
            Assert.IsTrue(f.Probability < 0.05);

        }

    }
}
