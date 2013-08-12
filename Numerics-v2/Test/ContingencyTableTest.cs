using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Statistics;

namespace Test {

    [TestClass]
    public class ContingencyTableTest {

        [TestMethod]
        public void ContingencyTableOperations () {

            ContingencyTable t = new ContingencyTable(4, 3);
            Assert.IsTrue(t.RowCount == 4);
            Assert.IsTrue(t.ColumnCount == 3);

            Assert.IsTrue(t.RowTotal(2) == 0);
            Assert.IsTrue(t.ColumnTotal(1) == 0);
            Assert.IsTrue(t.Total == 0);

            t[1, 1] = 2;
            Assert.IsTrue(t.RowTotal(2) == 0);
            Assert.IsTrue(t.ColumnTotal(1) == 2);
            Assert.IsTrue(t.Total == 2);

        }

        private void ChooseRandomCell (double[,] pp, double p, out int r, out int c) {
            double ps = 0.0;
            r = 0; c = 0;
            while (r < pp.GetLength(0)) {
                c = 0;
                while (c < pp.GetLength(1)) {
                    ps += pp[r, c];
                    if (ps >= p) return;
                    c++;
                }
                r++;
            }
        }

        [TestMethod]
        public void ContingencyTableProbabilitiesAndUncertainties () {

            // start with an underlying population
            double[,] pp = new double[,]
                { { 1.0 / 45.0, 2.0 / 45.0, 3.0 / 45.0 },
                  { 4.0 / 45.0, 5.0 / 45.0, 6.0 / 45.0 },
                  { 7.0 / 45.0, 8.0 / 45.0, 9.0 / 45.0 } };

            // form 50 contingency tables, each with N = 50
            Random rng = new Random(314159);
            BivariateSample p22s = new BivariateSample();
            BivariateSample pr0s = new BivariateSample();
            BivariateSample pc1s = new BivariateSample();
            BivariateSample pr2c0s = new BivariateSample();
            BivariateSample pc1r2s = new BivariateSample();
            for (int i = 0; i < 50; i++) {

                ContingencyTable T = new ContingencyTable(3, 3);
                for (int j = 0; j < 50; j++) {
                    int r, c;
                    ChooseRandomCell(pp, rng.NextDouble(), out r, out c);
                    T.Increment(r, c);
                }

                Assert.IsTrue(T.Total == 50);

                // for each contingency table, compute estimates of various population quantities

                UncertainValue p22 = T.Probability(2, 2);
                UncertainValue pr0 = T.ProbabilityOfRow(0);
                UncertainValue pc1 = T.ProbabilityOfColumn(1);
                UncertainValue pr2c0 = T.ProbabilityOfRowConditionalOnColumn(2, 0);
                UncertainValue pc1r2 = T.ProbabilityOfColumnConditionalOnRow(1, 2);
                p22s.Add(p22.Value, p22.Uncertainty);
                pr0s.Add(pr0.Value, pr0.Uncertainty);
                pc1s.Add(pc1.Value, pc1.Uncertainty);
                pr2c0s.Add(pr2c0.Value, pr2c0.Uncertainty);
                pc1r2s.Add(pc1r2.Value, pc1r2.Uncertainty);

            }

            // the estimated population mean of each probability should include the correct probability in the underlyting distribution
            Assert.IsTrue(p22s.X.PopulationMean.ConfidenceInterval(0.95).ClosedContains(9.0 / 45.0));
            Assert.IsTrue(pr0s.X.PopulationMean.ConfidenceInterval(0.95).ClosedContains(6.0 / 45.0));
            Assert.IsTrue(pc1s.X.PopulationMean.ConfidenceInterval(0.95).ClosedContains(15.0 / 45.0));
            Assert.IsTrue(pr2c0s.X.PopulationMean.ConfidenceInterval(0.95).ClosedContains(7.0 / 12.0));
            Assert.IsTrue(pc1r2s.X.PopulationMean.ConfidenceInterval(0.95).ClosedContains(8.0 / 24.0));

            // the estimated uncertainty for each population parameter should be the standard deviation across independent measurements
            // since the reported uncertainly changes each time, we use the mean value for comparison
            Assert.IsTrue(p22s.X.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(p22s.Y.Mean));
            Assert.IsTrue(pr0s.X.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(pr0s.Y.Mean));
            Assert.IsTrue(pc1s.X.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(pc1s.Y.Mean));
            Assert.IsTrue(pr2c0s.X.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(pr2c0s.Y.Mean));
            Assert.IsTrue(pc1r2s.X.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(pc1r2s.Y.Mean));

        }

#if FUTURE

        [TestMethod]
        public void ContingencyTableNamedOperations () {

            ContingencyTable t = new ContingencyTable(2, 3);
            t.RowNames[0] = "Male";
            t.RowNames[1] = "Female";
            t.ColumnNames[0] = "Party 1";
            t.ColumnNames[1] = "Party 2";
            t.ColumnNames[2] = "Party 3";
            t["Male", "Party 1"] = 10;
            t["Male", "Party 2"] = 20;
            t["Male", "Party 3"] = 30;
            t["Female", "Party 1"] = 30;
            t["Female", "Party 2"] = 20;
            t["Female", "Party 3"] = 10;

        }

#endif

    }
}
