using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Statistics;

namespace Test {

    [TestClass]
    public class HistogramTest {

        private static void AssertBinCounts (Histogram histogram, int[] counts) {

            Assert.IsTrue(histogram.BelowRangeBin.Counts == counts[0]);
            for (int i = 0; i < histogram.Bins.Count; i++) {
                Assert.IsTrue(histogram.Bins[i].Counts == counts[i + 1]);
            }
            Assert.IsTrue(histogram.AboveRangeBin.Counts == counts[histogram.Bins.Count + 1]);

        }

        [TestMethod]
        public void HistogramOperations () {

            Histogram h = new Histogram(new double[] { 1.0, 3.0, 5.0, 7.0 });

            Assert.IsTrue(h.Bins.Count == 3);
            Assert.IsTrue(h.Range == Interval.FromEndpoints(1.0, 7.0));

            Assert.IsTrue(h.BelowRangeBin.Range == Interval.FromEndpoints(Double.NegativeInfinity, 1.0));
            Assert.IsTrue(h.Bins[0].Range == Interval.FromEndpoints(1.0, 3.0));
            Assert.IsTrue(h.Bins[1].Range == Interval.FromEndpoints(3.0, 5.0));
            Assert.IsTrue(h.Bins[2].Range == Interval.FromEndpoints(5.0, 7.0));
            Assert.IsTrue(h.AboveRangeBin.Range == Interval.FromEndpoints(7.0, Double.PositiveInfinity));

            AssertBinCounts(h, new int[] { 0, 0, 0, 0, 0 });

            h.Add(-1.0);
            AssertBinCounts(h, new int[] { 1, 0, 0, 0, 0 });

            h.Add(0.0);
            AssertBinCounts(h, new int[] { 2, 0, 0, 0, 0 });

            h.Add(1.0);
            AssertBinCounts(h, new int[] { 2, 1, 0, 0, 0 });

            h.Add(2.0);
            AssertBinCounts(h, new int[] { 2, 2, 0, 0, 0 });

            h.Add(3.0);
            AssertBinCounts(h, new int[] { 2, 2, 1, 0, 0 });

            h.Add(4.0);
            AssertBinCounts(h, new int[] { 2, 2, 2, 0, 0 });

            h.Add(5.0);
            AssertBinCounts(h, new int[] { 2, 2, 2, 1, 0 });

            h.Add(6.0);
            AssertBinCounts(h, new int[] { 2, 2, 2, 2, 0 });

            h.Add(7.0);
            AssertBinCounts(h, new int[] { 2, 2, 2, 2, 1 });

            h.Add(8.0);
            AssertBinCounts(h, new int[] { 2, 2, 2, 2, 2 });

            Assert.IsTrue(h.TotalCounts == 10);
            //Assert.IsTrue(h.OutOfRangeCounts == 4);

            h.Clear();
            AssertBinCounts(h, new int[] { 0, 0, 0, 0, 0 });
            Assert.IsTrue(h.TotalCounts == 0);

            foreach (HistogramBin bin in h.Bins) {
                bin.Increment();
            }
            AssertBinCounts(h, new int[] { 0, 1, 1, 1, 0 });

            foreach (HistogramBin bin in h.Bins) {
                bin.Deincrement();
            }
            AssertBinCounts(h, new int[] { 0, 0, 0, 0, 0 });

            //Assert.IsTrue(h.OutOfRangeCounts == 0);

        }

    }
}
