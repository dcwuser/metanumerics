using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace Test {

    [TestClass]
    public class StreamingSampleSummaryTest {

        [TestMethod]
        public void StreamSampleSummaryAgreement () {

            // Streaming properties should give same answers as list methods.

            Random rng = new Random(2);
            List<double> sample = new List<double>(TestUtilities.CreateDataSample(rng, new UniformDistribution(Interval.FromEndpoints(-4.0, 3.0)), 32));

            SummaryStatistics summary = new SummaryStatistics(sample);
            Assert.IsTrue(summary.Count == sample.Count);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(sample.Mean(), summary.Mean));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(sample.Variance(), summary.Variance));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(sample.PopulationMean(), summary.PopulationMean));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(sample.PopulationVariance(), summary.PopulationVariance));

            Assert.IsTrue(sample.Minimum() == summary.Minimum);
            Assert.IsTrue(sample.Maximum() == summary.Maximum);

        }

        [TestMethod]
        public void StreamingSampleSummaryCombination () {

            // Combining partial summaries should give same answer as full summary 

            Random rng = new Random(1);
            List<double> sample = new List<double>(TestUtilities.CreateDataSample(rng, new UniformDistribution(Interval.FromEndpoints(-4.0, 3.0)), 64));

            SummaryStatistics summary = new SummaryStatistics(sample);
            Assert.IsTrue(summary.Count == sample.Count);

            for (int i = 0; i < 4; i++) {

                // Pick a split point in the data
                int m = rng.Next(0, sample.Count);

                // Create a summary of the first part.
                SummaryStatistics summary1 = new SummaryStatistics(sample.Take(m));
                Assert.IsTrue(summary1.Count == m);

                // Create a summary of the second part.
                SummaryStatistics summary2 = new SummaryStatistics(sample.Skip(m));
                Assert.IsTrue(summary2.Count == sample.Count - m);

                // Combine them. Their summary statistics should agree with the original summary.
                SummaryStatistics combined = SummaryStatistics.Combine(summary1, summary2);
                Assert.IsTrue(combined.Count == summary.Count);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(combined.Mean, summary.Mean));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(combined.Variance, summary.Variance));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(combined.StandardDeviation, summary.StandardDeviation));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(combined.Skewness, summary.Skewness));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(combined.PopulationMean, summary.PopulationMean));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(combined.PopulationVariance, summary.PopulationVariance));

                Assert.IsTrue(combined.Minimum == summary.Minimum);
                Assert.IsTrue(combined.Maximum == summary.Maximum);
            }
        }

    }
}
