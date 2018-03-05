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

            SampleSummary summary = new SampleSummary(sample);
            Assert.IsTrue(summary.Count == sample.Count);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(sample.Mean(), summary.Mean));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(sample.RawMoment(2), summary.SecondRawMoment));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(sample.RawMoment(3), summary.ThirdRawMoment));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(sample.Variance(), summary.Variance));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(sample.CentralMoment(3), summary.ThirdCentralMoment));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(sample.CentralMoment(4), summary.FourthCentralMoment));
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

            SampleSummary summary = new SampleSummary(sample);
            Assert.IsTrue(summary.Count == sample.Count);

            for (int i = 0; i < 4; i++) {

                // Pick a split point in the data
                int m = rng.Next(0, sample.Count);

                // Create a summary of the first part.
                SampleSummary summary1 = new SampleSummary(sample.Take(m));
                Assert.IsTrue(summary1.Count == m);

                // Create a summary of the second part.
                SampleSummary summary2 = new SampleSummary(sample.Skip(m));
                Assert.IsTrue(summary2.Count == sample.Count - m);

                // Combine them. Their summary statistics should agree with the original summary.
                SampleSummary combined = SampleSummary.Combine(summary1, summary2);
                Assert.IsTrue(combined.Count == summary.Count);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(combined.Mean, summary.Mean));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(combined.Variance, summary.Variance));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(combined.StandardDeviation, summary.StandardDeviation));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(combined.SecondRawMoment, summary.SecondRawMoment));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(combined.ThirdCentralMoment, summary.ThirdCentralMoment));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(combined.FourthCentralMoment, summary.FourthCentralMoment));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(combined.PopulationMean, summary.PopulationMean));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(combined.PopulationVariance, summary.PopulationVariance));

                Assert.IsTrue(combined.Minimum == summary.Minimum);
                Assert.IsTrue(combined.Maximum == summary.Maximum);
            }
        }

    }
}
