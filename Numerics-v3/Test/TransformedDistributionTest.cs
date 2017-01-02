using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Statistics.Distributions;

namespace Test {

#if FUTURE

    [TestClass]
    public class TransformedDistributionTest {

        [TestMethod]
        public void TestMethod1 () {

            Distribution n0 = new TransformedDistribution(new NormalDistribution(), -2.0, 3.0);
            Distribution n1 = new NormalDistribution(-2.0, 3.0);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(n0.Mean, n1.Mean));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(n0.Variance, n1.Variance));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(n0.StandardDeviation, n1.StandardDeviation));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(n0.Skewness, n1.Skewness));

            for (int k = 0; k < 8; k++) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(n0.Moment(k), n1.Moment(k)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(n0.MomentAboutMean(k), n1.MomentAboutMean(k)));
            }

            foreach (double x in TestUtilities.GenerateUniformRealValues(-8.0, 8.0, 8)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(n0.ProbabilityDensity(x), n1.ProbabilityDensity(x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(n0.LeftProbability(x), n1.LeftProbability(x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(n0.RightProbability(x), n1.RightProbability(x)));
            }

            foreach (double P in TestUtilities.GenerateRealValues(1.0E-4, 1.0, 4)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(n0.InverseLeftProbability(P), n1.InverseLeftProbability(P)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(n0.InverseRightProbability(P), n1.InverseRightProbability(P)));
            }

        }
    }

#endif

}
