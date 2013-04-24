using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace Test {

    [TestClass]
    public class NullDistributionTests {

        [TestMethod]
        public void KolmogorovNullDistributionTest () {

            // The distribution is irrelevent; pick one at random 
            Distribution sampleDistribution = new LognormalDistribution();

            // Loop over various sample sizes
            foreach (int n in TestUtilities.GenerateIntegerValues(2, 128, 16)) {

                // Create a sample to hold the KS statistics
                Sample testStatistics = new Sample();
                // and a variable to hold the claimed null distribution, which should be the same for each test
                Distribution nullDistribution = null;

                // Create a bunch of samples, each with n data points
                for (int i = 0; i < 256; i++) {

                    // Just use n+i as a seed in order to get different points each time
                    Sample sample = TestUtilities.CreateSample(sampleDistribution, n, 512 * n + i + 1);

                    // Do a KS test of the sample against the distribution each time
                    TestResult r1 = sample.KolmogorovSmirnovTest(sampleDistribution);

                    // Record the test statistic value and the claimed null distribution
                    testStatistics.Add(r1.Statistic);
                    nullDistribution = r1.Distribution;

                }

                // Do a Kuiper test of our sample of KS statistics against the claimed null distribution
                // We could use a KS test here instead, which would be way cool and meta, but we picked Kuiper instead for variety
                TestResult r2 = testStatistics.KuiperTest(nullDistribution);
                Console.WriteLine("{0} {1} {2}", n, r2.Statistic, r2.LeftProbability);
                Assert.IsTrue(r2.RightProbability > 0.05);

                // Test moment matches, too
                Console.WriteLine(" {0} {1}", testStatistics.PopulationMean, nullDistribution.Mean);
                Console.WriteLine(" {0} {1}", testStatistics.PopulationVariance, nullDistribution.Variance);
                Assert.IsTrue(testStatistics.PopulationMean.ConfidenceInterval(0.99).ClosedContains(nullDistribution.Mean));
                Assert.IsTrue(testStatistics.PopulationVariance.ConfidenceInterval(0.99).ClosedContains(nullDistribution.Variance));

            }

        }

        [TestMethod]
        public void KuiperNullDistributionTest () {

            // The distribution is irrelevent; pick one at random 
            Distribution sampleDistribution = new NormalDistribution();

            // Loop over various sample sizes
            foreach (int n in TestUtilities.GenerateIntegerValues(2, 128, 16)) {

                // Create a sample to hold the KS statistics
                Sample testStatistics = new Sample();
                // and a variable to hold the claimed null distribution, which should be the same for each test
                Distribution nullDistribution = null;

                // Create a bunch of samples, each with n+1 data points
                // We pick n+1 instead of n just to have different sample size values than in the KS test case
                for (int i = 0; i < 256; i++) {

                    // Just use n+i as a seed in order to get different points each time
                    Sample sample = TestUtilities.CreateSample(sampleDistribution, n + 1, 512 * n + i + 2);

                    // Do a Kuiper test of the sample against the distribution each time
                    TestResult r1 = sample.KuiperTest(sampleDistribution);

                    // Record the test statistic value and the claimed null distribution
                    testStatistics.Add(r1.Statistic);
                    nullDistribution = r1.Distribution;

                }

                // Do a KS test of our sample of Kuiper statistics against the claimed null distribution
                // We could use a Kuiper test here instead, which would be way cool and meta, but we picked KS instead for variety
                TestResult r2 = testStatistics.KolmogorovSmirnovTest(nullDistribution);
                Console.WriteLine("{0} {1} {2}", n, r2.Statistic, r2.LeftProbability);
                Assert.IsTrue(r2.RightProbability > 0.01);

                // Test moment matches, too
                Console.WriteLine(" {0} {1}", testStatistics.PopulationMean, nullDistribution.Mean);
                Console.WriteLine(" {0} {1}", testStatistics.PopulationVariance, nullDistribution.Variance);
                Assert.IsTrue(testStatistics.PopulationMean.ConfidenceInterval(0.99).ClosedContains(nullDistribution.Mean));
                Assert.IsTrue(testStatistics.PopulationVariance.ConfidenceInterval(0.99).ClosedContains(nullDistribution.Variance));

            }

        }

        [TestMethod]
        public void SpearmanNullDistributionTest () {

            // pick independent distributions for x and y, which needn't be normal and needn't be related
            Distribution xDistrubtion = new UniformDistribution();
            Distribution yDistribution = new CauchyDistribution();
            Random rng = new Random(1);

            // generate bivariate samples of various sizes
            foreach (int n in TestUtilities.GenerateIntegerValues(4, 64, 8)) {
 
                Sample testStatistics = new Sample();
                Distribution testDistribution = null;

                for (int i = 0; i < 128; i++) {

                    BivariateSample sample = new BivariateSample();
                    for (int j = 0; j < n; j++) {
                        sample.Add(xDistrubtion.GetRandomValue(rng), yDistribution.GetRandomValue(rng));
                    }

                    TestResult result = sample.SpearmanRhoTest();
                    testStatistics.Add(result.Statistic);
                    testDistribution = result.Distribution;
                }

                TestResult r2 = testStatistics.KuiperTest(testDistribution);
                Console.WriteLine("n={0} P={1}", n, r2.LeftProbability);
                Assert.IsTrue(r2.RightProbability > 0.05);
            }

        }

        [TestMethod]
        public void KendallNullDistributionTest () {

            // pick independent distributions for x and y, which needn't be normal and needn't be related
            Distribution xDistrubtion = new LogisticDistribution();
            Distribution yDistribution = new ExponentialDistribution();
            Random rng = new Random(1000000000);

            // generate bivariate samples of various sizes
            //int n = 64; {
            foreach (int n in TestUtilities.GenerateIntegerValues(4, 64, 8)) {

                Sample testStatistics = new Sample();
                Distribution testDistribution = null;

                for (int i = 0; i < 128; i++) {

                    BivariateSample sample = new BivariateSample();
                    for (int j = 0; j < n; j++) {
                        sample.Add(xDistrubtion.GetRandomValue(rng), yDistribution.GetRandomValue(rng));
                    }

                    TestResult result = sample.KendallTauTest();
                    testStatistics.Add(result.Statistic);
                    testDistribution = result.Distribution;
                }

                TestResult r2 = testStatistics.KolmogorovSmirnovTest(testDistribution);
                Console.WriteLine("n={0} P={1}", n, r2.LeftProbability);
                //Assert.IsTrue(r2.RightProbability > 0.05);
            }

        }

    }
}
