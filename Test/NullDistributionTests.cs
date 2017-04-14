using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Functions;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace Test {

    [TestClass]
    public class NullDistributionTests {

        [TestMethod]
        public void KolmogorovNullDistributionTest () {

            // The distribution is irrelevent; pick one at random 
            ContinuousDistribution sampleDistribution = new LognormalDistribution();

            // Loop over various sample sizes
            foreach (int n in TestUtilities.GenerateIntegerValues(2, 128, 16)) {

                // Create a sample to hold the KS statistics
                Sample testStatistics = new Sample();
                // and a variable to hold the claimed null distribution, which should be the same for each test
                ContinuousDistribution nullDistribution = null;

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
            ContinuousDistribution sampleDistribution = new NormalDistribution();

            // Loop over various sample sizes
            foreach (int n in TestUtilities.GenerateIntegerValues(2, 128, 16)) {

                // Create a sample to hold the KS statistics
                Sample testStatistics = new Sample();
                // and a variable to hold the claimed null distribution, which should be the same for each test
                ContinuousDistribution nullDistribution = null;

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
        public void TwoSampleKolmogorovNullDistributionTest () {

            ContinuousDistribution population = new ExponentialDistribution();

            int[] sizes = new int[] { 23, 30, 175 };

            foreach (int na in sizes) {
                foreach (int nb in sizes) {
                    Console.WriteLine("{0} {1}", na, nb);

                    Sample d = new Sample();
                    ContinuousDistribution nullDistribution = null;
                    for (int i = 0; i < 128; i++) {

                        Sample a = TestUtilities.CreateSample(population, na, 31415 + na + i);
                        Sample b = TestUtilities.CreateSample(population, nb, 27182 + nb + i);

                        TestResult r = Sample.KolmogorovSmirnovTest(a, b);
                        d.Add(r.Statistic);
                        nullDistribution = r.Distribution;

                    }
                    // Only do full KS test if the number of bins is larger than the sample size, otherwise we are going to fail
                    // because the KS test detects the granularity of the distribution
                    TestResult mr = d.KolmogorovSmirnovTest(nullDistribution);
                    Console.WriteLine(mr.LeftProbability);
                    if (AdvancedIntegerMath.LCM(na, nb) > d.Count) Assert.IsTrue(mr.LeftProbability < 0.99);
                    // But always test that mean and standard deviation are as expected
                    Console.WriteLine("{0} {1}", nullDistribution.Mean, d.PopulationMean.ConfidenceInterval(0.99));
                    Assert.IsTrue(d.PopulationMean.ConfidenceInterval(0.99).ClosedContains(nullDistribution.Mean));
                    Console.WriteLine("{0} {1}", nullDistribution.StandardDeviation, d.PopulationStandardDeviation.ConfidenceInterval(0.99));
                    Assert.IsTrue(d.PopulationStandardDeviation.ConfidenceInterval(0.99).ClosedContains(nullDistribution.StandardDeviation));
                    Console.WriteLine("{0} {1}", nullDistribution.CentralMoment(3), d.PopulationCentralMoment(3).ConfidenceInterval(0.99));
                    //Assert.IsTrue(d.PopulationMomentAboutMean(3).ConfidenceInterval(0.99).ClosedContains(nullDistribution.MomentAboutMean(3)));

                    //Console.WriteLine("m {0} {1}", nullDistribution.Mean, d.PopulationMean);
                }
            }

        }

        [TestMethod]
        public void SpearmanNullDistributionTest () {

            // pick independent distributions for x and y, which needn't be normal and needn't be related
            ContinuousDistribution xDistrubtion = new UniformDistribution();
            ContinuousDistribution yDistribution = new CauchyDistribution();
            Random rng = new Random(1);

            // generate bivariate samples of various sizes
            foreach (int n in TestUtilities.GenerateIntegerValues(4, 64, 8)) {
 
                Sample testStatistics = new Sample();
                ContinuousDistribution testDistribution = null;

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

                Assert.IsTrue(testStatistics.PopulationMean.ConfidenceInterval(0.99).ClosedContains(testDistribution.Mean));
                Assert.IsTrue(testStatistics.PopulationVariance.ConfidenceInterval(0.99).ClosedContains(testDistribution.Variance));

            }

        }

        [TestMethod]
        public void KendallNullDistributionTest () {

            // pick independent distributions for x and y, which needn't be normal and needn't be related
            ContinuousDistribution xDistrubtion = new LogisticDistribution();
            ContinuousDistribution yDistribution = new ExponentialDistribution();
            Random rng = new Random(314159265);

            // generate bivariate samples of various sizes
            //int n = 64; {
            foreach (int n in TestUtilities.GenerateIntegerValues(4, 64, 8)) {

                Sample testStatistics = new Sample();
                ContinuousDistribution testDistribution = null;

                for (int i = 0; i < 128; i++) {

                    BivariateSample sample = new BivariateSample();
                    for (int j = 0; j < n; j++) {
                        sample.Add(xDistrubtion.GetRandomValue(rng), yDistribution.GetRandomValue(rng));
                    }

                    TestResult result = sample.KendallTauTest();
                    testStatistics.Add(result.Statistic);
                    testDistribution = result.Distribution;
                }

                //TestResult r2 = testStatistics.KolmogorovSmirnovTest(testDistribution);
                //Console.WriteLine("n={0} P={1}", n, r2.LeftProbability);
                //Assert.IsTrue(r2.RightProbability > 0.05);

                Console.WriteLine("{0} {1}", testStatistics.PopulationVariance, testDistribution.Variance);
                Assert.IsTrue(testStatistics.PopulationMean.ConfidenceInterval(0.95).ClosedContains(testDistribution.Mean));
                Assert.IsTrue(testStatistics.PopulationVariance.ConfidenceInterval(0.95).ClosedContains(testDistribution.Variance));
            }

        }

    }
}
