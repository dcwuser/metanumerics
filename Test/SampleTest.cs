using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.Generic;
using System.Collections;

using Meta.Numerics;
using Meta.Numerics.Statistics;

namespace Test {
    
    
    /// <summary>
    ///This is a test class for SampleTest and is intended
    ///to contain all SampleTest Unit Tests
    ///</summary>
    [TestClass()]
    public class SampleTest {


        private TestContext testContextInstance;

        /// <summary>
        ///Gets or sets the test context which provides
        ///information about and functionality for the current test run.
        ///</summary>
        public TestContext TestContext {
            get {
                return testContextInstance;
            }
            set {
                testContextInstance = value;
            }
        }

        #region Additional test attributes
        // 
        //You can use the following additional attributes as you write your tests:
        //
        //Use ClassInitialize to run code before running the first test in the class
        //[ClassInitialize()]
        //public static void MyClassInitialize(TestContext testContext)
        //{
        //}
        //
        //Use ClassCleanup to run code after all tests in a class have run
        //[ClassCleanup()]
        //public static void MyClassCleanup()
        //{
        //}
        //
        //Use TestInitialize to run code before running each test
        //[TestInitialize()]
        //public void MyTestInitialize()
        //{
        //}
        //
        //Use TestCleanup to run code after each test has run
        //[TestCleanup()]
        //public void MyTestCleanup()
        //{
        //}
        //
        #endregion

        [TestMethod]
        public void SampleManipulationsTest () {
            
            // create a sample
            double[] data = new double[] { -1.1, 2.2, -3.3, 4.4 };
            Sample sample = new Sample(data);

            // check the length
            Assert.IsTrue(sample.Count == data.Length);

            // add a datum and check the length
            sample.Add(5.5);
            Assert.IsTrue(sample.Count == data.Length + 1);

            // check wether an elements exists, remove it, check the length, check that it doesn't exist
            Assert.IsTrue(sample.Contains(2.2));
            Assert.IsTrue(sample.Remove(2.2));
            Assert.IsTrue(sample.Count == data.Length);
            Assert.IsFalse(sample.Contains(2.2));

            // clear the sample and check the length
            sample.Clear();
            Assert.IsTrue(sample.Count == 0);

        }

        private static Sample CreateSample (Distribution distribution, int count) {
            return (CreateSample(distribution, count, 1));
        }

        private static Sample CreateSample (Distribution distribution, int count, int seed) {

            Sample sample = new Sample();

            Random rng = new Random(seed);
            for (int i = 0; i < count; i++) {
                double x = distribution.InverseLeftProbability(rng.NextDouble());
                sample.Add(x);
            }

            return (sample);
        }

        private Distribution[] distributions = new Distribution[] {
            new UniformDistribution(Interval.FromEndpoints(-2.0,1.0)),
            new NormalDistribution(3.0,2.0),
            new ExponentialDistribution(2.0),
        };

        [TestMethod]
        public void SampleMomentsTest () {
            foreach (Distribution distribution in distributions) {
                //Distribution distribution = new NormalDistribution(3.0, 2.0);
                Sample sample = CreateSample(distribution, 100);

                Assert.IsTrue(sample.Count == 100);

                UncertainValue m = sample.PopulationMean;
                Interval mi = m.ConfidenceInterval(0.95);
                Assert.IsTrue(mi.ClosedContains(distribution.Mean));

                UncertainValue s = sample.PopulationStandardDeviation;
                Interval si = s.ConfidenceInterval(0.95);
                Assert.IsTrue(si.ClosedContains(distribution.StandardDeviation));

                for (int n = 3; n < 4; n++) {
                    UncertainValue c = sample.PopulationMomentAboutMean(3);
                    Interval ci = c.ConfidenceInterval(0.95);
                    Assert.IsTrue(ci.ClosedContains(distribution.MomentAboutMean(n)));

                    UncertainValue r = sample.PopulationMoment(n);
                    Interval ri = r.ConfidenceInterval(0.95);
                    Assert.IsTrue(ri.ClosedContains(distribution.Moment(n)));
                }
            }

        }

        [TestMethod]
        public void SampleInterquartileRangeTest () {
            foreach (Distribution distribution in distributions) {
                //Distribution distribution = new NormalDistribution(3.0, 2.0);
                Sample sample = CreateSample(distribution, 100);

                Interval iqr = sample.InterquartileRange;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(iqr.LeftEndpoint, sample.InverseLeftProbability(0.25)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(iqr.RightEndpoint, sample.InverseLeftProbability(0.75)));
            }

        }

        [TestMethod]
        public void SampleFitNormalTest () {

            // pick mu >> sigma so that we get no negative values;
            // otherwise the attempt to fit to an exponential will fail
            Distribution distribution = new NormalDistribution(5.0, 2.0);
            Sample sample = CreateSample(distribution, 100);

            // fit to normal should be good
            FitResult nfit = sample.FitToNormalDistribution();
            Console.WriteLine("P_n = {0}", nfit.GoodnessOfFit.LeftProbability);
            Assert.IsTrue(nfit.GoodnessOfFit.LeftProbability < 0.95, String.Format("P_n = {0}", nfit.GoodnessOfFit.LeftProbability));

            // fit to exponential should be bad
            FitResult efit = sample.FitToExponentialDistribution();
            Console.WriteLine("P_e = {0}", efit.GoodnessOfFit.LeftProbability);
            Assert.IsTrue(efit.GoodnessOfFit.LeftProbability > 0.95, String.Format("P_e = {0}", efit.GoodnessOfFit.LeftProbability));

        }

        [TestMethod]
        public void SampleFitExponentialTest () {

            Distribution distribution = new ExponentialDistribution(5.0);
            Sample sample = CreateSample(distribution, 100);

            // fit to normal should be bad
            FitResult nfit = sample.FitToNormalDistribution();
            Console.WriteLine("P_n = {0}", nfit.GoodnessOfFit.LeftProbability);
            Assert.IsTrue(nfit.GoodnessOfFit.LeftProbability > 0.95, String.Format("P_n = {0}", nfit.GoodnessOfFit.LeftProbability));

            // fit to exponential should be good
            FitResult efit = sample.FitToExponentialDistribution();
            Console.WriteLine("P_e = {0}", efit.GoodnessOfFit.LeftProbability);
            Assert.IsTrue(efit.GoodnessOfFit.LeftProbability < 0.95, String.Format("P_e = {0}", efit.GoodnessOfFit.LeftProbability));

        }

        [TestMethod]
        public void SampleFitExponentialUncertaintyTest () {

            // check that the uncertainty in reported fit parameters is actually meaningful
            // it should be the standard deviation of fit parameter values in a sample of many fits

            // define a population distribution 
            Distribution distribution = new ExponentialDistribution(4.0);

            // draw a lot of samples from it; fit each sample and
            // record the reported parameter value and error of each
            Sample values = new Sample();
            Sample uncertainties = new Sample();
            for (int i = 0; i < 50; i++) {
                Sample sample = CreateSample(distribution, 10, i);
                FitResult fit = sample.FitToExponentialDistribution();
                UncertainValue lambda = fit.Parameter(0);
                values.Add(lambda.Value);
                uncertainties.Add(lambda.Uncertainty);
            }

            Console.WriteLine(uncertainties.Mean);
            Console.WriteLine(values.PopulationStandardDeviation);

            // the reported errors should agree with the standard deviation of the reported parameters
            Assert.IsTrue(values.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(uncertainties.Mean));

        }

        [TestMethod]
        public void SampleFitChiSquaredTest () {

            Distribution distribution = new ChiSquaredDistribution(4);
            Sample sample = CreateSample(distribution, 100);

            // fit to normal should be bad
            // this is harder than others, because a chi^2 isn't so very different from a normal; to help, increse N or decrease vu
            FitResult nfit = sample.FitToNormalDistribution();
            Console.WriteLine("P_n = {0}", nfit.GoodnessOfFit.LeftProbability);
            Assert.IsTrue(nfit.GoodnessOfFit.LeftProbability > 0.95, String.Format("P_n = {0}", nfit.GoodnessOfFit.LeftProbability));

            // fit to exponential should also be bad
            FitResult efit = sample.FitToExponentialDistribution();
            Console.WriteLine("P_e = {0}", efit.GoodnessOfFit.LeftProbability);
            Assert.IsTrue(efit.GoodnessOfFit.LeftProbability > 0.95, String.Format("P_e = {0}", efit.GoodnessOfFit.LeftProbability));

        }

        [TestMethod]
        public void SampleTTestTest () {

            // start with a normally distributed population
            Distribution xDistribution = new NormalDistribution(2.0, 3.0);

            // draw 100 samples from it and compute the t statistic for each
            Sample tSample = new Sample();
            for (int i = 0; i < 50; i++) {
                Sample xSample = CreateSample(xDistribution, 10, i);
                TestResult tResult = xSample.StudentTTest(2.0);
                double t = tResult.Statistic;
                Console.WriteLine("t = {0}", t);
                tSample.Add(t);
            }

            // sanity check our sample of t's
            Assert.IsTrue(tSample.Count == 50);

            // check that the t statistics are distributed as expected
            Distribution tDistribution = new StudentDistribution(10);

            // check on the mean
            Console.WriteLine("m = {0} vs. {1}", tSample.PopulationMean, tDistribution.Mean);
            Assert.IsTrue(tSample.PopulationMean.ConfidenceInterval(0.95).ClosedContains(tDistribution.Mean), String.Format("{0} vs. {1}", tSample.PopulationMean, tDistribution.Mean));

            // check on the standard deviation
            Console.WriteLine("s = {0} vs. {1}", tSample.PopulationStandardDeviation, tDistribution.StandardDeviation);
            Assert.IsTrue(tSample.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(tDistribution.StandardDeviation));

            // do a KS test
            TestResult ksResult = tSample.KolmogorovSmirnovTest(tDistribution);
            Assert.IsTrue(ksResult.LeftProbability < 0.95);

            // check that we can distinguish the t distribution from a normal distribution?
        }

        [TestMethod]
        public void SampleComparisonTest () {

            // create one set of samples from our distributions
            Sample[] aSamples = new Sample[distributions.Length];
            for (int i = 0; i < distributions.Length; i++) {
                aSamples[i] = CreateSample(distributions[i], 40, 1);
            }

            // create another set
            Sample[] bSamples = new Sample[distributions.Length];
            for (int i = 0; i < distributions.Length; i++) {
                bSamples[i] = CreateSample(distributions[i], 80, 2);
            }

            KolmogorovDistribution kd = new KolmogorovDistribution();
            Console.WriteLine("P={0} => D={1}", 0.50, kd.InverseLeftProbability(0.50));
            Console.WriteLine("P={0} => D={1}", 0.90, kd.InverseLeftProbability(0.90));
            Console.WriteLine("P={0} => D={1}", 0.95, kd.InverseLeftProbability(0.95));
            Console.WriteLine("P={0} => D={1}", 0.99, kd.InverseLeftProbability(0.99));

            // cross-test using KS; like samples should agree and unlike samples should be distinguished
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {

                    //aSamples[0] = new Sample(new double[] { 10, 19, 15, 20, 12, 8, 15, 21 });
                    //bSamples[0] = new Sample(new double[] { 15, 22, 17, 9, 12, 10, 29, 11, 25, 31 });

                    //foreach (double datum in aSamples[i]) Console.WriteLine("a={0}", datum);
                    //foreach (double datum in bSamples[j]) Console.WriteLine("b={0}", datum);

                    TestResult result = aSamples[i].KolmogorovSmirnovTest(bSamples[j]);
                    Console.WriteLine("{0} v. {1}: D={2} P={3}", i, j, result.Statistic, result.LeftProbability);
                    if (i == j) {
                        Assert.IsTrue(result.LeftProbability < 0.90);
                    } else {
                        Assert.IsTrue(result.LeftProbability > 0.90);
                    }

                    // the order shouldn't matter
                    TestResult reverse = bSamples[j].KolmogorovSmirnovTest(aSamples[i]);
                    Assert.IsTrue(reverse.Statistic == result.Statistic);
                    Assert.IsTrue(reverse.RightProbability == result.RightProbability);

                }
            }

        }

        [TestMethod]
        public void SampleKolmogorovSmirnovTest () {

            // this test has a whiff of meta-statistics about it:
            // we want to make sure that the KS test statistic D is distributed according to the KS
            // distribution; to do this, we create a distributions of D statistics and KS test it
            // against the KS distribution

            // start with any 'ol underlying distribution
            Distribution distribution = new UniformDistribution(Interval.FromEndpoints(-2.0, 4.0));

            // generate some samples from it, and for each one get a D statistic from a KS test
            Sample DSample = new Sample();
            Distribution DDistribution = null;
            for (int i = 0; i < 50; i++) {
                Sample sample = CreateSample(distribution, 10, i);
                TestResult ks = sample.KolmogorovSmirnovTest(distribution);
                double D = ks.Statistic;
                Console.WriteLine("D = {0}", D);
                DSample.Add(D);
                DDistribution = ks.Distribution;
            }

            // check that the Ds are distributed as expected
            //Distribution DDistribution = new KolmogorovDistribution();

            // check on the mean
            Console.WriteLine("m = {0} vs. {1}", DSample.PopulationMean, DDistribution.Mean);
            Assert.IsTrue(DSample.PopulationMean.ConfidenceInterval(0.95).ClosedContains(DDistribution.Mean), String.Format("{0} vs. {1}", DSample.PopulationMean, DDistribution.Mean));

            // check on the standard deviation
            Console.WriteLine("s = {0} vs. {1}", DSample.PopulationStandardDeviation, DDistribution.StandardDeviation);
            Assert.IsTrue(DSample.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(DDistribution.StandardDeviation));

        }


    }
}
