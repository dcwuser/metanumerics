using System;
using System.Collections.Generic;
using System.Linq;

using TestClassAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestClassAttribute;
using TestMethodAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestMethodAttribute;
using ExpectedExceptionAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.ExpectedExceptionAttribute;
using Assert = Microsoft.VisualStudio.TestTools.UnitTesting.Assert;

using Meta.Numerics;
using Meta.Numerics.Data;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace Test {
    
    [TestClass]
    public class SampleTest {

        [TestMethod]
        public void SampleManipulations () {
            
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
            Assert.IsFalse(sample.Remove(2.2));

            // clear the sample and check the length
            sample.Clear();
            Assert.IsTrue(sample.Count == 0);

        }

        [TestMethod]
        public void SampleCopy () {

            // test independency of copy

            Sample sample1 = new Sample();
            sample1.Add(1.0, 2.0);

            Sample sample2 = sample1.Copy();
            sample2.Add(3.0, 4.0);

            Assert.IsTrue(sample1.Count == 2);
            Assert.IsTrue(sample2.Count == 4);

        }

        [TestMethod]
        public void LowSampleMoments () {

            // this is designed to test that means and standard deviations change properly when values
            // are added and removed

            Sample s = new Sample();
            s.Add(2.0);
            s.Add(4.0);
            s.Add(6.0);

            // set is 2, 4, 6: M = 4, V = (2^2 + 2^2) / 3
            Assert.IsTrue(s.Count == 3);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(s.Mean, 4.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(s.Variance, 8.0 / 3.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(s.StandardDeviation, Math.Sqrt(8.0 / 3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(s.CorrectedStandardDeviation, Math.Sqrt(8.0 / 2.0)));

            s.Remove(2.0);

            // set is 4, 6: M = 5, V = (1^2 + 1^2) / 2
            Assert.IsTrue(s.Count == 2);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(s.Mean, 5.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(s.Variance, 2.0 / 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(s.StandardDeviation, Math.Sqrt(2.0 / 2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(s.CorrectedStandardDeviation, Math.Sqrt(2.0 / 1.0)));

        }

        public static Sample CreateSample (ContinuousDistribution distribution, int count) {
            return (CreateSample(distribution, count, 1));
        }

        public static Sample CreateSample (ContinuousDistribution distribution, int count, int seed) {

            Sample sample = new Sample();

            Random rng = new Random(seed);
            for (int i = 0; i < count; i++) {
                double x = distribution.GetRandomValue(rng);
                sample.Add(x);
            }

            return (sample);
        }

        private ContinuousDistribution[] distributions = new ContinuousDistribution[] {
            new UniformDistribution(Interval.FromEndpoints(-2.0,1.0)),
            new NormalDistribution(3.0,2.0),
            new ExponentialDistribution(2.0),
        };

        [TestMethod]
        public void SampleMoments () {
            foreach (ContinuousDistribution distribution in distributions) {

                int n = 256;
                Random rng = new Random(1);
                List<double> sample = new List<double>(distribution.GetRandomValues(rng, n));

                Assert.IsTrue(sample.Count == n);

                Assert.IsTrue(sample.PopulationMean().ConfidenceInterval(0.99).ClosedContains(distribution.Mean));

                Assert.IsTrue(sample.PopulationStandardDeviation().ConfidenceInterval(0.99).ClosedContains(distribution.StandardDeviation));

                for (int r = 0; r < 8; r++) {
                    Assert.IsTrue(sample.PopulationCentralMoment(r).ConfidenceInterval(0.99).ClosedContains(distribution.CentralMoment(r)));
                    Assert.IsTrue(sample.PopulationRawMoment(r).ConfidenceInterval(0.99).ClosedContains(distribution.RawMoment(r)));
                }
            }

        }


        [TestMethod]
        public void SamplePopulationMomentEstimateVariances () {

            ContinuousDistribution d = new LognormalDistribution();

            // for various sample sizes...
            foreach (int n in TestUtilities.GenerateIntegerValues(4, 32, 8)) {

                Console.WriteLine("n={0}", n);

                // we are going to store values for a bunch of estimators and their uncertainties
                MultivariateSample estimates = new MultivariateSample("M1", "C2", "C3", "C4");
                MultivariateSample variances = new MultivariateSample("M1", "C2", "C3", "C4");

                // create a bunch of samples
                for (int i = 0; i < 4096; i++) {

                    Sample s = TestUtilities.CreateSample(d, n, 16384 * n + i + 1);

                    UncertainValue M1 = s.PopulationMean;
                    UncertainValue C2 = s.PopulationVariance;
                    UncertainValue C3 = s.PopulationCentralMoment(3);
                    UncertainValue C4 = s.PopulationCentralMoment(4);
                    estimates.Add(M1.Value, C2.Value, C3.Value, C4.Value);
                    variances.Add(MoreMath.Sqr(M1.Uncertainty), MoreMath.Sqr(C2.Uncertainty), MoreMath.Sqr(C3.Uncertainty), MoreMath.Sqr(C4.Uncertainty));

                }

                // the claimed variance should agree with the measured variance of the estimators
                for (int c = 0; c < estimates.Dimension; c++) {
                    Console.WriteLine("{0} {1} {2}", estimates.Column(c).Name, estimates.Column(c).PopulationVariance, variances.Column(c).Mean);
                    Assert.IsTrue(estimates.Column(c).PopulationVariance.ConfidenceInterval(0.95).ClosedContains(variances.Column(c).Mean));
                }

            }

        }

        [TestMethod]
        public void SampleMedian () {
            Sample sample = new Sample();


            sample.Add(2.0, 1.0);
            Assert.IsTrue(sample.Minimum == 1.0);
            Assert.IsTrue(sample.Median == 1.5);
            Assert.IsTrue(sample.Maximum == 2.0);

            sample.Add(3.0);
            Assert.IsTrue(sample.Minimum == 1.0);
            Assert.IsTrue(sample.Median == 2.0);
            Assert.IsTrue(sample.Maximum == 3.0);
        }

        [TestMethod]
        public void SampleInterquartileRange () {
            foreach (ContinuousDistribution distribution in distributions) {
                
                Sample sample = CreateSample(distribution, 100);

                Interval iqr = sample.InterquartileRange;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(iqr.LeftEndpoint, sample.InverseLeftProbability(0.25)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(iqr.RightEndpoint, sample.InverseLeftProbability(0.75)));
            }

        }


        [TestMethod]
        public void SampleTrimean () {

            // Examples from https://explorable.com/trimean
            // Indexes 8 * 0.25 = 2, 8 * 0.50 = 4, 8 * 0.75 = 6
            // For d1, (161 + 2 * 166 + 171) / 4 = 166
            // For d2, (163 + 2 * 166 + 181) / 4 = 169

            double[] d1 = new double[] { 155, 158, 161, 162, 166, 170, 171, 174, 179 };
            double[] d2 = new double[] { 162, 162, 163, 165, 166, 175, 181, 186, 192 };

            Assert.IsTrue(TestUtilities.IsNearlyEqual(d1.Trimean(), 166));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(d2.Trimean(), 169));
        }

        [TestMethod]
        public void ZeroSampleProperties () {
            double[] sample = new double[] { };
            Assert.IsTrue(Double.IsNaN(sample.Minimum()));
            Assert.IsTrue(Double.IsNaN(sample.Maximum()));
            Assert.IsTrue(Double.IsNaN(sample.Mean()));
            //Assert.IsTrue(sample.RawMoment(0) == 1.0);
            //Assert.IsTrue(sample.CentralMoment(0) == 1.0);
        }

        [TestMethod]
        public void OneSampleProperties () {
            double[] sample = new double[] { 2.0 };
            Assert.IsTrue(sample.RawMoment(0) == 1.0);
            Assert.IsTrue(sample.CentralMoment(0) == 1.0);
            Assert.IsTrue(sample.RawMoment(1) == sample[0]);
            Assert.IsTrue(sample.CentralMoment(1) == 0.0);
            Assert.IsTrue(sample.Mean() == sample[0]);
            Assert.IsTrue(sample.Minimum() == sample[0]);
            Assert.IsTrue(sample.Maximum() == sample[0]);
            Assert.IsTrue(sample.Median() == sample[0]);
        }

        [TestMethod]
        public void LognormalFit () {

            LognormalDistribution distribution = new LognormalDistribution(1.0, 2.0);
            Sample sample = CreateSample(distribution, 100);

            // fit to normal should be bad
            NormalFitResult nfit = NormalDistribution.FitToSample(sample);
            Assert.IsTrue(nfit.GoodnessOfFit.Probability < 0.05);
            ValidateParameterCollection(nfit.Parameters);

            // fit to lognormal should be good
            LognormalFitResult efit = LognormalDistribution.FitToSample(sample);
            Assert.IsTrue(efit.GoodnessOfFit.Probability > 0.05);
            Assert.IsTrue(efit.Mu.ConfidenceInterval(0.95).ClosedContains(distribution.Mu));
            Assert.IsTrue(efit.Sigma.ConfidenceInterval(0.95).ClosedContains(distribution.Sigma));
            ValidateParameterCollection(efit.Parameters);

        }

        [TestMethod]
        public void BetaFit () {
            // Create a very non-normal Beta distribution
            BetaDistribution B = new BetaDistribution(0.25, 3.0);

            // Sample from it
            Sample S = CreateSample(B, 128);

            // Fit the sample to a beta distribution
            // The fit should agree with the population and the fit should be good
            BetaFitResult RB = BetaDistribution.FitToSample(S);
            Assert.IsTrue(RB.Alpha.ConfidenceInterval(0.99).ClosedContains(B.Alpha));
            Assert.IsTrue(RB.Beta.ConfidenceInterval(0.99).ClosedContains(B.Beta));
            Assert.IsTrue(RB.GoodnessOfFit.Probability > 0.05);
            ValidateParameterCollection(RB.Parameters);

            // Fit to a normal distribution should be bad
            NormalFitResult RN = NormalDistribution.FitToSample(S);
            Assert.IsTrue(RN.GoodnessOfFit.Probability < 0.05);
            ValidateParameterCollection(RN.Parameters);
        }

        [TestMethod]
        public void BetaFitUncertainty () {

            // check that the uncertainty in reported fit parameters is actually meaningful
            // it should be the standard deviation of fit parameter values in a sample of many fits

            // define a population distribution 
            ContinuousDistribution distribution = new BetaDistribution(1.0 / 3.0, 2.0);

            // draw a lot of samples from it; fit each sample and
            // record the reported parameter value and error of each
            BivariateSample values = new BivariateSample();
            BivariateSample uncertainties = new BivariateSample();
            for (int i = 0; i < 50; i++) {
                Sample sample = CreateSample(distribution, 10, i);
                BetaFitResult fit = BetaDistribution.FitToSample(sample);
                UncertainValue a = fit.Alpha;
                UncertainValue b = fit.Beta;
                values.Add(a.Value, b.Value);
                uncertainties.Add(a.Uncertainty, b.Uncertainty);
            }

            // the reported errors should agree with the standard deviation of the reported parameters
            Assert.IsTrue(values.X.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(uncertainties.X.Mean));
            Assert.IsTrue(values.Y.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(uncertainties.Y.Mean));

        }

        [TestMethod]
        public void GammaFit () {

            // Create a gamma distribution
            GammaDistribution B = new GammaDistribution(2.75, 2.0);

            // Sample from it
            Sample S = CreateSample(B, 100);

            // Fit the sample to a gamma distribution
            // The fit should agree with the population and the fit should be good
            GammaFitResult R = GammaDistribution.FitToSample(S);
            Assert.IsTrue(R.Shape.ConfidenceInterval(0.95).ClosedContains(B.Shape));
            Assert.IsTrue(R.Scale.ConfidenceInterval(0.95).ClosedContains(B.Scale));
            Assert.IsTrue(R.GoodnessOfFit.Probability > 0.05);
            ValidateParameterCollection(R.Parameters);

        }

        [TestMethod]
        public void GammaFitUncertainty () {

            // check that the uncertainty in reported fit parameters is actually meaningful
            // it should be the standard deviation of fit parameter values in a sample of many fits

            // define a population distribution 
            ContinuousDistribution distribution = new GammaDistribution(1.5, 2.0);

            // draw a lot of samples from it; fit each sample and
            // record the reported parameter value and error of each
            BivariateSample values = new BivariateSample();
            BivariateSample uncertainties = new BivariateSample();
            for (int i = 0; i < 100; i++) {
                Sample sample = CreateSample(distribution, 50, i);
                GammaFitResult fit = GammaDistribution.FitToSample(sample);
                UncertainValue a = fit.Parameters["Shape"].Estimate;
                UncertainValue b = fit.Parameters["Scale"].Estimate;
                values.Add(a.Value, b.Value);
                uncertainties.Add(a.Uncertainty, b.Uncertainty);
            }

            // the reported errors should agree with the standard deviation of the reported parameters
            Assert.IsTrue(values.X.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(uncertainties.X.Mean));
            Assert.IsTrue(values.Y.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(uncertainties.Y.Mean));

        }



        [TestMethod]
        public void SampleFitChiSquaredTest () {
            ContinuousDistribution distribution = new ChiSquaredDistribution(4);
            Sample sample = CreateSample(distribution, 100);

            // fit to normal should be bad
            // this is harder than others, because a chi^2 isn't so very different from a normal; to help, increase N or decrease nu
            NormalFitResult nfit = NormalDistribution.FitToSample(sample);
            Assert.IsTrue(nfit.GoodnessOfFit.Probability < 0.05);
            ValidateParameterCollection(nfit.Parameters);

            // fit to exponential should also be bad
            ExponentialFitResult efit = ExponentialDistribution.FitToSample(sample);
            Assert.IsTrue(efit.GoodnessOfFit.Probability < 0.05);
            ValidateParameterCollection(efit.Parameters);
        }

        [TestMethod]
        public void WeibullFitSimple () {

            foreach (double alpha in new double[] { 0.031, 0.13, 1.3, 3.1, 13.1, 33.1, 313.1, 1313.1 }) {
                WeibullDistribution W = new WeibullDistribution(2.2, alpha);
                Sample S = CreateSample(W, 50);
                WeibullFitResult R = WeibullDistribution.FitToSample(S);
                Assert.IsTrue(R.Scale.ConfidenceInterval(0.99).ClosedContains(W.Scale));
                Assert.IsTrue(R.Shape.ConfidenceInterval(0.99).ClosedContains(W.Shape));
                Assert.IsTrue(R.GoodnessOfFit.Probability > 0.01);
                ValidateParameterCollection(R.Parameters);
            }

        }

        [TestMethod]
        public void WeibullRngProblem () {

            Random rng = new Random(1);

            Sample s = new Sample();
            for (int i = 0; i < 50; i++) {
                s.Add(rng.NextDouble());
            }

            WeibullDistribution w = new WeibullDistribution(2.2, 0.031);

            s.Transform(x => w.Scale * Math.Pow(-Math.Log(1.0 - x), 1.0 / w.Shape));

        }

        [TestMethod]
        public void WeibullFitUncertainties () {

            // check that the uncertainty in reported fit parameters is actually meaningful
            // it should be the standard deviation of fit parameter values in a sample of many fits

            // define a population distribution 
            WeibullDistribution distribution = new WeibullDistribution(2.5, 1.5);

            // draw a lot of samples from it; fit each sample and
            // record the reported parameter value and error of each
            BivariateSample values = new BivariateSample();
            MultivariateSample uncertainties = new MultivariateSample(3);
            FrameTable table = new FrameTable();
            table.AddColumn<int>("Id");
            table.AddColumns<double>("Scale", "Shape", "ScaleVariance", "ShapeVariance", "ScaleShapeCovariance");

            for (int i = 0; i < 50; i++) {
                Sample sample = CreateSample(distribution, 10, i);
                WeibullFitResult fit = WeibullDistribution.FitToSample(sample);
                UncertainValue a = fit.Scale;
                UncertainValue b = fit.Shape;
                values.Add(a.Value, b.Value);
                uncertainties.Add(a.Uncertainty, b.Uncertainty, fit.Parameters.CovarianceMatrix[0,1]);
                table.AddRow(i, fit.Scale.Value, fit.Shape.Value,
                    fit.Parameters.CovarianceOf("Scale", "Scale"),
                    fit.Parameters.CovarianceOf("Shape", "Shape"),
                    fit.Parameters.CovarianceOf("Scale", "Shape")
                );
            }

            Assert.IsTrue(table["Scale"].As<double>().PopulationMean().ConfidenceInterval(0.99).ClosedContains(distribution.Scale));
            Assert.IsTrue(table["Shape"].As<double>().PopulationMean().ConfidenceInterval(0.99).ClosedContains(distribution.Shape));

            Assert.IsTrue(table["Scale"].As<double>().PopulationVariance().ConfidenceInterval(0.99).ClosedContains(table["ScaleVariance"].As<double>().Median()));
            Assert.IsTrue(table["Shape"].As<double>().PopulationVariance().ConfidenceInterval(0.99).ClosedContains(table["ShapeVariance"].As<double>().Median()));
            Assert.IsTrue(table["Scale"].As<double>().PopulationCovariance(table["Shape"].As<double>()).ConfidenceInterval(0.99).ClosedContains(table["ScaleShapeCovariance"].As<double>().Median()));

            // the reported errors should agree with the standard deviation of the reported parameters
            Assert.IsTrue(values.X.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(uncertainties.Column(0).Mean));
            Assert.IsTrue(values.Y.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(uncertainties.Column(1).Mean));
            //Console.WriteLine("{0} {1}", values.PopulationCovariance, uncertainties.Column(2).Mean);
            //Assert.IsTrue(values.PopulationCovariance.ConfidenceInterval(0.95).ClosedContains(uncertainties.Column(2).Mean));

        }

        [TestMethod]
        public void TTestDistribution () {

            // start with a normally distributed population
            ContinuousDistribution xDistribution = new NormalDistribution(2.0, 3.0);
            Random rng = new Random(1);

            // draw 100 samples from it and compute the t statistic for each
            Sample tSample = new Sample();
            for (int i = 0; i < 100; i++) {

                // each sample has 9 values
                List<double> xSample = TestUtilities.CreateDataSample(rng, xDistribution, 9).ToList();
                TestResult tResult = xSample.StudentTTest(2.0);
                Assert.IsTrue(tResult.Statistic.Name == "t");
                double t = tResult.Statistic.Value;
                tSample.Add(t);
            }

            // sanity check our sample of t's
            Assert.IsTrue(tSample.Count == 100);

            // check that the t statistics are distributed as expected
            ContinuousDistribution tDistribution = new StudentDistribution(9);

            // check on the mean
            Assert.IsTrue(tSample.PopulationMean.ConfidenceInterval(0.95).ClosedContains(tDistribution.Mean));

            // check on the standard deviation
            Assert.IsTrue(tSample.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(tDistribution.StandardDeviation));

            // do a KS test
            TestResult ksResult = tSample.KolmogorovSmirnovTest(tDistribution);
            Assert.IsTrue(ksResult.Probability > 0.05);

            // check that we can distinguish the t distribution from a normal distribution?
        }

        [TestMethod]
        public void SignTestDistribution () {

            // start with a non-normally distributed population
            ContinuousDistribution xDistribution = new ExponentialDistribution();
            Random rng = new Random(1);

            // draw 100 samples from it and compute the t statistic for each
            Sample wSample = new Sample();
            for (int i = 0; i < 100; i++) {

                // each sample has 8 observations
                Sample xSample = new Sample();
                for (int j = 0; j < 8; j++) { xSample.Add(xDistribution.GetRandomValue(rng)); }
                //Sample xSample = CreateSample(xDistribution, 8, i);
                TestResult wResult = xSample.SignTest(xDistribution.Median);
                double W = wResult.Statistic.Value;
                //Console.WriteLine("W = {0}", W);
                wSample.Add(W);
            }

            // sanity check our sample of t's
            Assert.IsTrue(wSample.Count == 100);

            // check that the t statistics are distributed as expected
            DiscreteDistribution wDistribution = new BinomialDistribution(0.5, 8);

            // check on the mean
            Console.WriteLine("m = {0} vs. {1}", wSample.PopulationMean, wDistribution.Mean);
            Assert.IsTrue(wSample.PopulationMean.ConfidenceInterval(0.95).ClosedContains(wDistribution.Mean));

            // check on the standard deviation
            Console.WriteLine("s = {0} vs. {1}", wSample.PopulationStandardDeviation, wDistribution.StandardDeviation);
            Assert.IsTrue(wSample.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(wDistribution.StandardDeviation));

            // check on the skew
            Console.WriteLine("t = {0} vs. {1}", wSample.PopulationCentralMoment(3), wDistribution.CentralMoment(3));
            Assert.IsTrue(wSample.PopulationCentralMoment(3).ConfidenceInterval(0.95).ClosedContains(wDistribution.CentralMoment(3)));

            // check on the kuritosis
            Console.WriteLine("u = {0} vs. {1}", wSample.PopulationCentralMoment(4), wDistribution.CentralMoment(4));
            Assert.IsTrue(wSample.PopulationCentralMoment(4).ConfidenceInterval(0.95).ClosedContains(wDistribution.CentralMoment(4)));

            // KS tests are only for continuous distributions            

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
            
            // cross-test using KS; like samples should agree and unlike samples should be distinguished
            for (int i = 0; i < aSamples.Length; i++) {
                for (int j = 0; j < bSamples.Length; j++) {

                    TestResult result = Sample.KolmogorovSmirnovTest(aSamples[i], bSamples[j]);
                    if (i == j) {
                        Assert.IsTrue(result.Probability > 0.05);
                    } else {
                        Assert.IsTrue(result.Probability < 0.05);
                    }

                    // the order shouldn't matter
                    TestResult reverse = Sample.KolmogorovSmirnovTest(bSamples[j], aSamples[i]);
                    Assert.IsTrue(reverse.Statistic.Value == result.Statistic.Value);
                    Assert.IsTrue(reverse.Probability == result.Probability);

                }
            }

        }


        [TestMethod]
        public void SampleKolmogorovSmirnovTest () {

            // This test has a whiff of "meta" about it.
            // We want to make sure that the KS test statistic D is distributed according to the Kolmogorov
            // distribution. To do this, we create a sample of D statistics and do KS/Kuiper tests
            // comparing it to the claimed Kolmogorov distribution.

            // Start with any old underlying distribution
            ContinuousDistribution distribution = new UniformDistribution(Interval.FromEndpoints(-2.0, 4.0));

            // Generate some samples from it, and for each one get a D statistic from a KS test comparing
            // the sample to the source distribution.
            // The sample sizes must be such that our approximation to the D-distribution is good.
            // At the moment we are good when the size is small enough that we can compute the
            // exact distribution and when the size is large enough that our asymptotic approximation
            // is accurate. In the intermediate range it gets iffy.
            Random rng = new Random(1);
            double[] sizes = new double[] { 7, 55, 333 };
            foreach (int n in sizes) {

                List<double> DSample = new List<double>();
                ContinuousDistribution DDistribution = null;
                for (int i = 0; i < 32; i++) {
                    List<double> sample = TestUtilities.CreateDataSample(rng, distribution, n).ToList();
                    TestResult ks = sample.KolmogorovSmirnovTest(distribution);
                    double D = ks.Statistic.Value;
                    DSample.Add(D);
                    DDistribution = ks.Statistic.Distribution;
                }

                // check on the mean
                Assert.IsTrue(DSample.PopulationMean().ConfidenceInterval(0.99).ClosedContains(DDistribution.Mean));

                // check on the standard deviation
                 Assert.IsTrue(DSample.PopulationStandardDeviation().ConfidenceInterval(0.99).ClosedContains(DDistribution.StandardDeviation));

                // do a KS test comparing the sample to the expected distribution
                TestResult kst = DSample.KolmogorovSmirnovTest(DDistribution);
                Assert.IsTrue(kst.Probability > 0.01);

                // do a Kuiper test comparing the sample to the expected distribution
                TestResult kut = DSample.KuiperTest(DDistribution);
                Assert.IsTrue(kut.Probability > 0.01);
            }
        }


        [TestMethod]
        public void SampleKuiperTest () {

            // This test has a whiff of "meta" about it.
            // We want to make sure that the KS test statistic D is distributed according to the Kolmogorov
            // distribution. To do this, we create a sample of D statistics and do KS/Kuiper tests
            // comparing it to the claimed Kolmogorov distribution.

            // Start with any old underlying distribution
            ContinuousDistribution distribution = new UniformDistribution(Interval.FromEndpoints(-2.0, 4.0));

            // Generate some samples from it, and for each one get a D statistic from a KS test comparing
            // the sample to the source distribution.
            // The sample sizes must be such that our approximation to the D-distribution is good.
            // At the moment we are good when the size is small enough that we can compute the
            // exact distribution and when the size is large enough that our asymptotic approximation
            // is accurate. In the intermediate range it gets iffy.
            Random rng = new Random(7);
            double[] sizes = new double[] { 6, 54, 321 };
            foreach (int n in sizes) {

                List<double> VSample = new List<double>();
                ContinuousDistribution VDistribution = null;
                for (int i = 0; i < 32; i++) {
                    List<double> sample = TestUtilities.CreateDataSample(rng, distribution, n).ToList();
                    TestResult kp = sample.KuiperTest(distribution);
                    double V = kp.Statistic.Value;
                    VSample.Add(V);
                    VDistribution = kp.Statistic.Distribution;
                }

                // check on the mean
                Assert.IsTrue(VSample.PopulationMean().ConfidenceInterval(0.99).ClosedContains(VDistribution.Mean));

                // check on the standard deviation
                Assert.IsTrue(VSample.PopulationStandardDeviation().ConfidenceInterval(0.99).ClosedContains(VDistribution.StandardDeviation));

                // do a KS test comparing the sample to the expected distribution
                TestResult kst = VSample.KolmogorovSmirnovTest(VDistribution);
                Assert.IsTrue(kst.Probability > 0.01);

                // do a Kuiper test comparing the sample to the expected distribution
                TestResult kut = VSample.KuiperTest(VDistribution);
                Assert.IsTrue(kut.Probability > 0.01);
            }
        }

        [TestMethod]
        public void MaximumLikelihoodNormalAgreement () {

            // create a normal sample
            double mu = -1.0;
            double sigma = 2.0;
            ContinuousDistribution dist = new NormalDistribution(mu, sigma);
            Random rng = new Random(1);
            List<double> sample = TestUtilities.CreateDataSample(rng, dist, 256).ToList();

            // do a maximum likelihood fit to the normal form
            Func<IReadOnlyDictionary<string, double>, ContinuousDistribution> factory = dict => new NormalDistribution(dict["Mu"], dict["Sigma"]);
            Dictionary<string, double> start = new Dictionary<string, double>() { { "Mu", 0.0 }, { "Sigma", 1.0 } };
            DistributionFitResult<ContinuousDistribution> result  = sample.MaximumLikelihoodFit(factory, start);

            // the fit should find the right parameters
            Assert.IsTrue(result.GoodnessOfFit.Probability > 0.05);
            Assert.IsTrue(result.Parameters.Count == 2);
            Assert.IsTrue(result.Parameters["Mu"].Estimate.ConfidenceInterval(0.99).ClosedContains(mu));
            Assert.IsTrue(result.Parameters["Sigma"].Estimate.ConfidenceInterval(0.99).ClosedContains(sigma));
        }

        [TestMethod]
        public void MaximumLikelihoodExponentialAgreement () {
            Sample sample = CreateSample(new ExponentialDistribution(3.0), 256);

            // Do an explicit exponential fit
            ExponentialFitResult fit1 = ExponentialDistribution.FitToSample(sample);

            // Do a maximum likelihood fit
            Func<IReadOnlyList<double>, ContinuousDistribution> factory = p => new ExponentialDistribution(p[0]);
            DistributionFitResult<ContinuousDistribution> fit2 = sample.MaximumLikelihoodFit(factory, new double[] { 1.0 });

            // Results should at least approximately agree
            Assert.IsTrue(TestUtilities.IsNearlyEqual(fit1.Parameters[0].Estimate.Value, fit2.Parameters[0].Estimate.Value, 1.0E-4));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(fit1.Parameters[0].Estimate.Uncertainty, fit2.Parameters[0].Estimate.Uncertainty, 1.0E-2));
        }

        /*
        [TestMethod]
        public void SampleMaximumLikelihoodFit () {

            // normal distribution
            Console.WriteLine("normal");

            double mu = -1.0;
            double sigma = 2.0;
            ContinuousDistribution nd = new NormalDistribution(mu, sigma);
            Sample ns = CreateSample(nd, 500);

            DistributionFitResult<ContinuousDistribution> nr = ns.MaximumLikelihoodFit((IReadOnlyList<double> p) => new NormalDistribution(p[0], p[1]), new double[] { mu + 1.0, sigma + 1.0 });

            Assert.IsTrue(nr.Dimension == 2);
            Assert.IsTrue(nr.Parameter(0).ConfidenceInterval(0.95).ClosedContains(mu));
            Assert.IsTrue(nr.Parameter(1).ConfidenceInterval(0.95).ClosedContains(sigma));

            NormalFitResult nr2 = NormalDistribution.FitToSample(ns);

            Console.WriteLine(nr.Covariance(0,1));

            // test analytic expression
            Assert.IsTrue(TestUtilities.IsNearlyEqual(nr.Parameter(0).Value, ns.Mean, Math.Sqrt(TestUtilities.TargetPrecision)));
            // we don't expect to be able to test sigma against analytic expression because ML value has known bias for finite sample size

            
            // exponential distribution
            
            Console.WriteLine("exponential");
            double em = 3.0;
            ContinuousDistribution ed = new ExponentialDistribution(em);
            Sample es = CreateSample(ed, 100);
            //FitResult er = es.MaximumLikelihoodFit(new ExponentialDistribution(em + 1.0));
            FitResult er = es.MaximumLikelihoodFit((IReadOnlyList<double> p) => new ExponentialDistribution(p[0]), new double[] { em + 1.0 });

            Console.WriteLine(er.Parameter(0));

            Assert.IsTrue(er.Dimension == 1);
            Assert.IsTrue(er.Parameter(0).ConfidenceInterval(0.95).ClosedContains(em));

            // test against analytic expression
            Assert.IsTrue(TestUtilities.IsNearlyEqual(er.Parameter(0).Value, es.Mean, Math.Sqrt(TestUtilities.TargetPrecision)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(er.Parameter(0).Uncertainty, es.Mean / Math.Sqrt(es.Count), Math.Sqrt(Math.Sqrt(TestUtilities.TargetPrecision))));
            
            

            // lognormal distribution
            Console.WriteLine("lognormal");

            double l1 = -4.0;
            double l2 = 5.0;

            ContinuousDistribution ld = new LognormalDistribution(l1, l2);
            Sample ls = CreateSample(ld, 100);
            //FitResult lr = ls.MaximumLikelihoodFit(new LognormalDistribution(l1 + 1.0, l2 + 1.0));
            FitResult lr = ls.MaximumLikelihoodFit((IReadOnlyList<double> p) => new LognormalDistribution(p[0], p[1]), new double[] { l1 + 1.0, l2 + 1.0 });

            Console.WriteLine(lr.Parameter(0));
            Console.WriteLine(lr.Parameter(1));
            Console.WriteLine(lr.Covariance(0, 1));

            Assert.IsTrue(lr.Dimension == 2);
            Assert.IsTrue(lr.Parameter(0).ConfidenceInterval(0.99).ClosedContains(l1));
            Assert.IsTrue(lr.Parameter(1).ConfidenceInterval(0.99).ClosedContains(l2));

            // weibull distribution
            Console.WriteLine("weibull");

            double w_scale = 4.0;
            double w_shape = 2.0;
            WeibullDistribution w_d = new WeibullDistribution(w_scale, w_shape);
            Sample w_s = CreateSample(w_d, 20);
            //FitResult w_r = w_s.MaximumLikelihoodFit(new WeibullDistribution(1.0, 0.5));
            FitResult w_r = w_s.MaximumLikelihoodFit((IReadOnlyList<double> p) => new WeibullDistribution(p[0], p[1]), new double[] { 2.0, 2.0 });

            Console.WriteLine(w_r.Parameter(0));
            Console.WriteLine(w_r.Parameter(1));
            Console.WriteLine(w_r.Covariance(0, 1));

            Assert.IsTrue(w_r.Parameter(0).ConfidenceInterval(0.95).ClosedContains(w_d.Scale));
            Assert.IsTrue(w_r.Parameter(1).ConfidenceInterval(0.95).ClosedContains(w_d.Shape));

            // logistic distribution
            Console.WriteLine("logistic");

            double logistic_m = -3.0;
            double logistic_s = 2.0;
            ContinuousDistribution logistic_distribution = new LogisticDistribution(logistic_m, logistic_s);
            Sample logistic_sample = CreateSample(logistic_distribution, 100);
            //FitResult logistic_result = logistic_sample.MaximumLikelihoodFit(new LogisticDistribution());
            FitResult logistic_result = logistic_sample.MaximumLikelihoodFit((IReadOnlyList<double> p) => new LogisticDistribution(p[0], p[1]), new double[] { 2.0, 3.0 });

            Console.WriteLine(logistic_result.Parameter(0));
            Console.WriteLine(logistic_result.Parameter(1));

            Assert.IsTrue(logistic_result.Dimension == 2);
            Assert.IsTrue(logistic_result.Parameter(0).ConfidenceInterval(0.95).ClosedContains(logistic_m));
            Assert.IsTrue(logistic_result.Parameter(1).ConfidenceInterval(0.95).ClosedContains(logistic_s));

        }
        */
        /*
        [TestMethod]
        public void SampleMannWhitneyComputationTest () {

            Sample sample1 = new Sample();
            sample1.Add(new double[] { 1.2, 1.3, 1.7, 2.4, 5.3 });

            Sample sample2 = new Sample();
            sample2.Add(new double[] { 1.0, 2.0, 3.0, 4.0 });

            TestResult result = sample1.MannWhitneyTest(sample2);
            Console.WriteLine(result.Statistic);
            Console.WriteLine(result.LeftProbability);
            Console.WriteLine(result.RightProbability);

        }
        */

        [TestMethod]
        public void SampleFisherFTest () {


            // Create 3 samples.
            // 1 and 2 have the same mean but different variances, the F test should catch the difference.
            // 1 and 3 have different means but the same variance, the F test should rule them equivalent.
            Sample sample1 = CreateSample(new NormalDistribution(1.0, 1.0), 20, 1);
            Sample sample2 = CreateSample(new NormalDistribution(1.0, 2.0), 20, 2);
            Sample sample3 = CreateSample(new NormalDistribution(3.0, 1.0), 20, 3);

            TestResult f12 = Sample.FisherFTest(sample1, sample2);
            TestResult f21 = Sample.FisherFTest(sample2, sample1);

            // sample 1 has a smaller variance
            Assert.IsTrue(f12.Statistic.Value < 1.0);

            // 1/2 is the inverse of 2/1
            Assert.IsTrue(TestUtilities.IsNearlyEqual(f12.Statistic, 1.0 / f21.Statistic));

            // the F test detects the difference between the variance of 1 and 2
            Assert.IsTrue(f12.Probability < 0.05);

            // the F test detects no difference between the variance of 1 and 3
            TestResult f13 = Sample.FisherFTest(sample1, sample3);
            Assert.IsTrue(f13.Probability > 0.05);
        }

        [TestMethod]
        public void SampleMannWhitneyTest () {

            // define two non-normal distributions
            ContinuousDistribution d1 = new ExponentialDistribution(2.0);
            ContinuousDistribution d2 = new ExponentialDistribution(3.0);

            // create three samples from them
            Sample s1a = CreateSample(d1, 15, 1);
            Sample s1b = CreateSample(d1, 30, 2);
            Sample s2 = CreateSample(d2, 60, 3);

            // Mann-Whitney test 1a vs. 1b; they should not be distinguished
            TestResult rab = Sample.MannWhitneyTest(s1a, s1b);
            Assert.IsTrue(rab.Probability > 0.05);

            // Mann-Whitney test 1 vs. 2; they should be distinguished
            // with 1 consistently less than 2, so U abnormally small
            TestResult r12 = Sample.MannWhitneyTest(s1b, s2);
            Assert.IsTrue(r12.Probability < 0.05);

        }

        
        [TestMethod]
        public void AnovaTest () {

            Sample A = new Sample();
            A.Add(new double[] { 25, 30, 20, 32 });

            Sample B = new Sample();
            B.Add(new double[] { 30, 33, 29, 40, 36 });

            Sample C = new Sample();
            C.Add(new double[] { 32, 39, 35, 41, 44 });

            OneWayAnovaResult result = Sample.OneWayAnovaTest(A, B, C);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Factor.SumOfSquares + result.Residual.SumOfSquares, result.Total.SumOfSquares));
            Assert.IsTrue(result.Factor.DegreesOfFreedom + result.Residual.DegreesOfFreedom == result.Total.DegreesOfFreedom);

            Assert.IsTrue(result.Total.DegreesOfFreedom == A.Count + B.Count + C.Count - 1);

            Assert.IsTrue(result.Result.Statistic.Value == result.Factor.Result.Statistic.Value);

        }
        

        [TestMethod]
        public void AnovaDistribution () {

            ContinuousDistribution sDistribution = new NormalDistribution();
            Random rng = new Random(1);

            Sample fSample = new Sample();

            // do 100 ANOVAs
            for (int t = 0; t < 100; t++) {
                // each ANOVA has 4 groups
                List<Sample> groups = new List<Sample>();
                for (int g = 0; g < 4; g++) {
                    // each group has 3 data points
                    Sample group = new Sample();
                    for (int i = 0; i < 3; i++) {
                        group.Add(sDistribution.GetRandomValue(rng));
                    }
                    groups.Add(group);
                }

                OneWayAnovaResult result = Sample.OneWayAnovaTest(groups);
                fSample.Add(result.Factor.Result.Statistic);
            }

            // compare the distribution of F statistics to the expected distribution
            ContinuousDistribution fDistribution = new FisherDistribution(3, 8);
            TestResult kResult = fSample.KolmogorovSmirnovTest(fDistribution);
            Assert.IsTrue(kResult.Probability > 0.05);

        }

        [TestMethod]
        public void AnovaStudentAgreement () {

            // Create two samples.
            Sample A = new Sample(new double[] { 10, 20, 30 });
            Sample B = new Sample(new double[] { 15, 25, 35, 45 });

            // Do a Student t-test and a one-way ANOVA comparing them.
            TestResult ts = Sample.StudentTTest(A, B);
            TestResult ta = Sample.OneWayAnovaTest(A, B).Factor.Result;

            // The results should agree, with F = t^2 and same probability.
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ta.Statistic, MoreMath.Sqr(ts.Statistic)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ts.Probability, ta.Probability));
        }

        [TestMethod]
        public void ZTestDistribution () {
            Random rng = new Random(26);

            // define the sampling population (which must be normal for a z-test)
            ContinuousDistribution population = new NormalDistribution(2.0, 3.0);

            // Do a lot of z-tests.
            Sample zSample = new Sample();
            ContinuousDistribution zDistribution = null;
            for (int i = 0; i < 128; i++) {

                // We will do each z-test on a 4-point sample. We pick a small number so that any
                // non-normality for small sample sizes would show up.
                double[] sample = TestUtilities.CreateDataSample(rng, population, 4).ToArray();

                // For each sample, do a z-test against the population.
                TestResult zResult = sample.ZTest(population.Mean, population.StandardDeviation);
                zSample.Add(zResult.Statistic);
                zDistribution = zResult.Statistic.Distribution;
            }

            // The z's should be distributed normally, in accordance with the claimed distribution.

            TestResult result = zSample.KolmogorovSmirnovTest(zDistribution);
            Assert.IsTrue(result.Probability > 0.05);

        }

        [TestMethod]
        public void KruskalWallis () {

            // Generate four samples drawn from the same distribution,
            // and one sample drawn from a different distribution.
            Sample s1 = CreateSample(new ExponentialDistribution(1.0), 20, 1);
            Sample s2 = CreateSample(new ExponentialDistribution(1.0), 30, 2);
            Sample s3 = CreateSample(new ExponentialDistribution(1.0), 50, 3);
            Sample s4 = CreateSample(new ExponentialDistribution(1.0), 80, 4);
            Sample t4 = CreateSample(new ExponentialDistribution(2.0), 80, 4);

            // The test using the different sample should reject the null hypothesis.
            TestResult ssst = Sample.KruskalWallisTest(s1, s2, s3, t4);
            Assert.IsTrue(ssst.Probability < 0.05);

            // The test using the same samples should not reject the null hypothesis.
            TestResult ssss = Sample.KruskalWallisTest(s1, s2, s3, s4);
            Assert.IsTrue(ssss.Probability > 0.05);

        }

        [TestMethod]
        public void SampleTransform () {

            Sample S1 = new Sample(1.0, 2.0, 3.0, 4.0);

            Sample S2 = new Sample(1.0, 4.0, 9.0, 16.0);
            S2.Transform(Math.Sqrt);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(S1.Mean, S2.Mean));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(S1.Variance, S2.Variance));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(S1.Minimum, S2.Minimum));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(S1.Maximum, S2.Maximum));

        }

        [TestMethod]
        public void BivariateLogisticRegression () {

            double[] c = new double[] { -0.1, 1.0 };

            Random rng = new Random(1);
            UniformDistribution pointDistribution = new UniformDistribution(Interval.FromEndpoints(-4.0, 4.0));

            BivariateSample sample1 = new BivariateSample();
            MultivariateSample sample2 = new MultivariateSample(2);
            for (int k = 0; k < 1000; k++) {
                double x = pointDistribution.GetRandomValue(rng);
                double z = c[0] * x + c[1];
                double ez = Math.Exp(z);
                double p = ez / (1.0 + ez);
                double y = (rng.NextDouble() < p) ? 1.0 : 0.0;
                sample1.Add(x, y);
                sample2.Add(x, y);
            }

            Console.WriteLine(sample1.Covariance / sample1.X.Variance / sample1.Y.Mean / (1.0 - sample1.Y.Mean));
            Console.WriteLine(sample1.Covariance / sample1.X.Variance / sample1.Y.Variance);

            LinearLogisticRegressionResult result1 = sample1.LinearLogisticRegression();
            LinearLogisticRegressionResult result2 = sample2.TwoColumns(0, 1).LinearLogisticRegression();
            MultiLinearLogisticRegressionResult result3 = sample2.LogisticLinearRegression(1);

            //LinearLogisticRegressionResult result4 = Bivariate.LinearLogisticRegression(sample1.X.ToList(), sample1.Y.Select(x => (x == 1.0)).ToList());
            //MultiLinearLogisticRegressionResult result5 = sample1.Y.Select(x => (x == 1.0)).ToList().MultiLinearLogisticRegression(new IReadOnlyList<double>[] { sample1.X });
 
            // Fits should give same result
            /*
            for (int i = 0; i < result1.Dimension; i++) {
                Console.WriteLine("{0} {1} {2}", i, result1.Parameter(i), result3.Parameter(i) );
            }
            */
        }

        [TestMethod]
        public void MultivariateLinearRegression () {

            int outputIndex = 2;

            double[] c = new double[] { -1.0, 2.0, -3.0, 4.0 };

            Random rng = new Random(1001110000);
            UniformDistribution pointDistribution = new UniformDistribution(Interval.FromEndpoints(-4.0, 4.0));

            MultivariateSample sample = new MultivariateSample(c.Length);

            for (int k = 0; k < 1000; k++) {
                double[] row = new double[sample.Dimension];
                double z = 0.0;
                for (int i = 0; i < row.Length; i++) {
                    if (i == outputIndex) {
                        z += c[i];
                    } else {
                        row[i] = pointDistribution.GetRandomValue(rng);
                        z += row[i] * c[i];
                    }
                }
                double ez = Math.Exp(z);
                double p = ez / (1.0 + ez);
                row[outputIndex] = (rng.NextDouble() < p) ? 1.0 : 0.0;
                sample.Add(row);
            }

            MultiLinearLogisticRegressionResult result = sample.LogisticLinearRegression(outputIndex);

            for (int i = 0; i < result.Parameters.Count; i++) {
                Console.WriteLine(result.Parameters[i].Estimate);
                Assert.IsTrue(result.Parameters[i].Estimate.ConfidenceInterval(0.99).ClosedContains(c[i]));
            }
        }

        public void ValidateParameterCollection(ParameterCollection parameters) {
            for (int i = 0; i < parameters.Count; i++) {
                string iName = parameters[i].Name;
                Assert.IsTrue(parameters[i].Estimate.Value == parameters.ValuesVector[i]);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(parameters[i].Estimate.Uncertainty, Math.Sqrt(parameters.CovarianceMatrix[i, i])));
                Assert.IsTrue(parameters.VarianceOf(iName) == parameters.CovarianceMatrix[i, i]);
                foreach (Parameter jParameter in parameters) {
                    Assert.IsTrue(parameters.CovarianceOf(iName, jParameter.Name) == parameters.CovarianceMatrix[i, jParameter.Index]);
                }
            }
        }

    }
}
