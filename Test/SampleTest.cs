using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.Generic;
using System.Collections;

using Meta.Numerics;
using Meta.Numerics.Analysis;
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

            // set is 4, 6: M = 5, V = (1^2 + 1^2) / 2
            s.Remove(2.0);
            Assert.IsTrue(s.Count == 2);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(s.Mean, 5.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(s.Variance, 1.0));

        }

        private static Sample CreateSample (Distribution distribution, int count) {
            return (CreateSample(distribution, count, 1));
        }

        private static Sample CreateSample (Distribution distribution, int count, int seed) {

            Sample sample = new Sample();

            Random rng = new Random(seed);
            for (int i = 0; i < count; i++) {
                double x = distribution.GetRandomValue(rng);
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
        public void SampleMoments () {
            foreach (Distribution distribution in distributions) {

                Console.WriteLine(distribution.GetType().Name);
                
                Sample sample = CreateSample(distribution, 128);

                Assert.IsTrue(sample.Count == 128);

                UncertainValue m = sample.PopulationMean;
                Interval mi = m.ConfidenceInterval(0.95);
                Console.WriteLine("mean {0} {1}", mi, distribution.Mean);
                Assert.IsTrue(mi.ClosedContains(distribution.Mean));

                UncertainValue s = sample.PopulationStandardDeviation;
                Interval si = s.ConfidenceInterval(0.95);
                Console.WriteLine("stddev {0} {1}", si, distribution.StandardDeviation);
                Assert.IsTrue(si.ClosedContains(distribution.StandardDeviation));

                for (int n = 1; n < 8; n++) {
                    UncertainValue c = sample.PopulationMomentAboutMean(n);
                    Interval ci = c.ConfidenceInterval(0.95);
                    Console.WriteLine("C{0} {1} {2}", n, ci, distribution.MomentAboutMean(n));
                    Assert.IsTrue(ci.ClosedContains(distribution.MomentAboutMean(n)));

                    UncertainValue r = sample.PopulationMoment(n);
                    Interval ri = r.ConfidenceInterval(0.95);
                    Console.WriteLine("M{0} {1} {2}", n, ri, distribution.Moment(n));
                    Assert.IsTrue(ri.ClosedContains(distribution.Moment(n)));
                }
            }

        }


        [TestMethod]
        public void SamplePopulationMomentEstimateVariances () {

            Distribution d = new LognormalDistribution();

            // for various sample sizes...
            foreach (int n in TestUtilities.GenerateIntegerValues(4, 32, 8)) {

                Console.WriteLine("n={0}", n);

                // we are going to store values for a bunch of estimators and their uncertainties
                MultivariateSample estimates = new MultivariateSample("M1", "C2", "C3", "C4");
                MultivariateSample variances = new MultivariateSample("M1", "C2", "C3", "C4");

                // create a bunch of samples
                for (int i = 0; i < 256; i++) {

                    Sample s = TestUtilities.CreateSample(d, n, 512 * n + i + 1);

                    UncertainValue M1 = s.PopulationMean;
                    UncertainValue C2 = s.PopulationVariance;
                    UncertainValue C3 = s.PopulationMomentAboutMean(3);
                    UncertainValue C4 = s.PopulationMomentAboutMean(4);
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
            foreach (Distribution distribution in distributions) {
                
                Sample sample = CreateSample(distribution, 100);

                Interval iqr = sample.InterquartileRange;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(iqr.LeftEndpoint, sample.InverseLeftProbability(0.25)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(iqr.RightEndpoint, sample.InverseLeftProbability(0.75)));
            }

        }

        [TestMethod]
        public void NormalFit () {

            // pick mu >> sigma so that we get no negative values;
            // otherwise the attempt to fit to an exponential will fail
            Distribution distribution = new NormalDistribution(6.0, 2.0);
            Sample sample = CreateSample(distribution, 100);

            // fit to normal should be good
            FitResult nfit = NormalDistribution.FitToSample(sample);
            Console.WriteLine("P_n = {0}", nfit.GoodnessOfFit.LeftProbability);
            Assert.IsTrue(nfit.GoodnessOfFit.LeftProbability < 0.95);
            Assert.IsTrue(nfit.Parameter(0).ConfidenceInterval(0.95).ClosedContains(distribution.Mean));
            Assert.IsTrue(nfit.Parameter(1).ConfidenceInterval(0.95).ClosedContains(distribution.StandardDeviation));

            // fit to exponential should be bad
            FitResult efit = ExponentialDistribution.FitToSample(sample);
            Console.WriteLine("P_e = {0}", efit.GoodnessOfFit.LeftProbability);
            Assert.IsTrue(efit.GoodnessOfFit.LeftProbability > 0.95);

        }

        [TestMethod]
        public void NormalFitUncertainties () {

            NormalDistribution N = new NormalDistribution(-1.0, 2.0);

            // Create a bivariate sample to hold our fitted best mu and sigma values
            // so we can determine their covariance as well as their means and variances
            BivariateSample fits = new BivariateSample();

            double cmm = 0.0, css = 0.0, cms = 0.0;

            // A bunch of times, create a normal sample
            for (int i = 0; i < 64; i++) {

                // we will use small samples so the variation in mu and sigma will be more substantial
                Sample s = TestUtilities.CreateSample(N, 16, i);

                // fit each sample to a normal distribution
                FitResult fit = NormalDistribution.FitToSample(s);

                // and record the mu and sigma values from the fit into our bivariate sample
                fits.Add(fit.Parameter(0).Value, fit.Parameter(1).Value);

                // also record the claimed covariances among these parameters
                cmm += fit.Covariance(0, 0); css += fit.Covariance(1, 1); cms += fit.Covariance(0, 1);

            }

            cmm /= fits.Count; css /= fits.Count; cms /= fits.Count;


            // the mean fit values should agree with the population distribution
            Console.WriteLine("{0} {1}", fits.X.PopulationMean, N.Mean);
            Assert.IsTrue(fits.X.PopulationMean.ConfidenceInterval(0.95).ClosedContains(N.Mean));
            Console.WriteLine("{0} {1}", fits.Y.PopulationMean, N.StandardDeviation);
            Assert.IsTrue(fits.Y.PopulationMean.ConfidenceInterval(0.95).ClosedContains(N.StandardDeviation));

            // but also the covariances of those fit values should agree with the claimed covariances
            Console.WriteLine("{0} {1}", fits.X.PopulationVariance, cmm);
            Assert.IsTrue(fits.X.PopulationVariance.ConfidenceInterval(0.95).ClosedContains(cmm));
            Console.WriteLine("{0} {1}", fits.Y.PopulationVariance, css);
            Assert.IsTrue(fits.Y.PopulationVariance.ConfidenceInterval(0.95).ClosedContains(css));
            Console.WriteLine("{0} {1}", fits.PopulationCovariance, cms);
            Assert.IsTrue(fits.PopulationCovariance.ConfidenceInterval(0.95).ClosedContains(cms));

            /*
            Random rng = new Random(2718281);
            BivariateSample P = new BivariateSample();
            double cmm = 0.0;
            double css = 0.0;
            double cms = 0.0;
            for (int i = 0; i < 64; i++) {
                Sample s = new Sample();
                for (int j = 0; j < 16; j++) {
                    s.Add(N.GetRandomValue(rng));
                }
                FitResult r = NormalDistribution.FitToSample(s);
                P.Add(r.Parameter(0).Value, r.Parameter(1).Value);
                cmm += r.Covariance(0, 0);
                css += r.Covariance(1, 1);
                cms += r.Covariance(0, 1);
            }
            cmm /= P.Count;
            css /= P.Count;
            cms /= P.Count;

            Console.WriteLine("{0} {1}", P.X.PopulationMean, P.Y.PopulationMean);

            Assert.IsTrue(P.X.PopulationMean.ConfidenceInterval(0.95).ClosedContains(N.Mean));
            Assert.IsTrue(P.Y.PopulationMean.ConfidenceInterval(0.95).ClosedContains(N.StandardDeviation));

            Assert.IsTrue(P.X.PopulationVariance.ConfidenceInterval(0.95).ClosedContains(cmm));
            Assert.IsTrue(P.Y.PopulationVariance.ConfidenceInterval(0.95).ClosedContains(css));
            Assert.IsTrue(P.PopulationCovariance.ConfidenceInterval(0.95).ClosedContains(cms));
            */

        }

        [TestMethod]
        public void ExponentialFit () {

            ExponentialDistribution distribution = new ExponentialDistribution(5.0);
            Sample sample = CreateSample(distribution, 100);

            // fit to normal should be bad
            FitResult nfit = NormalDistribution.FitToSample(sample);
            Console.WriteLine("P_n = {0}", nfit.GoodnessOfFit.LeftProbability);
            Assert.IsTrue(nfit.GoodnessOfFit.LeftProbability > 0.95);

            // fit to exponential should be good
            FitResult efit = ExponentialDistribution.FitToSample(sample);
            Console.WriteLine("P_e = {0}", efit.GoodnessOfFit.LeftProbability);
            Assert.IsTrue(efit.GoodnessOfFit.LeftProbability < 0.95);
            Assert.IsTrue(efit.Parameter(0).ConfidenceInterval(0.95).ClosedContains(distribution.Mean));

        }

        [TestMethod]
        public void ExponentialFitUncertainty () {

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
                FitResult fit = ExponentialDistribution.FitToSample(sample);
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
        public void LognormalFit () {

            LognormalDistribution distribution = new LognormalDistribution(1.0, 2.0);
            Sample sample = CreateSample(distribution, 100);

            // fit to normal should be bad
            FitResult nfit = NormalDistribution.FitToSample(sample);
            Console.WriteLine("P_n = {0}", nfit.GoodnessOfFit.LeftProbability);
            Assert.IsTrue(nfit.GoodnessOfFit.LeftProbability > 0.95);

            // fit to lognormal should be good
            FitResult efit = LognormalDistribution.FitToSample(sample);
            Console.WriteLine("P_e = {0}", efit.GoodnessOfFit.LeftProbability);
            Assert.IsTrue(efit.GoodnessOfFit.LeftProbability < 0.95);
            Assert.IsTrue(efit.Parameter(0).ConfidenceInterval(0.95).ClosedContains(distribution.Mu));
            Assert.IsTrue(efit.Parameter(1).ConfidenceInterval(0.95).ClosedContains(distribution.Sigma));

        }

        [TestMethod]
        public void BetaFit () {
            // Create a very non-normal Beta distribution
            BetaDistribution B = new BetaDistribution(0.25, 3.0);

            // Sample from it
            Sample S = CreateSample(B, 128);

            // Fit the sample to a beta distribution
            // The fit should agree with the population and the fit should be good
            FitResult RB = BetaDistribution.FitToSample(S);
            Assert.IsTrue(RB.Parameter(0).ConfidenceInterval(0.99).ClosedContains(B.Alpha));
            Assert.IsTrue(RB.Parameter(1).ConfidenceInterval(0.99).ClosedContains(B.Beta));
            Assert.IsTrue(RB.GoodnessOfFit.LeftProbability < 0.99);

            // Fit to a normal distribution should be bad
            FitResult RN = NormalDistribution.FitToSample(S);
            Assert.IsTrue(RN.GoodnessOfFit.LeftProbability > 0.99);
        }

        [TestMethod]
        public void BetaFitUncertainty () {

            // check that the uncertainty in reported fit parameters is actually meaningful
            // it should be the standard deviation of fit parameter values in a sample of many fits

            // define a population distribution 
            Distribution distribution = new BetaDistribution(1.0 / 3.0, 2.0);

            // draw a lot of samples from it; fit each sample and
            // record the reported parameter value and error of each
            BivariateSample values = new BivariateSample();
            BivariateSample uncertainties = new BivariateSample();
            for (int i = 0; i < 50; i++) {
                Sample sample = CreateSample(distribution, 10, i);
                FitResult fit = BetaDistribution.FitToSample(sample);
                UncertainValue a = fit.Parameter(0);
                UncertainValue b = fit.Parameter(1);
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
            FitResult R = GammaDistribution.FitToSample(S);
            Assert.IsTrue(R.Parameter(0).ConfidenceInterval(0.95).ClosedContains(B.ShapeParameter));
            Assert.IsTrue(R.Parameter(1).ConfidenceInterval(0.95).ClosedContains(B.ScaleParameter));
            Assert.IsTrue(R.GoodnessOfFit.LeftProbability < 0.95);

        }

        [TestMethod]
        public void GammaFitUncertainty () {

            // check that the uncertainty in reported fit parameters is actually meaningful
            // it should be the standard deviation of fit parameter values in a sample of many fits

            // define a population distribution 
            Distribution distribution = new GammaDistribution(1.5, 2.0);

            // draw a lot of samples from it; fit each sample and
            // record the reported parameter value and error of each
            BivariateSample values = new BivariateSample();
            BivariateSample uncertainties = new BivariateSample();
            for (int i = 0; i < 100; i++) {
                Sample sample = CreateSample(distribution, 50, i);
                FitResult fit = GammaDistribution.FitToSample(sample);
                UncertainValue a = fit.Parameter(0);
                UncertainValue b = fit.Parameter(1);
                values.Add(a.Value, b.Value);
                uncertainties.Add(a.Uncertainty, b.Uncertainty);
            }

            // the reported errors should agree with the standard deviation of the reported parameters
            Assert.IsTrue(values.X.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(uncertainties.X.Mean));
            Assert.IsTrue(values.Y.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(uncertainties.Y.Mean));

        }

        [TestMethod]
        public void WaldFitUncertainties () {

            WaldDistribution wald = new WaldDistribution(3.5, 2.5);

            Random rng = new Random(314159);
            BivariateSample P = new BivariateSample();
            double cmm = 0.0;
            double css = 0.0;
            double cms = 0.0;
            for (int i = 0; i < 50; i++) {
                Sample s = new Sample();
                for (int j = 0; j < 50; j++) {
                    s.Add(wald.GetRandomValue(rng));
                }
                FitResult r = WaldDistribution.FitToSample(s);
                P.Add(r.Parameter(0).Value, r.Parameter(1).Value);
                cmm += r.Covariance(0, 0);
                css += r.Covariance(1, 1);
                cms += r.Covariance(0, 1);
            }
            cmm /= P.Count;
            css /= P.Count;
            cms /= P.Count;

            Console.WriteLine("{0} {1}", P.X.PopulationMean, P.Y.PopulationMean);

            Assert.IsTrue(P.X.PopulationMean.ConfidenceInterval(0.95).ClosedContains(wald.Mean));
            Assert.IsTrue(P.Y.PopulationMean.ConfidenceInterval(0.95).ClosedContains(wald.ShapeParameter));
            // the ML shape parameter estimate appears to be asymptoticly unbiased, as it must be according to ML fit theory,
            // but detectably upward biased for small n. we now correct for this.

            Console.WriteLine("{0} {1} {2}", P.X.PopulationVariance, P.Y.PopulationVariance, P.PopulationCovariance);
            Console.WriteLine("{0} {1} {2}", cmm, css, cms);

            Assert.IsTrue(P.X.PopulationVariance.ConfidenceInterval(0.95).ClosedContains(cmm));
            Assert.IsTrue(P.Y.PopulationVariance.ConfidenceInterval(0.95).ClosedContains(css));
            Assert.IsTrue(P.PopulationCovariance.ConfidenceInterval(0.95).ClosedContains(cms));

        }

        [TestMethod]
        public void SampleFitChiSquaredTest () {

            Distribution distribution = new ChiSquaredDistribution(4);
            Sample sample = CreateSample(distribution, 100);

            // fit to normal should be bad
            // this is harder than others, because a chi^2 isn't so very different from a normal; to help, increse N or decrease vu
            FitResult nfit = NormalDistribution.FitToSample(sample);
            Console.WriteLine("P_n = {0}", nfit.GoodnessOfFit.LeftProbability);
            Assert.IsTrue(nfit.GoodnessOfFit.LeftProbability > 0.95, String.Format("P_n = {0}", nfit.GoodnessOfFit.LeftProbability));

            // fit to exponential should also be bad
            FitResult efit = ExponentialDistribution.FitToSample(sample);
            Console.WriteLine("P_e = {0}", efit.GoodnessOfFit.LeftProbability);
            Assert.IsTrue(efit.GoodnessOfFit.LeftProbability > 0.95, String.Format("P_e = {0}", efit.GoodnessOfFit.LeftProbability));

        }

        [TestMethod]
        public void WeibullFit () {

            foreach (double alpha in new double[] { 0.03, 0.3, 1.3, 3.3 }) {
                WeibullDistribution W = new WeibullDistribution(2.2, alpha);
                Sample S = CreateSample(W, 50);
                FitResult R = WeibullDistribution.FitToSample(S);
                Assert.IsTrue(R.Parameter(0).ConfidenceInterval(0.95).ClosedContains(W.ScaleParameter));
                Assert.IsTrue(R.Parameter(1).ConfidenceInterval(0.95).ClosedContains(W.ShapeParameter));
            }

        }

        [TestMethod]
        public void WeibullFitUncertainties () {

            // check that the uncertainty in reported fit parameters is actually meaningful
            // it should be the standard deviation of fit parameter values in a sample of many fits

            // define a population distribution 
            Distribution distribution = new WeibullDistribution(2.5, 1.5);

            // draw a lot of samples from it; fit each sample and
            // record the reported parameter value and error of each
            BivariateSample values = new BivariateSample();
            MultivariateSample uncertainties = new MultivariateSample(3);
            for (int i = 0; i < 50; i++) {
                Sample sample = CreateSample(distribution, 10, i);
                FitResult fit = WeibullDistribution.FitToSample(sample);
                UncertainValue a = fit.Parameter(0);
                UncertainValue b = fit.Parameter(1);
                values.Add(a.Value, b.Value);
                uncertainties.Add(a.Uncertainty, b.Uncertainty, fit.Covariance(0,1));
            }

            // the reported errors should agree with the standard deviation of the reported parameters
            Assert.IsTrue(values.X.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(uncertainties.Column(0).Mean));
            Assert.IsTrue(values.Y.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(uncertainties.Column(1).Mean));
            //Console.WriteLine("{0} {1}", values.PopulationCovariance, uncertainties.Column(2).Mean);
            //Assert.IsTrue(values.PopulationCovariance.ConfidenceInterval(0.95).ClosedContains(uncertainties.Column(2).Mean));

        }

        [TestMethod]
        public void TTestDistribution () {

            // start with a normally distributed population
            Distribution xDistribution = new NormalDistribution(2.0, 3.0);
            Random rng = new Random(1);

            // draw 100 samples from it and compute the t statistic for each
            Sample tSample = new Sample();
            for (int i = 0; i < 100; i++) {

                // each sample has 9 values
                Sample xSample = new Sample();
                for (int j = 0; j < 9; j++) {
                    xSample.Add(xDistribution.GetRandomValue(rng));
                }
                //Sample xSample = CreateSample(xDistribution, 10, i);
                TestResult tResult = xSample.StudentTTest(2.0);
                double t = tResult.Statistic;
                Console.WriteLine("t = {0}", t);
                tSample.Add(t);
            }

            // sanity check our sample of t's
            Assert.IsTrue(tSample.Count == 100);

            // check that the t statistics are distributed as expected
            Distribution tDistribution = new StudentDistribution(9);

            // check on the mean
            Console.WriteLine("m = {0} vs. {1}", tSample.PopulationMean, tDistribution.Mean);
            Assert.IsTrue(tSample.PopulationMean.ConfidenceInterval(0.95).ClosedContains(tDistribution.Mean), String.Format("{0} vs. {1}", tSample.PopulationMean, tDistribution.Mean));

            // check on the standard deviation
            Console.WriteLine("s = {0} vs. {1}", tSample.PopulationStandardDeviation, tDistribution.StandardDeviation);
            Assert.IsTrue(tSample.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(tDistribution.StandardDeviation));

            // do a KS test
            TestResult ksResult = tSample.KolmogorovSmirnovTest(tDistribution);
            Assert.IsTrue(ksResult.LeftProbability < 0.95);
            Console.WriteLine("D = {0}", ksResult.Statistic);

            // check that we can distinguish the t distribution from a normal distribution?
        }

        [TestMethod]
        public void SignTestDistribution () {

            // start with a non-normally distributed population
            Distribution xDistribution = new ExponentialDistribution();
            Random rng = new Random(1);

            // draw 100 samples from it and compute the t statistic for each
            Sample wSample = new Sample();
            for (int i = 0; i < 100; i++) {

                // each sample has 8 observations
                Sample xSample = new Sample();
                for (int j = 0; j < 8; j++) { xSample.Add(xDistribution.GetRandomValue(rng)); }
                //Sample xSample = CreateSample(xDistribution, 8, i);
                TestResult wResult = xSample.SignTest(xDistribution.Median);
                double W = wResult.Statistic;
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
            Console.WriteLine("t = {0} vs. {1}", wSample.PopulationMomentAboutMean(3), wDistribution.MomentAboutMean(3));
            Assert.IsTrue(wSample.PopulationMomentAboutMean(3).ConfidenceInterval(0.95).ClosedContains(wDistribution.MomentAboutMean(3)));

            // check on the kuritosis
            Console.WriteLine("u = {0} vs. {1}", wSample.PopulationMomentAboutMean(4), wDistribution.MomentAboutMean(4));
            Assert.IsTrue(wSample.PopulationMomentAboutMean(4).ConfidenceInterval(0.95).ClosedContains(wDistribution.MomentAboutMean(4)));

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

                    TestResult result = Sample.KolmogorovSmirnovTest(aSamples[i], bSamples[j]);
                    Console.WriteLine("{0} v. {1}: D={2} P={3}", i, j, result.Statistic, result.LeftProbability);
                    if (i == j) {
                        Assert.IsTrue(result.LeftProbability < 0.90);
                    } else {
                        Assert.IsTrue(result.LeftProbability > 0.90);
                    }

                    // the order shouldn't matter
                    TestResult reverse = Sample.KolmogorovSmirnovTest(bSamples[j], aSamples[i]);
                    Assert.IsTrue(reverse.Statistic == result.Statistic);
                    Assert.IsTrue(reverse.RightProbability == result.RightProbability);

                }
            }

        }

        public void TestMoments (Distribution d) {

            // the support gives the limits of integration
            Interval support = d.Support;

            // raw moments
            double[] M = new double[6];
            for (int n = 0; n < 6; n++) {
                // define x^n p(x)
                Func<double, double> raw = delegate(double x) {
                    return (Math.Pow(x, n) * d.ProbabilityDensity(x));
                }; 
                // integrate it
                M[n] = FunctionMath.Integrate(raw, support);
                // compare with the claimed result
                Console.WriteLine("M{0} {1} v. {2}", n, M[n], d.Moment(n));
            }

            // central moments
            double[] C = new double[6];
            for (int n = 0; n < 6; n++) {
                // define (x-m)^n p(x)
                Func<double, double> central = delegate(double x) {
                    return (Math.Pow(x - M[1], n) * d.ProbabilityDensity(x));
                };
                // integrate it
                C[n] = FunctionMath.Integrate(central, support);
                // compare with the claimed result
                Console.WriteLine("C{0} {1} v. {2}", n, C[n], d.MomentAboutMean(n));
            }

            Console.WriteLine("Mean {0} v. {1}", M[1], d.Mean);
            Console.WriteLine("Standard Deviation {0} v. {1}", Math.Sqrt(C[2]), d.StandardDeviation);

        }

        [TestMethod]
        public void SampleKolmogorovSmirnovTest () {

            // this test has a whiff of meta-statistics about it
            // we want to make sure that the KS test statistic D is distributed according to the Kolmogorov
            // distribution; to do this, we create a sample of D statistics and do KS/Kuiper tests
            // comparing it to the claimed Kolmogorov distribution

            // start with any 'ol underlying distribution
            Distribution distribution = new UniformDistribution(Interval.FromEndpoints(-2.0, 4.0));

            // generate some samples from it, and for each one get a D statistic from a KS test
            Sample DSample = new Sample();
            Distribution DDistribution = null;
            for (int i = 0; i < 25; i++) {
                // the sample size must be large enough that the asymptotic assumptions are satistifed
                // at the moment this test fails if we make the sample size much smaller; we should
                // be able shrink this number when we expose the finite-sample distributions
                Sample sample = CreateSample(distribution, 250, i);
                TestResult ks = sample.KolmogorovSmirnovTest(distribution);
                double D = ks.Statistic;
                Console.WriteLine("D = {0}", D);
                DSample.Add(D);
                DDistribution = ks.Distribution;
            }

            // check on the mean
            Console.WriteLine("m = {0} vs. {1}", DSample.PopulationMean, DDistribution.Mean);
            Assert.IsTrue(DSample.PopulationMean.ConfidenceInterval(0.95).ClosedContains(DDistribution.Mean), String.Format("{0} vs. {1}", DSample.PopulationMean, DDistribution.Mean));

            // check on the standard deviation
            Console.WriteLine("s = {0} vs. {1}", DSample.PopulationStandardDeviation, DDistribution.StandardDeviation);
            Assert.IsTrue(DSample.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(DDistribution.StandardDeviation));

            // do a KS test comparing the sample to the expected distribution
            TestResult kst = DSample.KolmogorovSmirnovTest(DDistribution);
            Console.WriteLine("D = {0}, P = {1}", kst.Statistic, kst.LeftProbability);
            Assert.IsTrue(kst.LeftProbability < 0.95);

            // do a Kuiper test comparing the sample to the expected distribution
            TestResult kut = DSample.KuiperTest(DDistribution);
            Console.WriteLine("V = {0}, P = {1}", kut.Statistic, kut.LeftProbability);
            Assert.IsTrue(kut.LeftProbability < 0.95);


        }


        [TestMethod]
        public void SampleKuiperTest () {

            // this test has a whiff of meta-statistics about it
            // we want to make sure that the Kuiper test statistic V is distributed according to the Kuiper
            // distribution; to do this, we create a sample of V statistics and do KS/Kuiper tests
            // comparing it to the claimed Kuiper distribution
            
            // start with any 'ol underlying distribution
            Distribution distribution = new ExponentialDistribution(2.0);

            // generate some samples from it, and for each one get a V statistic from a KS test
            Sample VSample = new Sample();
            Distribution VDistribution = null;
            for (int i = 0; i < 25; i++) {
                // the sample size must be large enough that the asymptotic assumptions are satistifed
                // at the moment this test fails if we make the sample size much smaller; we should
                // be able shrink this number when we expose the finite-sample distributions
                Sample sample = CreateSample(distribution, 250, i);
                TestResult kuiper = sample.KuiperTest(distribution);
                double V = kuiper.Statistic;
                Console.WriteLine("V = {0}", V);
                VSample.Add(V);
                VDistribution = kuiper.Distribution;
            }

            // check on the mean
            Console.WriteLine("m = {0} vs. {1}", VSample.PopulationMean, VDistribution.Mean);
            Assert.IsTrue(VSample.PopulationMean.ConfidenceInterval(0.95).ClosedContains(VDistribution.Mean));

            // check on the standard deviation
            Console.WriteLine("s = {0} vs. {1}", VSample.PopulationStandardDeviation, VDistribution.StandardDeviation);
            Assert.IsTrue(VSample.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(VDistribution.StandardDeviation));

            // do a KS test comparing the sample to the expected distribution
            TestResult kst = VSample.KolmogorovSmirnovTest(VDistribution);
            Console.WriteLine("D = {0}, P = {1}", kst.Statistic, kst.LeftProbability);
            Assert.IsTrue(kst.LeftProbability < 0.95);

            // do a Kuiper test comparing the sample to the expected distribution
            TestResult kut = VSample.KuiperTest(VDistribution);
            Console.WriteLine("V = {0}, P = {1}", kut.Statistic, kut.LeftProbability);
            Assert.IsTrue(kut.LeftProbability < 0.95);

        }

        [TestMethod]
        public void SampleMaximumLikelihoodFit () {

            // normal distriubtion

            double mu = -1.0;
            double sigma = 2.0;
            Distribution nd = new NormalDistribution(mu, sigma);
            Sample ns = CreateSample(nd, 500);
            FitResult nr = ns.MaximumLikelihoodFit(new NormalDistribution(mu + 1.0, sigma + 1.0));

            Console.WriteLine(nr.Parameter(0));
            Console.WriteLine(nr.Parameter(1));

            Assert.IsTrue(nr.Dimension == 2);
            Assert.IsTrue(nr.Parameter(0).ConfidenceInterval(0.95).ClosedContains(mu));
            Assert.IsTrue(nr.Parameter(1).ConfidenceInterval(0.95).ClosedContains(sigma));

            Console.WriteLine(nr.Covariance(0,1));

            // test analytic expression
            Assert.IsTrue(TestUtilities.IsNearlyEqual(nr.Parameter(0).Value, ns.Mean, Math.Sqrt(TestUtilities.TargetPrecision)));

            // exponential distribution
            double em = 3.0;
            Distribution ed = new ExponentialDistribution(em);
            Sample es = CreateSample(ed, 100);
            FitResult er = es.MaximumLikelihoodFit(new ExponentialDistribution(em + 1.0));

            Console.WriteLine(er.Parameter(0));

            Assert.IsTrue(er.Dimension == 1);
            Assert.IsTrue(er.Parameter(0).ConfidenceInterval(0.95).ClosedContains(em));

            // test against analytic expression
            Assert.IsTrue(TestUtilities.IsNearlyEqual(er.Parameter(0).Value, es.Mean, Math.Sqrt(TestUtilities.TargetPrecision)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(er.Parameter(0).Uncertainty, es.Mean / Math.Sqrt(es.Count), Math.Sqrt(TestUtilities.TargetPrecision)));

            // lognormal distribution

            double l1 = -4.0;
            double l2 = 5.0;

            Distribution ld = new LognormalDistribution(l1, l2);
            Sample ls = CreateSample(ld, 100);
            FitResult lr = ls.MaximumLikelihoodFit(new LognormalDistribution(l1 + 1.0, l2 + 1.0));

            Console.WriteLine(lr.Parameter(0));
            Console.WriteLine(lr.Parameter(1));
            Console.WriteLine(lr.Covariance(0, 1));

            // weibull distribution

            double w_scale = 4.0;
            double w_shape = 2.0;
            Distribution w_d = new WeibullDistribution(w_scale, w_shape);
            Sample w_s = CreateSample(w_d, 20);
            FitResult w_r = w_s.MaximumLikelihoodFit(new WeibullDistribution(1.0, 0.5));

            Console.WriteLine(w_r.Parameter(0));
            Console.WriteLine(w_r.Parameter(1));
            Console.WriteLine(w_r.Covariance(0, 1));

            Assert.IsTrue(w_r.Parameter(0).ConfidenceInterval(0.95).ClosedContains(w_scale));
            Assert.IsTrue(w_r.Parameter(1).ConfidenceInterval(0.95).ClosedContains(w_shape));

            // logistic distribution
            double logistic_m = -3.0;
            double logistic_s = 2.0;
            Distribution logistic_distribution = new LogisticDistribution(logistic_m, logistic_s);
            Sample logistic_sample = CreateSample(logistic_distribution, 100);
            FitResult logistic_result = logistic_sample.MaximumLikelihoodFit(new LogisticDistribution());

            Console.WriteLine("Logistic:");
            Console.WriteLine(logistic_result.Parameter(0));
            Console.WriteLine(logistic_result.Parameter(1));

            Assert.IsTrue(logistic_result.Dimension == 2);
            Assert.IsTrue(logistic_result.Parameter(0).ConfidenceInterval(0.95).ClosedContains(logistic_m));
            Assert.IsTrue(logistic_result.Parameter(1).ConfidenceInterval(0.95).ClosedContains(logistic_s));


            // beta distribution
            // not yet!
            /*
            double beta_alpha = 0.5;
            double beta_beta = 2.0;
            Distribution beta_distribution = new BetaDistribution(beta_alpha, beta_beta);
            Sample beta_sample = CreateSample(beta_distribution, 100);
            FitResult beta_result = beta_sample.MaximumLikelihoodFit(new BetaDistribution(1.0, 1.0));

            Console.WriteLine("Beta:");
            Console.WriteLine(beta_result.Parameter(0));
            Console.WriteLine(beta_result.Parameter(1));

            Assert.IsTrue(beta_result.Dimension == 2);
            Assert.IsTrue(beta_result.Parameter(0).ConfidenceInterval(0.95).ClosedContains(beta_alpha));
            Assert.IsTrue(beta_result.Parameter(1).ConfidenceInterval(0.95).ClosedContains(beta_beta));
            */
        }

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


            // create 3 samles
            // 1 and 2 have the same mean but different variances, the F test should catch the difference
            // 1 and 3 have different means but the same variance, the F test should rule them equivilent
            Sample sample1 = CreateSample(new NormalDistribution(1.0, 1.0), 20, 1);
            Sample sample2 = CreateSample(new NormalDistribution(1.0, 2.0), 20, 2);
            Sample sample3 = CreateSample(new NormalDistribution(3.0, 1.0), 20, 3);

            TestResult f12 = Sample.FisherFTest(sample1, sample2);
            TestResult f21 = Sample.FisherFTest(sample2, sample1);

            // sample 1 has a smaller variance
            Console.WriteLine(f12.Statistic);
            Assert.IsTrue(f12.Statistic < 1.0);

            // 1/2 is the inverse of 2/1
            Assert.IsTrue(TestUtilities.IsNearlyEqual(f12.Statistic, 1.0 / f21.Statistic));
            
            // the F test detects the difference between the variance of 1 and 2
            Console.WriteLine(f12.LeftProbability);
            Assert.IsTrue(f12.RightProbability > 0.95);

            // the F test detects no difference between the variance of 1 and 3
            TestResult f13 = Sample.FisherFTest(sample1, sample3);
            Console.WriteLine(f13.Statistic);
            Console.WriteLine(f13.LeftProbability);
            Assert.IsTrue((f13.LeftProbability > 0.05) && (f13.RightProbability > 0.05));

        }

        [TestMethod]
        public void SampleMannWhitneyTest () {

            // define two non-normal distributions
            Distribution d1 = new ExponentialDistribution(2.0);
            Distribution d2 = new ExponentialDistribution(3.0);

            // create three samples from them
            Sample s1a = CreateSample(d1, 20, 1);
            Sample s1b = CreateSample(d1, 30, 2);
            Sample s2 = CreateSample(d2, 40, 3);

            // Mann-Whitney test 1a vs. 1b; they should not be distinguished
            TestResult rab = Sample.MannWhitneyTest(s1a, s1b);
            Console.WriteLine("{0} {1}", rab.Statistic, rab.LeftProbability);
            Assert.IsTrue((rab.LeftProbability < 0.95) && (rab.RightProbability < 0.95));

            // Mann-Whitney test 1 vs. 2; they should be distinguished
            // with 1 consistently less than 2, so U abnormally small
            TestResult r12 = Sample.MannWhitneyTest(s1b, s2);
            Console.WriteLine("{0} {1}", r12.Statistic, r12.LeftProbability);
            Assert.IsTrue(r12.RightProbability > 0.95);

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

            Assert.IsTrue(result.Result.Statistic == result.Factor.Result.Statistic);

        }
        

        [TestMethod]
        public void AnovaDistribution () {

            Distribution sDistribution = new NormalDistribution();
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
            Distribution fDistribution = new FisherDistribution(3, 8);
            Console.WriteLine("m={0} s={1}", fSample.PopulationMean, fSample.PopulationStandardDeviation);
            TestResult kResult = fSample.KolmogorovSmirnovTest(fDistribution);
            Console.WriteLine(kResult.LeftProbability);
            Assert.IsTrue(kResult.LeftProbability < 0.95);

        }

        [TestMethod]
        public void AnovaStudentAgreement () {

            // create two random samples
            Sample A = new Sample(new double[] { 10, 20, 30 });
            Sample B = new Sample(new double[] { 15, 25, 35, 45 });
            Console.WriteLine("A m={0} v={1}", A.Mean, A.Variance);
            Console.WriteLine("B m={0} v={1}", B.Mean, B.Variance);

            // do a Student t-test and a one-way ANOVA
            //TestResult ts = Sample.StudentTTest(A, B);
            TestResult ts = Sample.StudentTTest(A, B);
            TestResult ta = Sample.OneWayAnovaTest(A, B).Factor.Result;


            // the results should agree, with F = t^2 and same probability
            Console.WriteLine("t F {0} {1}", MoreMath.Pow(ts.Statistic, 2), ta.Statistic);
            Console.WriteLine("P P {0} {1}", ts.LeftProbability, ta.LeftProbability);

            StudentDistribution ds = (StudentDistribution)ts.Distribution;
            FisherDistribution da = (FisherDistribution)ta.Distribution;
            Console.WriteLine(ds.DegreesOfFreedom);
            Console.WriteLine("{0} {1}", da.NumeratorDegreesOfFreedom, da.DenominatorDegreesOfFreedom);

        }

        [TestMethod]
        public void ZTestDistribution () {

            Random rng = new Random(1);

            // define the sampling population (which must be normal for a z-test)
            Distribution population = new NormalDistribution(2.0, 3.0);

            // collect 100 samples
            Sample zSample = new Sample();
            for (int i = 0; i < 100; i++) {

                // each z-statistic is formed by making a 4-count sample from a normal distribution
                Sample sample = new Sample();
                for (int j = 0; j < 4; j++) {
                    sample.Add(population.GetRandomValue(rng));
                }

                // for each sample, do a z-test against the population
                TestResult zResult = sample.ZTest(population.Mean, population.StandardDeviation);
                zSample.Add(zResult.Statistic);

            }

            // the z's should be distrubuted normally

            TestResult result = zSample.KolmogorovSmirnovTest(new NormalDistribution());
            Console.WriteLine("{0} {1}", result.Statistic, result.LeftProbability);
            Assert.IsTrue((result.LeftProbability > 0.05) && (result.LeftProbability < 0.95));

        }

        [TestMethod]
        public void KruskalWallis () {

            // generate four samples drawn from the same distribution,
            // and a replacement sample for the last one drawn from a different distribution
            Sample s1 = CreateSample(new ExponentialDistribution(1.0), 20, 1);
            Sample s2 = CreateSample(new ExponentialDistribution(1.0), 30, 2);
            Sample s3 = CreateSample(new ExponentialDistribution(1.0), 50, 3);
            Sample s4 = CreateSample(new ExponentialDistribution(1.0), 80, 4);
            Sample t4 = CreateSample(new ExponentialDistribution(2.0), 80, 4);

            // the test using the different sample should reject the null hypothesis
            TestResult ssst = Sample.KruskalWallisTest(s1, s2, s3, t4);
            Console.WriteLine(ssst.RightProbability);
            Assert.IsTrue(ssst.RightProbability < 0.05);

            // the test using the same samples should not reject the null hypothesis
            TestResult ssss = Sample.KruskalWallisTest(s1, s2, s3, s4);
            Console.WriteLine(ssss.RightProbability);
            Assert.IsTrue(ssss.RightProbability > 0.05);

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

    }
}
