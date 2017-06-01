using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace Test {

    [TestClass]
    public class SampleFitTest {

        [TestMethod]
        public void NormalFit () {

            // pick mu >> sigma so that we get no negative values;
            // otherwise the attempt to fit to an exponential will fail
            ContinuousDistribution distribution = new NormalDistribution(6.0, 2.0);
            Sample sample = SampleTest.CreateSample(distribution, 100);

            // fit to normal should be good
            FitResult nfit = NormalDistribution.FitToSample(sample);
            Assert.IsTrue(nfit.GoodnessOfFit.Probability > 0.05);
            Assert.IsTrue(nfit.Parameter(0).ConfidenceInterval(0.95).ClosedContains(distribution.Mean));
            Assert.IsTrue(nfit.Parameter(1).ConfidenceInterval(0.95).ClosedContains(distribution.StandardDeviation));

            // fit to exponential should be bad
            FitResult efit = ExponentialDistribution.FitToSample(sample);
            Assert.IsTrue(efit.GoodnessOfFit.Probability < 0.05);
        }

        [TestMethod]
        public void NormalFitCovariances () {

            NormalDistribution N = new NormalDistribution(-1.0, 2.0);

            // Create a bivariate sample to hold our fitted best mu and sigma values
            // so we can determine their covariance as well as their means and variances
            BivariateSample parameters = new BivariateSample();
            MultivariateSample covariances = new MultivariateSample(3); 

            // A bunch of times, create a normal sample
            for (int i = 0; i < 128; i++) {

                // We use small samples so the variation in mu and sigma will be more substantial.
                Sample s = TestUtilities.CreateSample(N, 8, i);

                // Fit each sample to a normal distribution
                FitResult fit = NormalDistribution.FitToSample(s);

                // and record the mu and sigma values from the fit into our bivariate sample
                parameters.Add(fit.Parameter(0).Value, fit.Parameter(1).Value);

                // also record the claimed covariances among these parameters
                covariances.Add(fit.Covariance(0, 0), fit.Covariance(1, 1), fit.Covariance(0, 1));
            }

            // the mean fit values should agree with the population distribution
            Assert.IsTrue(parameters.X.PopulationMean.ConfidenceInterval(0.95).ClosedContains(N.Mean));
            Assert.IsTrue(parameters.Y.PopulationMean.ConfidenceInterval(0.95).ClosedContains(N.StandardDeviation));

            // but also the covariances of those fit values should agree with the claimed covariances
            Assert.IsTrue(parameters.X.PopulationVariance.ConfidenceInterval(0.95).ClosedContains(covariances.Column(0).Mean));
            Assert.IsTrue(parameters.Y.PopulationVariance.ConfidenceInterval(0.95).ClosedContains(covariances.Column(1).Mean));
            Assert.IsTrue(parameters.PopulationCovariance.ConfidenceInterval(0.95).ClosedContains(covariances.Column(2).Mean));

        }

        [TestMethod]
        public void ExponentialFit () {

            ExponentialDistribution distribution = new ExponentialDistribution(5.0);
            Sample sample = SampleTest.CreateSample(distribution, 100);

            // fit to normal should be bad
            FitResult nfit = NormalDistribution.FitToSample(sample);
            Assert.IsTrue(nfit.GoodnessOfFit.Probability < 0.05);

            // fit to exponential should be good
            FitResult efit = ExponentialDistribution.FitToSample(sample);
            Assert.IsTrue(efit.GoodnessOfFit.Probability > 0.05);
            Assert.IsTrue(efit.Parameter(0).ConfidenceInterval(0.95).ClosedContains(distribution.Mean));

        }

        [TestMethod]
        public void ExponentialFitUncertainty () {

            // check that the uncertainty in reported fit parameters is actually meaningful
            // it should be the standard deviation of fit parameter values in a sample of many fits

            // define a population distribution 
            ExponentialDistribution distribution = new ExponentialDistribution(4.0);

            // draw a lot of samples from it; fit each sample and
            // record the reported parameter value and error of each
            Sample values = new Sample();
            Sample uncertainties = new Sample();
            for (int i = 0; i < 128; i++) {
                Sample sample = SampleTest.CreateSample(distribution, 8, i);
                FitResult fit = ExponentialDistribution.FitToSample(sample);
                UncertainValue lambda = fit.Parameter(0);
                values.Add(lambda.Value);
                uncertainties.Add(lambda.Uncertainty);
            }

            // the reported values should agree with the source distribution
            Assert.IsTrue(values.PopulationMean.ConfidenceInterval(0.95).ClosedContains(distribution.Mean));

            // the reported errors should agree with the standard deviation of the reported parameters
            Assert.IsTrue(values.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(uncertainties.Mean));

        }

        [TestMethod]
        public void GumbelFit () {

            GumbelDistribution d = new GumbelDistribution(-1.0, 2.0);

            MultivariateSample parameters = new MultivariateSample(2);
            MultivariateSample variances = new MultivariateSample(3);

            // Do a bunch of fits, record reported parameters an variances
            for (int i = 0; i < 32; i++) {
                Sample s = SampleTest.CreateSample(d, 64, i);

                FitResult r = GumbelDistribution.FitToSample(s);
                parameters.Add(r.Parameters);
                variances.Add(r.Covariance(0, 0), r.Covariance(1, 1), r.Covariance(0, 1));

                Assert.IsTrue(r.GoodnessOfFit.Probability > 0.01);
            }

            // The reported parameters should agree with the underlying parameters
            Assert.IsTrue(parameters.Column(0).PopulationMean.ConfidenceInterval(0.99).ClosedContains(d.Location));
            Assert.IsTrue(parameters.Column(1).PopulationMean.ConfidenceInterval(0.99).ClosedContains(d.Scale));

            // The reported covariances should agree with the observed covariances
            Assert.IsTrue(parameters.Column(0).PopulationVariance.ConfidenceInterval(0.99).ClosedContains(variances.Column(0).Mean));
            Assert.IsTrue(parameters.Column(1).PopulationVariance.ConfidenceInterval(0.99).ClosedContains(variances.Column(1).Mean));
            Assert.IsTrue(parameters.TwoColumns(0, 1).PopulationCovariance.ConfidenceInterval(0.99).ClosedContains(variances.Column(2).Mean));

        }

        [TestMethod]
        public void WaldFit () {

            WaldDistribution wald = new WaldDistribution(3.5, 2.5);

            BivariateSample parameters = new BivariateSample();
            MultivariateSample variances = new MultivariateSample(3);

            for (int i = 0; i < 128; i++) {
                Sample s = SampleTest.CreateSample(wald, 16, i);

                FitResult r = WaldDistribution.FitToSample(s);
                parameters.Add(r.Parameters[0], r.Parameters[1]);
                variances.Add(r.Covariance(0, 0), r.Covariance(1, 1), r.Covariance(0, 1));

                Assert.IsTrue(r.GoodnessOfFit.Probability > 0.01);
            }

            Assert.IsTrue(parameters.X.PopulationMean.ConfidenceInterval(0.99).ClosedContains(wald.Mean));
            Assert.IsTrue(parameters.Y.PopulationMean.ConfidenceInterval(0.99).ClosedContains(wald.Shape));

            Assert.IsTrue(parameters.X.PopulationVariance.ConfidenceInterval(0.99).ClosedContains(variances.Column(0).Median));
            Assert.IsTrue(parameters.Y.PopulationVariance.ConfidenceInterval(0.99).ClosedContains(variances.Column(1).Median));
            Assert.IsTrue(parameters.PopulationCovariance.ConfidenceInterval(0.99).ClosedContains(variances.Column(2).Median));

        }

        [TestMethod]
        public void RayleighFit () {

            RayleighDistribution rayleigh = new RayleighDistribution(3.2);

            Sample parameter = new Sample();
            Sample variance = new Sample();

            for (int i = 0; i < 128; i++) {
                // We pick a quite-small sample, because we have a finite-n unbiased estimator.
                Sample s = SampleTest.CreateSample(rayleigh, 8, i);

                FitResult r = RayleighDistribution.FitToSample(s);
                parameter.Add(r.Parameters[0]);
                variance.Add(r.Covariance(0, 0));

                Assert.IsTrue(r.GoodnessOfFit.Probability > 0.01);
            }

            Assert.IsTrue(parameter.PopulationMean.ConfidenceInterval(0.99).ClosedContains(rayleigh.Scale));
            Assert.IsTrue(parameter.PopulationVariance.ConfidenceInterval(0.99).ClosedContains(variance.Median));

        }

    }
}
