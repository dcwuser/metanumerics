using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Data;
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
            NormalFitResult nfit = NormalDistribution.FitToSample(sample);
            Assert.IsTrue(nfit.GoodnessOfFit.Probability > 0.05);
            Assert.IsTrue(nfit.Mean.ConfidenceInterval(0.95).ClosedContains(distribution.Mean));
            Assert.IsTrue(nfit.StandardDeviation.ConfidenceInterval(0.95).ClosedContains(distribution.StandardDeviation));

            // fit to exponential should be bad
            ExponentialFitResult efit = ExponentialDistribution.FitToSample(sample);
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
                NormalFitResult fit = NormalDistribution.FitToSample(s);

                // and record the mu and sigma values from the fit into our bivariate sample
                parameters.Add(fit.Mean.Value, fit.StandardDeviation.Value);

                // also record the claimed covariances among these parameters
                covariances.Add(fit.Parameters.CovarianceMatrix[0, 0], fit.Parameters.CovarianceMatrix[1, 1], fit.Parameters.CovarianceMatrix[0, 1]);
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
            NormalFitResult nfit = NormalDistribution.FitToSample(sample);
            Assert.IsTrue(nfit.GoodnessOfFit.Probability < 0.05);

            // fit to exponential should be good
            ExponentialFitResult efit = ExponentialDistribution.FitToSample(sample);
            Assert.IsTrue(efit.GoodnessOfFit.Probability > 0.05);
            Assert.IsTrue(efit.Mean.ConfidenceInterval(0.95).ClosedContains(distribution.Mean));

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
                ExponentialFitResult fit = ExponentialDistribution.FitToSample(sample);
                UncertainValue lambda = fit.Parameters[0].Estimate;
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

            // Do a bunch of fits, record reported parameters and variances
            for (int i = 0; i < 32; i++) {
                Sample s = SampleTest.CreateSample(d, 64, i);

                GumbelFitResult r = GumbelDistribution.FitToSample(s);
                parameters.Add(r.Location.Value, r.Scale.Value);
                variances.Add(r.Parameters.CovarianceMatrix[0, 0], r.Parameters.CovarianceMatrix[1, 1], r.Parameters.CovarianceMatrix[0, 1]);

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

            FrameTable results = new FrameTable();
            results.AddColumns<double>("Mean", "Shape", "MeanVariance", "ShapeVariance", "MeanShapeCovariance");

            for (int i = 0; i < 128; i++) {
                Sample sample = SampleTest.CreateSample(wald, 16, i);

                WaldFitResult result = WaldDistribution.FitToSample(sample);
                Assert.IsTrue(result.Mean.Value == result.Parameters.ValuesVector[result.Parameters.IndexOf("Mean")]);
                Assert.IsTrue(result.Shape.Value == result.Parameters.ValuesVector[result.Parameters.IndexOf("Shape")]);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Parameters.VarianceOf("Mean"), MoreMath.Sqr(result.Mean.Uncertainty)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Parameters.VarianceOf("Shape"), MoreMath.Sqr(result.Shape.Uncertainty)));
                results.AddRow(
                    result.Mean.Value, result.Shape.Value,
                    result.Parameters.VarianceOf("Mean"), result.Parameters.VarianceOf("Shape"), result.Parameters.CovarianceOf("Mean", "Shape")
                );
            }

            Assert.IsTrue(results["Mean"].As<double>().PopulationMean().ConfidenceInterval(0.99).ClosedContains(wald.Mean));
            Assert.IsTrue(results["Shape"].As<double>().PopulationMean().ConfidenceInterval(0.99).ClosedContains(wald.Shape));

            Assert.IsTrue(results["Mean"].As<double>().PopulationVariance().ConfidenceInterval(0.99).ClosedContains(results["MeanVariance"].As<double>().Median()));
            Assert.IsTrue(results["Shape"].As<double>().PopulationVariance().ConfidenceInterval(0.99).ClosedContains(results["ShapeVariance"].As<double>().Median()));
            Assert.IsTrue(results["Mean"].As<double>().PopulationCovariance(results["Shape"].As<double>()).ConfidenceInterval(0.99).ClosedContains(results["MeanShapeCovariance"].As<double>().Median()));
        }

        [TestMethod]
        public void RayleighFit () {

            RayleighDistribution rayleigh = new RayleighDistribution(3.2);

            Sample parameter = new Sample();
            Sample variance = new Sample();

            for (int i = 0; i < 128; i++) {
                // We pick a quite-small sample, because we have a finite-n unbiased estimator.
                Sample s = SampleTest.CreateSample(rayleigh, 8, i);

                RayleighFitResult r = RayleighDistribution.FitToSample(s);
                parameter.Add(r.Scale.Value);
                variance.Add(r.Parameters.VarianceOf("Scale"));

                Assert.IsTrue(r.GoodnessOfFit.Probability > 0.01);
            }

            Assert.IsTrue(parameter.PopulationMean.ConfidenceInterval(0.99).ClosedContains(rayleigh.Scale));
            Assert.IsTrue(parameter.PopulationVariance.ConfidenceInterval(0.99).ClosedContains(variance.Median));

        }

    }
}
