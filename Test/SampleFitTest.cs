using System;
using System.Linq;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Data;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;
using System.Collections.Generic;

namespace Test {

    [TestClass]
    public class SampleFitTest {

        private IEnumerable<double> CreateSample(ContinuousDistribution distribution, int count) {
            return CreateSample(distribution, count, new Random(1));
        }

        private IEnumerable<double> CreateSample (ContinuousDistribution distribution, int count, Random rng) {
            for (int i = 0; i < count; i++) {
                yield return distribution.GetRandomValue(rng);
            }
        }

        [TestMethod]
        public void NormalFit () {

            // pick mu >> sigma so that we get no negative values;
            // otherwise the attempt to fit to an exponential will fail
            ContinuousDistribution distribution = new NormalDistribution(6.0, 2.0);
            List<double> sample = distribution.GetRandomValues(new Random(1), 100).ToList();

            // fit to normal should be good
            NormalFitResult nfit = NormalDistribution.FitToSample(sample);
            Assert.IsTrue(nfit.GoodnessOfFit.Probability > 0.05);
            Assert.IsTrue(nfit.Mean.ConfidenceInterval(0.95).Contains(distribution.Mean));
            Assert.IsTrue(nfit.StandardDeviation.ConfidenceInterval(0.95).Contains(distribution.StandardDeviation));

            // fit to exponential should be bad
            ExponentialFitResult efit = ExponentialDistribution.FitToSample(sample);
            Assert.IsTrue(efit.GoodnessOfFit.Probability < 0.05);
        }

        [TestMethod]
        public void NormalFitCovariances () {

            NormalDistribution N = new NormalDistribution(-1.0, 2.0);

            // Create a bivariate sample to hold our fitted best mu and sigma values
            // so we can determine their covariance as well as their means and variances
            FrameTable results = new FrameTable();
            results.AddColumns<double>("m", "s", "cmm", "css", "cms");

            // A bunch of times, create a normal sample
            Random rng = new Random(1);
            for (int i = 0; i < 128; i++) {

                // We use small samples so the variation in mu and sigma will be more substantial.
                List<double> s = N.GetRandomValues(rng, 8).ToList();

                // Fit each sample to a normal distribution
                NormalFitResult fit = NormalDistribution.FitToSample(s);

                results.AddRow(
                    fit.Mean.Value,
                    fit.StandardDeviation.Value,
                    fit.Parameters.CovarianceOf("Mean", "Mean"),
                    fit.Parameters.CovarianceOf("StandardDeviation", "StandardDeviation"),
                    fit.Parameters.CovarianceOf("Mean", "StandardDeviation")
                );
            }

            // The mean fit values should agree with the population distribution
            Assert.IsTrue(results["m"].As<double>().PopulationMean().ConfidenceInterval(0.95).Contains(N.Mean));
            Assert.IsTrue(results["s"].As<double>().PopulationMean().ConfidenceInterval(0.95).Contains(N.StandardDeviation));

            // But also the covariances of those fit values should agree with the claimed covariances
            Assert.IsTrue(results["m"].As<double>().PopulationVariance().ConfidenceInterval(0.95).Contains(results["cmm"].As<double>().Mean()));
            Assert.IsTrue(results["s"].As<double>().PopulationVariance().ConfidenceInterval(0.95).Contains(results["css"].As<double>().Mean()));
            Assert.IsTrue(Bivariate.PopulationCovariance(results["m"].As<double>(), results["s"].As<double>()).ConfidenceInterval(0.95).Contains(results["cms"].As<double>().Mean()));

        }

        [TestMethod]
        public void ExponentialFit () {

            ExponentialDistribution distribution = new ExponentialDistribution(5.0);
            List<double> sample = distribution.GetRandomValues(new Random(1), 100).ToList();

            // fit to normal should be bad
            NormalFitResult nfit = NormalDistribution.FitToSample(sample);
            Assert.IsTrue(nfit.GoodnessOfFit.Probability < 0.05);

            // fit to exponential should be good
            ExponentialFitResult efit = ExponentialDistribution.FitToSample(sample);
            Assert.IsTrue(efit.GoodnessOfFit.Probability > 0.05);
            Assert.IsTrue(efit.Mean.ConfidenceInterval(0.95).Contains(distribution.Mean));

        }

        [TestMethod]
        public void ExponentialFitUncertainty () {

            // check that the uncertainty in reported fit parameters is actually meaningful
            // it should be the standard deviation of fit parameter values in a sample of many fits

            // define a population distribution 
            ExponentialDistribution distribution = new ExponentialDistribution(4.0);

            // draw a lot of samples from it; fit each sample and
            // record the reported parameter value and error of each
            List<double> values = new List<double>();
            List<double> uncertainties = new List<double>();
            Random rng = new Random(1);
            for (int i = 0; i < 128; i++) {
                List<double> sample = distribution.GetRandomValues(rng, 8).ToList();
                ExponentialFitResult fit = ExponentialDistribution.FitToSample(sample);
                UncertainValue lambda = fit.Parameters[0].Estimate;
                values.Add(lambda.Value);
                uncertainties.Add(lambda.Uncertainty);
            }

            // the reported values should agree with the source distribution
            Assert.IsTrue(values.PopulationMean().ConfidenceInterval(0.95).ClosedContains(distribution.Mean));

            // the reported errors should agree with the standard deviation of the reported parameters
            Assert.IsTrue(values.PopulationStandardDeviation().ConfidenceInterval(0.95).ClosedContains(uncertainties.Mean()));

        }

        [TestMethod]
        public void GumbelFit () {

            GumbelDistribution gumbel = new GumbelDistribution(-1.0, 2.0);

            Random rng = new Random(22);

            FrameTable results = new FrameTable();
            results.AddColumns<double>("Location", "Scale", "LocationVariance", "ScaleVariance", "LocationScaleCovariance");

            // Do a bunch of fits, record reported parameters and variances
            for (int i = 0; i < 32; i++) {

                List<double> s = gumbel.GetRandomValues(rng, 64).ToList();

                GumbelFitResult r = GumbelDistribution.FitToSample(s);
                results.AddRow(
                    r.Parameters["Location"].Estimate.Value, r.Parameters["Scale"].Estimate.Value,
                    r.Parameters.CovarianceOf("Location", "Location"), r.Parameters.CovarianceOf("Scale", "Scale"), r.Parameters.CovarianceOf("Location", "Scale")
                );

                Assert.IsTrue(r.GoodnessOfFit.Probability > 0.01);
            }

            // The reported parameters should agree with the underlying parameters
            Assert.IsTrue(results["Location"].As<double>().PopulationMean().ConfidenceInterval(0.99).ClosedContains(gumbel.Location));
            Assert.IsTrue(results["Scale"].As<double>().PopulationMean().ConfidenceInterval(0.99).ClosedContains(gumbel.Scale));

            // The reported covariances should agree with the observed covariances
            Assert.IsTrue(results["Location"].As<double>().PopulationVariance().ConfidenceInterval(0.99).Contains(results["LocationVariance"].As<double>().Mean()));
            Assert.IsTrue(results["Scale"].As<double>().PopulationVariance().ConfidenceInterval(0.99).Contains(results["ScaleVariance"].As<double>().Mean()));
            Assert.IsTrue(Bivariate.PopulationCovariance(results["Location"].As<double>(), results["Scale"].As<double>()).ConfidenceInterval(0.99).Contains(results["LocationScaleCovariance"].As<double>().Mean()));

        }

        [TestMethod]
        public void WaldFit () {

            WaldDistribution wald = new WaldDistribution(3.5, 2.5);
            Random rng = new Random(33);

            FrameTable results = new FrameTable();
            results.AddColumns<double>("Mean", "Shape", "MeanVariance", "ShapeVariance", "MeanShapeCovariance");

            for (int i = 0; i < 128; i++) {
                List<double> sample = wald.GetRandomValues(rng, 16).ToList();
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
            Random rng = new Random(1);

            List<double> parameter = new List<double>();
            List<double> variance = new List<double>();

            for (int i = 0; i < 128; i++) {
                // We pick a quite-small sample, because we have a finite-n unbiased estimator.
                List<double> s = rayleigh.GetRandomValues(rng, 8).ToList();

                RayleighFitResult r = RayleighDistribution.FitToSample(s);
                parameter.Add(r.Scale.Value);
                variance.Add(r.Parameters.VarianceOf("Scale"));

                Assert.IsTrue(r.GoodnessOfFit.Probability > 0.01);
            }

            Assert.IsTrue(parameter.PopulationMean().ConfidenceInterval(0.99).ClosedContains(rayleigh.Scale));
            Assert.IsTrue(parameter.PopulationVariance().ConfidenceInterval(0.99).ClosedContains(variance.Median()));

        }

    }
}
