using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.ObjectModel;
using System.Text;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Analysis;
using Meta.Numerics.Matrices;
using Meta.Numerics.SignalProcessing;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace Test {

    [TestClass]
    public class TimeSeriesTest {

        // For red noise, use AR(1) with mu = 0

        private TimeSeries GenerateAR1TimeSeries (double alpha, double mu, double sigma, int count, int seed = 1) {
            TimeSeries series = new TimeSeries();
            Random rng = new Random(seed);
            NormalDistribution d = new NormalDistribution(0.0, sigma);
            double previousDeviation = d.GetRandomValue(rng) / Math.Sqrt(1.0 - alpha * alpha);
            for (int i = 0; i < count; i++) {
                double currentDeviation = alpha * previousDeviation + d.GetRandomValue(rng);
                series.Add(mu + currentDeviation);
                previousDeviation = currentDeviation;
            }
            return (series);
        }

        // For white noise, use MA(1) with beta = mu = 0

        private TimeSeries GenerateMA1TimeSeries (double beta, double mu, double sigma, int count, int seed = 1) {
            TimeSeries series = new TimeSeries();
            Random rng = new Random(seed);
            NormalDistribution eDist = new NormalDistribution(0.0, sigma);
            double uPrevious = eDist.GetRandomValue(rng);
            for (int i = 0; i < count; i++) {
                double u = eDist.GetRandomValue(rng);
                series.Add(mu + u + beta * uPrevious);
                uPrevious = u;
            }
            return (series);
        }

        [TestMethod]
        public void TimeSeriesBasics () {

            TimeSeries series = new TimeSeries();

            Assert.IsTrue(series.Count == 0);

            series.Add(3.0, 2.0, 1.0);

            Assert.IsTrue(series[0] == 3.0);
            Assert.IsTrue(series[1] == 2.0);
            Assert.IsTrue(series[2] == 1.0);

            Assert.IsTrue(series.Count == 3);
            Assert.IsFalse(series.Contains(0.0));
            Assert.IsTrue(series.Contains(1.0));
            Assert.IsTrue(series.IndexOf(0.0) == -1);
            Assert.IsTrue(series.IndexOf(1.0) == 2);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(series.Mean, 2.0));

            series[2] = 0.0;
            Assert.IsTrue(series.Count == 3);
            Assert.IsTrue(series.Contains(0.0));
            Assert.IsFalse(series.Contains(1.0));
            Assert.IsTrue(series.IndexOf(0.0) == 2);
            Assert.IsTrue(series.IndexOf(1.0) == -1);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(series.Mean, 5.0 / 3.0));
            
            series.Clear();
            Assert.IsTrue(series.Count == 0);
            Assert.IsFalse(series.Contains(0.0));
            Assert.IsTrue(series.IndexOf(0.0) == -1);

        }

        [TestMethod]
        public void TimeSeriesSampleAgreement () {

            TimeSeries series = new TimeSeries(7.0, 5.0, 3.0, 2.0);

            Sample sample = series.AsSample();

            Assert.IsTrue(series.Count == sample.Count);
            Assert.IsTrue(series.Mean == sample.Mean);
            Assert.IsTrue(series.Autocovariance(0) == sample.Variance);

            Assert.IsTrue(sample.Contains(series[1]));

        }

        [TestMethod]
        public void TimeSeriesAutocovariance () {

            TimeSeries series1 = new TimeSeries();
            series1.Add(0.0, -1.0, 2.0, -3.0);

            double[] acv1 = series1.Autocovariance();
            TestAutocovariance(acv1);

        }

        private void TestAutocovariance (IList<double> acv) {

            // Verify that variance is positive
            Assert.IsTrue(acv[0] >= 0.0);

            // Verify that correlations -1 <= r <= 1
            for (int i = 0; i < acv.Count; i++) {
                Assert.IsTrue(Math.Abs(acv[i]) <= acv[0]);
            }

            // Verify that autocovariance is positive definite
            SymmetricMatrix C = new SymmetricMatrix(acv.Count);
            for (int r = 0; r < acv.Count; r++) {
                for (int c = 0; c <= r; c++) {
                    C[r, c] = acv[r - c];
                }
            }
            CholeskyDecomposition CD = C.CholeskyDecomposition();
            Assert.IsTrue(CD != null);
        }

        [TestMethod]
        public void TimeSeriesPopulationStatistics () {

            double alpha = -0.8;
            double beta = 0.6;
            double mu = -0.4;
            double sigma = 0.2;

            TimeSeries ar1 = GenerateAR1TimeSeries(alpha, mu, sigma, 100, 314159);
            TimeSeriesPopulationStatistics ar1ps = ar1.PopulationStatistics();
            Assert.IsTrue(ar1ps.Mean.ConfidenceInterval(0.99).ClosedContains(mu));

            // For AR1, autocovariance decreases by a factor \alpha for each additional lag step.
            double ar1v = sigma * sigma / (1.0 - alpha * alpha);
            Assert.IsTrue(ar1ps.Autocovariance(0).ConfidenceInterval(0.99).ClosedContains(ar1v));
            Assert.IsTrue(ar1ps.Autocovariance(1).ConfidenceInterval(0.99).ClosedContains(alpha * ar1v));
            Assert.IsTrue(ar1ps.Autocovariance(2).ConfidenceInterval(0.99).ClosedContains(alpha * alpha * ar1v));
            Assert.IsTrue(ar1ps.Autocovariance(3).ConfidenceInterval(0.99).ClosedContains(alpha * alpha * alpha * ar1v));

            TimeSeries ma1 = GenerateMA1TimeSeries(beta, mu, sigma, 100, 314159);
            TimeSeriesPopulationStatistics ma1ps = ma1.PopulationStatistics();
            Assert.IsTrue(ma1ps.Mean.ConfidenceInterval(0.99).ClosedContains(mu));

            // For MA(1), autocovariance vanishes for all lags greater than 1.
            Assert.IsTrue(ma1ps.Autocovariance(0).ConfidenceInterval(0.99).ClosedContains((1.0 + beta * beta) * sigma * sigma));
            Assert.IsTrue(ma1ps.Autocovariance(1).ConfidenceInterval(0.99).ClosedContains(beta * sigma * sigma));
            Assert.IsTrue(ma1ps.Autocovariance(2).ConfidenceInterval(0.99).ClosedContains(0.0));
            Assert.IsTrue(ma1ps.Autocovariance(3).ConfidenceInterval(0.99).ClosedContains(0.0));

            // Note we do not yet test that reported uncertainties accurately predict
            // variances. That is because our uncertainty estimates are not very rigorous.
            // When we make them more rigorous, we should test this.

        }

        [TestMethod]
        public void TimeSeriesFitToMA1 () {

            double beta = -0.2;
            double mu = 0.4;
            double sigma = 0.6;
            int n = 100;

            // If we are going to strictly test parameter values and variances,
            // we can't pick n too small, because the formulas we use are only
            // asymptotically unbiased.

            MultivariateSample parameters = new MultivariateSample(3);
            MultivariateSample covariances = new MultivariateSample(6);
            for (int i = 0; i < 100; i++) {

                TimeSeries series = GenerateMA1TimeSeries(beta, mu, sigma, n, i + 314159);

                Debug.Assert(series.Count == n);

                FitResult result = series.FitToMA1();

                Assert.IsTrue(result.Dimension == 3);

                parameters.Add(result.Parameters);
                covariances.Add(
                    result.CovarianceMatrix[0, 0],
                    result.CovarianceMatrix[1, 1],
                    result.CovarianceMatrix[2, 2],
                    result.CovarianceMatrix[0, 1],
                    result.CovarianceMatrix[0, 2],
                    result.CovarianceMatrix[1, 2]
                );

            }

            Assert.IsTrue(parameters.Column(0).PopulationMean.ConfidenceInterval(0.99).ClosedContains(beta)); ;
            Assert.IsTrue(parameters.Column(1).PopulationMean.ConfidenceInterval(0.99).ClosedContains(mu));
            Assert.IsTrue(parameters.Column(2).PopulationMean.ConfidenceInterval(0.99).ClosedContains(sigma));

            Assert.IsTrue(parameters.Column(0).PopulationVariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(0).Mean));
            Assert.IsTrue(parameters.Column(1).PopulationVariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(1).Mean));
            Assert.IsTrue(parameters.Column(2).PopulationVariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(2).Mean));
            Assert.IsTrue(parameters.TwoColumns(0, 1).PopulationCovariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(3).Mean));
            Assert.IsTrue(parameters.TwoColumns(0, 2).PopulationCovariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(4).Mean));
            Assert.IsTrue(parameters.TwoColumns(1, 2).PopulationCovariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(5).Mean));

        }


        [TestMethod]
        public void TimeSeriesFitAR1() {

            double alpha = 0.3;
            double mu = 0.2;
            double sigma = 0.4;
            int n = 20;

            // For our fit to AR(1), we have incorporated bias correction (at least
            // for the most important parameter alpha), so we can do a small-n test.

            MultivariateSample parameters = new MultivariateSample(3);
            MultivariateSample covariances = new MultivariateSample(6);
            for (int i = 0; i < 100; i++) {

                TimeSeries series = GenerateAR1TimeSeries(alpha, mu, sigma, n, i + 314159);

                FitResult result = series.FitToAR1();

                parameters.Add(result.Parameters);
                covariances.Add(
                    result.CovarianceMatrix[0, 0],
                    result.CovarianceMatrix[1, 1],
                    result.CovarianceMatrix[2, 2],
                    result.CovarianceMatrix[0, 1],
                    result.CovarianceMatrix[0, 2],
                    result.CovarianceMatrix[1, 2]
                );

            }

            // Check that fit parameters agree with inputs
            Assert.IsTrue(parameters.Column(0).PopulationMean.ConfidenceInterval(0.99).ClosedContains(alpha));
            Assert.IsTrue(parameters.Column(1).PopulationMean.ConfidenceInterval(0.99).ClosedContains(mu));
            Assert.IsTrue(parameters.Column(2).PopulationMean.ConfidenceInterval(0.99).ClosedContains(sigma));

            // Check that reported variances agree with actual variances
            Assert.IsTrue(parameters.Column(0).PopulationVariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(0).Median));
            Assert.IsTrue(parameters.Column(1).PopulationVariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(1).Median));
            Assert.IsTrue(parameters.Column(2).PopulationVariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(2).Median));
            Assert.IsTrue(parameters.TwoColumns(0, 1).PopulationCovariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(3).Mean));
            Assert.IsTrue(parameters.TwoColumns(0, 2).PopulationCovariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(4).Mean));
            Assert.IsTrue(parameters.TwoColumns(1, 2).PopulationCovariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(5).Mean));

            // For small n, the fitted alpha can vary considerably, and the formula for var(m) varies
            // quite strongly with alpha, so the computed var(m) have a very long tail. This pushes the
            // mean computed var(m) quite a bit higher than a typical value, so we use medians instead
            // of means for our best guess for the predicted variance.

        }

        [TestMethod]
        public void LjungBoxNullDistribution () {

            Sample Qs = new Sample();
            Distribution d = null;
            for (int i = 0; i < 100; i++) {
                TimeSeries series = GenerateMA1TimeSeries(0.0, 1.0, 2.0, 10, i);
                TestResult lbResult = series.LjungBoxTest(5);
                Qs.Add(lbResult.Statistic);
                d = lbResult.Distribution;
            }

            TestResult kResult = Qs.KuiperTest(d);
            Assert.IsTrue(kResult.Probability > 0.05);

        }


        [TestMethod]
        public void LjungBoxRejection () {

            TimeSeries series = GenerateAR1TimeSeries(0.3, 0.2, 0.1, 100);
            TestResult result = series.LjungBoxTest(10);
            Assert.IsTrue(result.Probability < 0.05);

        }

        [TestMethod]
        public void WhiteNoisePowerSpectrum () {

            TimeSeries series = GenerateMA1TimeSeries(0.0, 0.0, 2.0, 20);

            double[] spectrum = series.PowerSpectrum();

        }

    }
}
