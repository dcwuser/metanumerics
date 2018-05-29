using System;
using System.Collections.Generic;
using System.Diagnostics;

using TestClassAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestClassAttribute;
using TestMethodAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestMethodAttribute;
using ExpectedExceptionAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.ExpectedExceptionAttribute;
using Assert = Microsoft.VisualStudio.TestTools.UnitTesting.Assert;

using Meta.Numerics;
using Meta.Numerics.Data;
using Meta.Numerics.Matrices;
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

            for (int i = 0; i < acv1.Length; i++) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(acv1[i], series1.Autocovariance(i)));
            }

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
            Sample tests = new Sample("p");
            for (int i = 0; i < 64; i++) {

                TimeSeries series = GenerateMA1TimeSeries(beta, mu, sigma, n, n * i + 314159);

                Debug.Assert(series.Count == n);

                MA1FitResult result = series.FitToMA1();

                //Assert.IsTrue(result.Dimension == 3);
                parameters.Add(result.Parameters.ValuesVector);
                covariances.Add(
                    result.Parameters.CovarianceMatrix[0, 0],
                    result.Parameters.CovarianceMatrix[1, 1],
                    result.Parameters.CovarianceMatrix[2, 2],
                    result.Parameters.CovarianceMatrix[0, 1],
                    result.Parameters.CovarianceMatrix[0, 2],
                    result.Parameters.CovarianceMatrix[1, 2]
                );
                tests.Add(result.GoodnessOfFit.Probability);

            }

            Assert.IsTrue(parameters.Column(0).PopulationMean.ConfidenceInterval(0.99).ClosedContains(mu)); ;
            Assert.IsTrue(parameters.Column(1).PopulationMean.ConfidenceInterval(0.99).ClosedContains(beta));
            Assert.IsTrue(parameters.Column(2).PopulationMean.ConfidenceInterval(0.99).ClosedContains(sigma));

            Assert.IsTrue(parameters.Column(0).PopulationVariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(0).Mean));
            Assert.IsTrue(parameters.Column(1).PopulationVariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(1).Mean));
            Assert.IsTrue(parameters.Column(2).PopulationVariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(2).Mean));
            Assert.IsTrue(parameters.TwoColumns(0, 1).PopulationCovariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(3).Mean));
            Assert.IsTrue(parameters.TwoColumns(0, 2).PopulationCovariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(4).Mean));
            Assert.IsTrue(parameters.TwoColumns(1, 2).PopulationCovariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(5).Mean));

            Assert.IsTrue(tests.KuiperTest(new UniformDistribution()).Probability > 0.01);

        }


        [TestMethod]
        public void TimeSeriesFitAR1() {

            double alpha = 0.3;
            double mu = 0.2;
            double sigma = 0.4;
            int n = 24;

            // For our fit to AR(1), we have incorporated bias correction (at least
            // for the most important parameter alpha), so we can do a small-n test.

            FrameTable data = new FrameTable();
            data.AddColumn<UncertainValue>("mu");
            data.AddColumn<UncertainValue>("alpha");
            data.AddColumn<UncertainValue>("sigma");
            data.AddColumn<SymmetricMatrix>("covariance");
            data.AddColumn<double>("p");

            for (int i = 0; i < 128; i++) {

                TimeSeries series = GenerateAR1TimeSeries(alpha, mu, sigma, n, n * i + 271828);

                AR1FitResult result = series.FitToAR1();

                data.AddRow(
                    result.Mu, result.Alpha, result.Sigma,
                    result.Parameters.CovarianceMatrix, result.GoodnessOfFit.Probability
                );

            }

            data.AddComputedColumn("alphaValue", r => ((UncertainValue) r["alpha"]).Value);
            data.AddComputedColumn("muValue", r => ((UncertainValue) r["mu"]).Value);

            // Check that fit parameters agree with inputs
            Assert.IsTrue(data["mu"].As((UncertainValue v) => v.Value).PopulationMean().ConfidenceInterval(0.99).ClosedContains(mu));
            Assert.IsTrue(data["alpha"].As((UncertainValue v) => v.Value).PopulationMean().ConfidenceInterval(0.99).ClosedContains(alpha));
            Assert.IsTrue(data["sigma"].As((UncertainValue v) => v.Value).PopulationMean().ConfidenceInterval(0.99).ClosedContains(sigma));

            // Check that reported variances agree with actual variances
            Assert.IsTrue(data["mu"].As((UncertainValue v) => v.Value).PopulationStandardDeviation().ConfidenceInterval(0.99).ClosedContains(data["mu"].As((UncertainValue v) => v.Uncertainty).Median()));
            Assert.IsTrue(data["alpha"].As((UncertainValue v) => v.Value).PopulationStandardDeviation().ConfidenceInterval(0.99).ClosedContains(data["alpha"].As((UncertainValue v) => v.Uncertainty).Median()));
            Assert.IsTrue(data["sigma"].As((UncertainValue v) => v.Value).PopulationStandardDeviation().ConfidenceInterval(0.99).ClosedContains(data["sigma"].As((UncertainValue v) => v.Uncertainty).Median()));

            // Check that reported co-variances agree with actual co-variances
            Assert.IsTrue(data["mu"].As((UncertainValue v) => v.Value).PopulationCovariance(data["alpha"].As((UncertainValue v) => v.Value)).ConfidenceInterval(0.99).ClosedContains(data["covariance"].As((SymmetricMatrix c) => c[0, 1]).Median()));

            // For small n, the fitted alpha can vary considerably, and the formula for var(m) varies
            // quite strongly with alpha, so the computed var(m) have a very long tail. This pushes the
            // mean computed var(m) quite a bit higher than a typical value, so we use medians instead
            // of means for our best guess for the predicted variance.

            TestResult ks = data["p"].As<double>().KolmogorovSmirnovTest(new UniformDistribution());
            Assert.IsTrue(ks.Probability > 0.05);

            // This is an onerous way to store values, but it does let us test how the data-frame machinery deals with
            // non-trivial storage types.

        }

        [TestMethod]
        public void TimeSeriesBadFit () {

            // Fit AR1 to MA1; the fit should be bad

            TimeSeries series = GenerateMA1TimeSeries(0.4, 0.3, 0.2, 1000);

            AR1FitResult result = series.FitToAR1();

            Assert.IsTrue(result.GoodnessOfFit.Probability < 0.01);

        }

        [TestMethod]
        public void LjungBoxNullDistribution () {

            Sample Qs = new Sample();
            ContinuousDistribution d = null;
            for (int i = 0; i < 100; i++) {
                TimeSeries series = GenerateMA1TimeSeries(0.0, 1.0, 2.0, 10, i);
                TestResult lbResult = series.LjungBoxTest(5);
                Qs.Add(lbResult.Statistic.Value);
                d = lbResult.Statistic.Distribution;
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
        public void DifferenceToStationarity () {

            int n = 100;
            // Generate white noise
            TimeSeries series = GenerateMA1TimeSeries(0.0, 0.1, 1.0, n);

            // Integrate it
            series.Integrate(0.0);
            //Assert.IsTrue(series.Count == n);
            // The result should be correlated
            TestResult result1 = series.LjungBoxTest();
            Assert.IsTrue(result1.Probability < 0.01);

            // Difference to return to original values
            series.Difference();
            //Assert.IsTrue(series.Count == (n - 1));
            // The result should be correlated
            TestResult result2 = series.LjungBoxTest();
            Assert.IsTrue(result2.Probability > 0.01);

        }

        [TestMethod]
        public void WhiteNoisePowerSpectrum () {

            TimeSeries series = GenerateMA1TimeSeries(0.0, 0.0, 2.0, 20);

            double[] spectrum = series.PowerSpectrum();

        }

    }
}
