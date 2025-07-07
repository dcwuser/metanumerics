using System;
using System.Collections.Generic;
using System.Linq;

using TestClassAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestClassAttribute;
using TestMethodAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestMethodAttribute;
using Assert = Microsoft.VisualStudio.TestTools.UnitTesting.Assert;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Data;
using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;
using System.Data.Common;
using System.Runtime.InteropServices;
using FluentAssertions;

namespace Test {

    [TestClass]
    public class BivariateSampleTest {

        [TestMethod]
        public void CovarianceCorrelationRelatioship () {

            // Make a random bivariate sample
            int n = 10;
            Random rng = new Random(1);
            ContinuousDistribution xDistribution = new NormalDistribution(1.0, 2.0);
            ContinuousDistribution yDistribution = new LognormalDistribution();
            List<double> x = xDistribution.GetRandomValues(rng, n).ToList();
            List<double> y = yDistribution.GetRandomValues(rng, n).ToList();

            // Compute correlation, covariance, and variances
            double xyCov = Bivariate.Covariance(x, y);
            double xVar = Univariate.Variance(x);
            double yVar = Univariate.Variance(y);
            double r = Bivariate.CorrelationCoefficient(x, y);

            r.Should().BeNearly(xyCov / Math.Sqrt(xVar * yVar));

        }

        [TestMethod]
        public void LinearLogisticRegressionSimple () {

            Polynomial m = Polynomial.FromCoefficients(-1.0, 2.0);

            FrameTable table = new FrameTable();
            table.AddColumn<double>("x");
            table.AddColumn<string>("z");

            Random rng = new Random(2);
            ContinuousDistribution xDistribution = new CauchyDistribution(4.0, 2.0);
            for (int i = 0; i < 24; i++) {
                double x = xDistribution.GetRandomValue(rng);
                double y = m.Evaluate(x);
                double p = 1.0 / (1.0 + Math.Exp(-y));
                bool z = (rng.NextDouble() < p);
                table.AddRow(x, z.ToString());
            }

            LinearLogisticRegressionResult fit = table["z"].As((string s) => Boolean.Parse(s)).LinearLogisticRegression(table["x"].As<double>());
            Assert.IsTrue(fit.Intercept.ConfidenceInterval(0.99).ClosedContains(m.Coefficient(0)));
            Assert.IsTrue(fit.Slope.ConfidenceInterval(0.99).ClosedContains(m.Coefficient(1)));
        }

        [TestMethod]
        public void LinearLogisticRegressionVariances () {

            // define model y = a + b0 * x0 + b1 * x1 + noise
            double a = -2.0;
            double b = 1.0;
            ContinuousDistribution xDistribution = new StudentDistribution(2.0);

            FrameTable data = new FrameTable();
            data.AddColumns<double>("a", "da", "b", "db", "abcov", "p", "dp");

            // draw a sample from the model
            Random rng = new Random(3);
            for (int j = 0; j < 32; j++) {
                List<double> xs = new List<double>();
                List<bool> ys = new List<bool>();

                for (int i = 0; i < 32; i++) {
                    double x = xDistribution.GetRandomValue(rng);
                    double t = a + b * x;
                    double p = 1.0 / (1.0 + Math.Exp(-t));
                    bool y = (rng.NextDouble() < p);
                    xs.Add(x);
                    ys.Add(y);
                }

                // do a linear regression fit on the model
                LinearLogisticRegressionResult result = ys.LinearLogisticRegression(xs);
                UncertainValue pp = result.Predict(1.0);

                data.AddRow(
                    result.Intercept.Value, result.Intercept.Uncertainty,
                    result.Slope.Value, result.Slope.Uncertainty,
                    result.Parameters.CovarianceMatrix[0, 1],
                    pp.Value, pp.Uncertainty
                );

            }

            // The estimated parameters should agree with the model that generated the data.

            // The variances of the estimates should agree with the claimed variances
            Assert.IsTrue(data["a"].As<double>().PopulationStandardDeviation().ConfidenceInterval(0.99).ClosedContains(data["da"].As<double>().Mean()));
            Assert.IsTrue(data["b"].As<double>().PopulationStandardDeviation().ConfidenceInterval(0.99).ClosedContains(data["db"].As<double>().Mean()));
            Assert.IsTrue(data["a"].As<double>().PopulationCovariance(data["b"].As<double>()).ConfidenceInterval(0.99).ClosedContains(data["abcov"].As<double>().Mean()));
            Assert.IsTrue(data["p"].As<double>().PopulationStandardDeviation().ConfidenceInterval(0.99).ClosedContains(data["dp"].As<double>().Mean()));

        }


        [TestMethod]
        public void LinearLogisticRegression () {

            // do a set of logistic regression fits
            // make sure not only that the fit parameters are what they should be, but that their variances/covariances are as returned

            Random rng = new Random(314159);

            // define logistic parameters
            double a0 = 1.0; double b0 = -0.5;
            //double a0 = -0.5; double b0 = 2.0;

            // keep track of sample of returned a and b fit parameters
            FrameTable t = new FrameTable();
            t.AddColumn<double>("a");
            t.AddColumn<double>("b");

            // also keep track of returned covariance estimates
            // since these vary slightly from fit to fit, we will average them
            double caa = 0.0;
            double cbb = 0.0;
            double cab = 0.0;

            // Do 50 fits with a sample size of 50 each.
            // We can re-use data storage array since they are done sequentially.
            double[] x = new double[50];
            bool[] y = new bool[50];
            for (int k = 0; k < 50; k++) {

                // Generate a synthetic data set
                for (int i = 0; i < x.Length; i++) {
                    x[i] = 2.0 * rng.NextDouble() - 1.0;
                    double ez = Math.Exp(a0 + b0 * x[i]);
                    double P = ez / (1.0 + ez);
                    y[i] = (rng.NextDouble() < P);
                }

                // Do the regression
                LinearLogisticRegressionResult r = y.LinearLogisticRegression(x);

                // record best fit parameters
                double a = r.Intercept.Value;
                double b = r.Slope.Value;
                t.AddRow(a, b);

                // Record estimated covariances
                caa += r.Parameters.CovarianceMatrix[0, 0];
                cbb += r.Parameters.CovarianceMatrix[1, 1];
                cab += r.Parameters.CovarianceMatrix[0, 1];

            }

            // The estimated covariances will vary somewhat, so average them.
            caa /= t.Rows.Count;
            cbb /= t.Rows.Count;
            cab /= t.Rows.Count;

            IReadOnlyList<double> aValues = t["a"].As<double>();
            IReadOnlyList<double> bValues = t["b"].As<double>();

            // check that mean parameter estimates are what they should be: the underlying population parameters
            Assert.IsTrue(aValues.PopulationMean().ConfidenceInterval(0.95).Contains(a0));
            Assert.IsTrue(bValues.PopulationMean().ConfidenceInterval(0.95).Contains(b0));

            // check that parameter covarainces are what they should be: the reported covariance estimates
            Assert.IsTrue(aValues.PopulationVariance().ConfidenceInterval(0.95).Contains(caa));
            Assert.IsTrue(bValues.PopulationVariance().ConfidenceInterval(0.95).Contains(cbb));
            Assert.IsTrue(Bivariate.PopulationCovariance(aValues, bValues).ConfidenceInterval(0.95).Contains(cab));

        }

        [TestMethod]
        public void LinearRegressionSimple () {

            double a = -1.0;
            double b = 2.0;

            ContinuousDistribution xDistribution = new CauchyDistribution();
            ContinuousDistribution eDistribution = new NormalDistribution();

            int n = 16;
            Random rng = new Random(1);
            double[] x = new double[n];
            double[] y = new double[n];
            for (int i = 0; i < 16; i++) {
                x[i] = xDistribution.GetRandomValue(rng);
                y[i] = a + b * x[i] + eDistribution.GetRandomValue(rng);
            }

            LinearRegressionResult result = y.LinearRegression(x);

            // Parameters should be right
            Assert.IsTrue(result.Intercept.ConfidenceInterval(0.95).ClosedContains(a));
            Assert.IsTrue(result.Slope.ConfidenceInterval(0.95).ClosedContains(b));

            // Reported values should be consistent
            Assert.IsTrue(result.Intercept == result.Parameters["Intercept"].Estimate);
            Assert.IsTrue(result.Intercept.Value == result.Parameters.ValuesVector[result.Parameters.IndexOf("Intercept")]);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Intercept.Uncertainty, Math.Sqrt(result.Parameters.VarianceOf("Intercept"))));
            Assert.IsTrue(result.Slope == result.Parameters["Slope"].Estimate);
            Assert.IsTrue(result.Slope.Value == result.Parameters.ValuesVector[result.Parameters.IndexOf("Slope")]);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Slope.Uncertainty, Math.Sqrt(result.Parameters.VarianceOf("Slope"))));

            // Residuals should agree with definition
            for (int i = 0; i < x.Length; i++) {
                double yp = result.Predict(x[i]).Value;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Residuals[i], y[i] - yp));
            }

            // R and R-squared agree
            Assert.IsTrue(TestUtilities.IsNearlyEqual(result.RSquared, MoreMath.Sqr(result.R.Statistic.Value)));

            // F-test and R-test agree
            Assert.IsTrue(TestUtilities.IsNearlyEqual(result.F.Probability, result.R.Probability));

            // ANOVA's sums of squares are correct
            double SST = y.Variance() * y.Length;
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SST, result.Anova.Total.SumOfSquares));
            double SSR = 0.0;
            foreach (double z in result.Residuals) SSR += z * z;
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SSR, result.Anova.Residual.SumOfSquares));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SSR, result.SumOfSquaredResiduals));

            // R is same as correlation coefficient
            Assert.IsTrue(TestUtilities.IsNearlyEqual(x.CorrelationCoefficient(y), result.R.Statistic.Value));

        }

        [TestMethod]
        public void LinearRegressionVariances () {

            // Do a set of logistic regression fits
            // Make sure not only that the fit parameters are what they should be, but that their variances/covariances are as returned
            // and that the variance of a predicted value is as returned.

            Random rng = new Random(314159);

            // define line parameters
            double a0 = 2.0; double b0 = -1.0;

            // do a lot of fits, recording results of each
            FrameTable data = new FrameTable();
            data.AddColumns<double>("a", "va", "b", "vb", "abCov", "p", "dp");

            for (int k = 0; k < 128; k++) {

                // we should be able to draw x's from any distribution; noise should be drawn from a normal distribution
                ContinuousDistribution xd = new LogisticDistribution();
                ContinuousDistribution nd = new NormalDistribution(0.0, 2.0);

                // generate a synthetic data set
                List<double> xs = new List<double>();
                List<double> ys = new List<double>();
                for (int i = 0; i < 12; i++) {
                    double x = xd.GetRandomValue(rng);
                    xs.Add(x);
                    double y = a0 + b0 * x + nd.GetRandomValue(rng);
                    ys.Add(y);
                }

                // do the regression
                LinearRegressionResult result = ys.LinearRegression(xs);

                // record result
                UncertainValue p = result.Predict(12.0);
                data.AddRow(new Dictionary<string, object>() {
                    {"a", result.Intercept.Value },
                    {"va", result.Parameters.VarianceOf("Intercept") },
                    {"b", result.Slope.Value},
                    {"vb", result.Parameters.VarianceOf("Slope") },
                    {"abCov", result.Parameters.CovarianceOf("Slope", "Intercept") },
                    {"p", p.Value},
                    {"dp", p.Uncertainty }
                });
            }

            // variances of parameters should agree with predictions
            Assert.IsTrue(data["a"].As<double>().PopulationVariance().ConfidenceInterval(0.99).Contains(data["va"].As<double>().Median()));
            Assert.IsTrue(data["b"].As<double>().PopulationVariance().ConfidenceInterval(0.99).Contains(data["vb"].As<double>().Median()));
            Assert.IsTrue(data["a"].As<double>().PopulationCovariance(data["b"].As<double>()).ConfidenceInterval(0.99).Contains(data["abCov"].As<double>().Median()));

            // variance of prediction should agree with claim
            Assert.IsTrue(data["p"].As<double>().PopulationStandardDeviation().ConfidenceInterval(0.99).Contains(data["dp"].As<double>().Median()));

        }

        [TestMethod]
        public void BivariateLinearRegressionNullDistribution () {

            // Create uncorrelated x and y values and do a linear fit.
            // The r-tests and F-test statistics returned by the linear fits
            // should agree and both test statistics should follow their claimed
            // distributions.

            Random rng = new Random(987654321);
            NormalDistribution xd = new NormalDistribution(1.0, 2.0);
            NormalDistribution yd = new NormalDistribution(-3.0, 4.0);

            List<double> rSample = new List<double>();
            ContinuousDistribution rDistribution = null;

            List<double> fSample = new List<double>();
            ContinuousDistribution fDistribution = null;

            for (int i = 0; i < 127; i++) {

                List<double> x = new List<double>();
                List<double> y = new List<double>();
                for (int j = 0; j < 7; j++) {
                    x.Add(xd.GetRandomValue(rng));
                    y.Add(yd.GetRandomValue(rng));
                }
                LinearRegressionResult result = y.LinearRegression(x);

                rSample.Add(result.R.Statistic.Value);
                rDistribution = result.R.Statistic.Distribution;

                fSample.Add(result.F.Statistic.Value);
                fDistribution = result.F.Statistic.Distribution;

                Assert.IsTrue(result.F.Statistic.Value == result.Anova.Result.Statistic.Value);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    result.R.Probability, result.F.Probability,
                    new EvaluationSettings() { RelativePrecision = 1.0E-13, AbsolutePrecision = 1.0E-16 }
                ));

            }

            Assert.IsTrue(rSample.KuiperTest(rDistribution).Probability > 0.05);
            Assert.IsTrue(fSample.KuiperTest(fDistribution).Probability > 0.05);

        }

        [TestMethod]
        public void BivariateRegressionRSquaredDistribution() {

            // Do a bunch of linear regressions. r^2 should be distributed as expected.

            double a0 = 1.0;
            double b0 = 0.0;

            Random rng = new Random(1001110000);
            ContinuousDistribution xDistribution = new UniformDistribution(Interval.FromEndpoints(-2.0, 4.0));
            ContinuousDistribution eDistribution = new NormalDistribution();

            List<double> rSquared = new List<double>();
            for (int i = 0; i < 500; i++) {
                double[] x = new double[10];
                double[] y = new double[10];
                for (int k = 0; k < x.Length; k++) {
                    x[k] = xDistribution.GetRandomValue(rng);
                    y[k] = a0 + b0 * x[k] + eDistribution.GetRandomValue(rng);
                }
                LinearRegressionResult fit = y.LinearRegression(x);
                rSquared.Add(fit.RSquared);
            }

            ContinuousDistribution rSquaredDistribution = new BetaDistribution((2 - 1) / 2.0, (10 - 2) / 2.0);
            TestResult ks = rSquared.KolmogorovSmirnovTest(rSquaredDistribution);
            Assert.IsTrue(ks.Probability > 0.05);
        }

        [TestMethod]
        public void PolynomialRegressionSimple () {

            // Pick a simple polynomial
            Polynomial p = Polynomial.FromCoefficients(3.0, -2.0, 1.0);

            // Use it to generate a data set
            Random rng = new Random(1);
            ContinuousDistribution xDistribution = new CauchyDistribution(1.0, 2.0);
            ContinuousDistribution errorDistribution = new NormalDistribution(0.0, 3.0);
            List<double> xs = new List<double>(TestUtilities.CreateDataSample(rng, xDistribution, 10));
            List<double> ys = new List<double>(xs.Select(x => p.Evaluate(x) + errorDistribution.GetRandomValue(rng)));

            PolynomialRegressionResult fit = Bivariate.PolynomialRegression(ys, xs, p.Degree);

            // Parameters should agree
            Assert.IsTrue(fit.Parameters.Count == p.Degree + 1);
            for (int k = 0; k <= p.Degree; k++) {
                Assert.IsTrue(fit.Coefficient(k).ConfidenceInterval(0.99).ClosedContains(p.Coefficient(k)));
            }

            // Residuals should agree
            Assert.IsTrue(fit.Residuals.Count == xs.Count);
            for (int i = 0; i < xs.Count; i++) {
                double z = ys[i] - fit.Predict(xs[i]).Value;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(z, fit.Residuals[i]));
            }

            // Intercept is same as coefficient of x^0
            Assert.IsTrue(fit.Intercept == fit.Coefficient(0));

        }


        [TestMethod]
        public void PolynomialRegressionCovariance () {

            // do a set of polynomial regression fits
            // make sure not only that the fit parameters are what they should be, but that their variances/covariances are as claimed

            Random rng = new Random(271828);

            // define logistic parameters
            double[] a = new double[] { 0.0, -1.0, 2.0, -3.0 };

            // keep track of sample of returned a and b fit parameters
            List<ColumnVector> ps = new List<ColumnVector>();

            // also keep track of returned covariance estimates
            // since these vary slightly from fit to fit, we will average them
            SymmetricMatrix C = new SymmetricMatrix(a.Length);

            // do 100 fits
            for (int k = 0; k < 100; k++) {

                // we should be able to draw x's from any distribution; noise should be drawn from a normal distribution
                ContinuousDistribution xd = new CauchyDistribution();
                ContinuousDistribution nd = new NormalDistribution(0.0, 4.0);

                // generate a synthetic data set

                List<double> xs = new List<double>();
                List<double> ys = new List<double>();
                for (int j = 0; j < 20; j++) {
                    double x = xd.GetRandomValue(rng);
                    double y = nd.GetRandomValue(rng);
                    for (int i = 0; i < a.Length; i++) {
                        y += a[i] * MoreMath.Pow(x, i);
                    }
                    xs.Add(x);
                    ys.Add(y);
                }

                // do the regression
                PolynomialRegressionResult r = ys.PolynomialRegression(xs, a.Length - 1);

                // record the fit parameters
                ps.Add(r.Parameters.ValuesVector);

                // record estimated covariances
                C += r.Parameters.CovarianceMatrix;

            }

            C /= ps.Count;

            // check that mean parameter estimates are what they should be: the underlying population parameters
            for (int i = 0; i < a.Length; i++) {
                Assert.IsTrue(ps.Select(t => t[i]).ToList().PopulationMean().ConfidenceInterval(0.95).ClosedContains(a[i]));
            }

            // check that parameter covarainces are what they should be: the reported covariance estimates
            for (int i = 0; i < a.Length; i++) {
                for (int j = i; j < a.Length; j++) {
                    Assert.IsTrue(Bivariate.PopulationCovariance(ps.Select(t => t[i]).ToList(), ps.Select(t => t[j]).ToList()).ConfidenceInterval(0.95).ClosedContains(C[i, j]));
                }
            }

        }

        public void PolynomialRegressionFNullDistribution () {

            // Null distribution assumes all coefficients zero;

            double[] x = TestUtilities.GenerateUniformRealValues(0.1, 10.0, 10);
            double[] y = new double[x.Length];

            Random rng = new Random(1);
            NormalDistribution err = new NormalDistribution(0.0, 3.0);
            List<double> F = new List<double>();
            ContinuousDistribution FDistribution = null;
            for (int i = 0; i < 100; i ++) {

                for (int j = 0; j < y.Length; j++) {
                    y[j] = -2.0 + err.GetRandomValue(rng);
                }

                PolynomialRegressionResult r = y.PolynomialRegression(x, 4);
                F.Add(r.F.Statistic.Value);
                FDistribution = r.F.Statistic.Distribution;

            }

            Assert.IsTrue(F.KolmogorovSmirnovTest(FDistribution).Probability > 0.05);

        }

        [TestMethod]
        public void PolynomialRegressionLinearRegressionAgreement () {

            // A degree-1 polynomial fit should give the same answer as a linear fit
            double[] x = new double[] { 0.0, 3.0, 1.0, 4.0, 2.0 };
            double[] y = new double[] { 5.0, 6.0, 7.0, 8.0, 9.0 };
            GeneralLinearRegressionResult PR = y.PolynomialRegression(x, 1);
            GeneralLinearRegressionResult LR = y.LinearRegression(x);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(PR.Parameters.ValuesVector, LR.Parameters.ValuesVector));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(PR.Parameters.CovarianceMatrix, LR.Parameters.CovarianceMatrix));

        }

        [TestMethod]
        public void BivariateNonlinearFitSimple () {

            double t0 = 3.0;
            double s0 = 1.0;

            ContinuousDistribution xDistribution = new CauchyDistribution(0.0, 2.0);
            ContinuousDistribution eDistribution = new NormalDistribution(0.0, 0.5);

            Random rng = new Random(5);
            List<double> x = TestUtilities.CreateDataSample(rng, xDistribution, 48).ToList();
            List<double> y = x.Select(z => Math.Sin(2.0 * Math.PI * z / t0 + s0) + eDistribution.GetRandomValue(rng)).ToList();

            Func<IReadOnlyDictionary<string, double>, double, double> fitFunction = (d, z) => {
                double t = d["Period"];
                double s = d["Phase"];
                return (Math.Sin(2.0 * Math.PI * z / t + s));
            };

            Dictionary<string, double> start = new Dictionary<string, double>() {
                {"Period", 2.5 }, {"Phase", 1.5 }
            };

            NonlinearRegressionResult result = y.NonlinearRegression(x, fitFunction, start);

            Assert.IsTrue(result.Parameters["Period"].Estimate.ConfidenceInterval(0.99).ClosedContains(t0));
            Assert.IsTrue(result.Parameters["Phase"].Estimate.ConfidenceInterval(0.99).ClosedContains(s0));

            for (int i = 0; i < x.Count; i++) {
                double yp = result.Predict(x[i]);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Residuals[i], y[i] - yp));
            }

        }

        [TestMethod]
        public void BivariateNonlinearFitVariances () {

            // Verify that we can fit a non-linear function,
            // that the estimated parameters do cluster around the true values,
            // and that the estimated parameter covariances do reflect the actually observed covariances

            double a = 2.7;
            double b = 3.1;

            ContinuousDistribution xDistribution = new ExponentialDistribution(2.0);
            ContinuousDistribution eDistribution = new NormalDistribution(0.0, 4.0);

            FrameTable parameters = new FrameTable();
            parameters.AddColumns<double>("a", "b");

            // We want to compare (co)-variances of parameters to claimed covariance matrix, but the covariance
            // matrix changes from fit to fit. Average it.
            SymmetricMatrix C = new SymmetricMatrix(2);

            for (int i = 0; i < 64; i++) {

                List<double> xs = new List<double>();
                List<double> ys = new List<double>();
                Random rng = new Random(i);
                for (int j = 0; j < 8; j++) {
                    double x = xDistribution.GetRandomValue(rng);
                    xs.Add(x);
                    double y = a * Math.Pow(x, b) + eDistribution.GetRandomValue(rng);
                    ys.Add(y);
                }

                NonlinearRegressionResult fit = ys.NonlinearRegression(xs,
                    (IReadOnlyList<double> p, double x) => p[0] * Math.Pow(x, p[1]),
                    new double[] { 1.0, 1.0 }
                );

                parameters.AddRow(fit.Parameters.ValuesVector);
                C += fit.Parameters.CovarianceMatrix;

            }

            C /= parameters.Rows.Count;

            Assert.IsTrue(parameters["a"].As<double>().PopulationMean().ConfidenceInterval(0.99).Contains(a));
            Assert.IsTrue(parameters["b"].As<double>().PopulationMean().ConfidenceInterval(0.99).Contains(b));

            Assert.IsTrue(parameters["a"].As<double>().PopulationVariance().ConfidenceInterval(0.99).Contains(C[0, 0]));
            Assert.IsTrue(parameters["b"].As<double>().PopulationVariance().ConfidenceInterval(0.99).Contains(C[1, 1]));
            Assert.IsTrue(parameters["a"].As<double>().PopulationCovariance(parameters["b"].As<double>()).ConfidenceInterval(0.99).Contains(C[0, 1]));
            Assert.IsTrue(Bivariate.PopulationCovariance(parameters["a"].As<double>(), parameters["b"].As<double>()).ConfidenceInterval(0.99).Contains(C[0, 1]));
        }

        [TestMethod]
        public void BivariateAssociationDiscreteNullDistribution () {
            Random rng = new Random(1);

            // Pick very non-normal distributions for our non-parameteric tests
            ContinuousDistribution xd = new FrechetDistribution(1.0);
            ContinuousDistribution yd = new CauchyDistribution();

            // Pick small sample sizes to get exact distributions
            foreach (int n in TestUtilities.GenerateIntegerValues(4, 24, 4)) {

                // Do a bunch of test runs, recording reported statistic for each.
                List<int> spearmanStatistics = new List<int>();
                List<int> kendallStatistics = new List<int>();
                DiscreteDistribution spearmanDistribution = null;
                DiscreteDistribution kendallDistribution = null;

                for (int i = 0; i < 512; i++) {
                    List<double> x = new List<double>();
                    List<double> y = new List<double>();
                    for (int j = 0; j < n; j++) {
                        x.Add(xd.GetRandomValue(rng));
                        y.Add(yd.GetRandomValue(rng));
                    }

                    DiscreteTestStatistic spearman = Bivariate.SpearmanRhoTest(x, y).UnderlyingStatistic;
                    if (spearman != null) {
                        spearmanStatistics.Add(spearman.Value);
                        spearmanDistribution = spearman.Distribution;
                    }
                    DiscreteTestStatistic kendall = Bivariate.KendallTauTest(x, y).UnderlyingStatistic;
                    if (kendall != null) {
                        kendallStatistics.Add(kendall.Value);
                        kendallDistribution = kendall.Distribution;
                    }
                }

                // Test whether statistics are actually distributed as claimed
                if (spearmanDistribution != null) {
                    TestResult spearmanChiSquared = spearmanStatistics.ChiSquaredTest(spearmanDistribution);
                    Assert.IsTrue(spearmanChiSquared.Probability > 0.01);
                }
                if (kendallDistribution != null) {
                    TestResult kendallChiSquared = kendallStatistics.ChiSquaredTest(kendallDistribution);
                    Assert.IsTrue(kendallChiSquared.Probability > 0.01);
                }

            }
        }

        [TestMethod]
        public void BivariateDimensionMismatch () {

            double[] x = new double[3];
            double[] y = new double[4];

            Assert.ThrowsException<DimensionMismatchException>(() => x.CorrelationCoefficient(y));
            Assert.ThrowsException<DimensionMismatchException>(() => x.Covariance(y));
            Assert.ThrowsException<DimensionMismatchException>(() => y.LinearRegression(x));

        }

    }
}
