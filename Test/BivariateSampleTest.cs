using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Data;
using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace Test {

    [TestClass]
    public class BivariateSampleTest {

        [TestMethod]
        public void BivariateSampleManipulations () {

            BivariateSample s = new BivariateSample();
            s.Add(1.0, 9.0);
            s.Add(new XY(2.0, 8.0));
            s.Add(new double[] { 3.0, 4.0 }, new double[] { 7.0, 6.0 });
            s.Add(new XY[] { new XY(5.0, 5.0), new XY(6.0, 4.0) });
            Assert.IsTrue(s.Count == 6);

            Assert.IsTrue(!s.X.Contains(9.0));
            s.TransposeXY();
            Assert.IsTrue(s.X.Contains(9.0));
            s.TransposeXY();

            Assert.IsTrue(s.Remove(2.0, 8.0));
            Assert.IsTrue(s.Count == 5);
            Assert.IsFalse(s.Remove(2.0, 8.0));
            Assert.IsTrue(s.Count == 5);
            Assert.IsTrue(s.Remove(new XY(6.0, 4.0)));
            Assert.IsTrue(s.Count == 4);

            Assert.IsTrue(s.Contains(1.0, 9.0));
            Assert.IsFalse(s.Contains(9.0, 1.0));
            Assert.IsTrue(s.Contains(new XY(4.0, 6.0)));

            s.Clear();
            Assert.IsTrue(s.Count == 0);

        }

        [TestMethod]
        public void BivariateSampleEnumerations () {

            List<XY> points = new List<XY>(new XY[] { new XY(1.0, 2.0), new XY(2.0, 3.0), new XY(3.0, 4.0) });

            BivariateSample sample = new BivariateSample();
            sample.Add(points);

            Assert.IsTrue(sample.Count == points.Count);

            foreach (XY point in sample) {
                Assert.IsTrue(points.Remove(point));
            }

            Assert.IsTrue(points.Count == 0);
        }

        [TestMethod]
        public void BivariateSampleCopy () {

            // test independency of copy

            BivariateSample sample1 = new BivariateSample();
            sample1.Add(1.0, 2.0);

            BivariateSample sample2 = sample1.Copy();
            sample2.Add(3.0, 4.0);

            Assert.IsTrue(sample1.Count == 1);
            Assert.IsTrue(sample2.Count == 2);

        }

        [TestMethod]
        public void LinearLogisticRegressionSimple () {

            Polynomial m = Polynomial.FromCoefficients(-1.0, 2.0);

            FrameTable table = new FrameTable(
                new ColumnDefinition<double>("x"),
                new ColumnDefinition<string>("z")
            );
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

            FrameTable data = new FrameTable(
                new ColumnDefinition<double>("a"),
                new ColumnDefinition<double>("da"),
                new ColumnDefinition<double>("b"),
                new ColumnDefinition<double>("db"),
                new ColumnDefinition<double>("abcov"),
                new ColumnDefinition<double>("p"),
                new ColumnDefinition<double>("dp")
            );

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
                    result.Parameters.Covariance[0, 1],
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
            double a0 = 1.0; double b0 = -1.0 / 2.0;
            //double a0 = -0.5; double b0 = 2.0;

            // keep track of sample of returned a and b fit parameters
            BivariateSample ps = new BivariateSample();

            // also keep track of returned covariance estimates
            // since these vary slightly from fit to fit, we will average them
            double caa = 0.0;
            double cbb = 0.0;
            double cab = 0.0;

            // do 50 fits
            for (int k = 0; k < 50; k++) {

                Console.WriteLine("k={0}", k);

                // generate a synthetic data set
                BivariateSample s = new BivariateSample();
                for (int i = 0; i < 50; i++) {
                    double x = 2.0 * rng.NextDouble() - 1.0;
                    double ez = Math.Exp(a0 + b0 * x);
                    double P = ez / (1.0 + ez);
                    if (rng.NextDouble() < P) {
                        s.Add(x, 1.0);
                    } else {
                        s.Add(x, 0.0);
                    }
                }

                // do the regression
                LinearLogisticRegressionResult r = s.LinearLogisticRegression();

                // record best fit parameters
                double a = r.Intercept.Value;
                double b = r.Slope.Value;
                ps.Add(a, b);

                Console.WriteLine("{0}, {1}", a, b);

                // record estimated covariances
                caa += r.Parameters.Covariance[0, 0];
                cbb += r.Parameters.Covariance[1, 1];
                cab += r.Parameters.Covariance[0, 1]; 

            }

            caa /= ps.Count;
            cbb /= ps.Count;
            cab /= ps.Count;

            // check that mean parameter estimates are what they should be: the underlying population parameters
            Assert.IsTrue(ps.X.PopulationMean.ConfidenceInterval(0.95).ClosedContains(a0));
            Assert.IsTrue(ps.Y.PopulationMean.ConfidenceInterval(0.95).ClosedContains(b0));

            // check that parameter covarainces are what they should be: the reported covariance estimates
            Assert.IsTrue(ps.X.PopulationVariance.ConfidenceInterval(0.95).ClosedContains(caa));
            Assert.IsTrue(ps.Y.PopulationVariance.ConfidenceInterval(0.95).ClosedContains(cbb));
            Assert.IsTrue(ps.PopulationCovariance.ConfidenceInterval(0.95).ClosedContains(cab));

        }

        [TestMethod]
        public void BivariateLinearRegression () {

            // do a set of logistic regression fits
            // make sure not only that the fit parameters are what they should be, but that their variances/covariances are as returned

            Random rng = new Random(314159);

            // define line parameters
            double a0 = 2.0; double b0 = -1.0;

            // keep track of sample of returned a and b fit parameters
            BivariateSample pSample = new BivariateSample();

            // also keep track of returned covariance estimates
            // since these vary slightly from fit to fit, we will average them
            double caa = 0.0;
            double cbb = 0.0;
            double cab = 0.0;

            // Record predictions for a new point
            double x0 = 12.0;
            Sample ySample = new Sample();
            double ySigma = 0.0;

            // do 100 fits
            for (int k = 0; k < 128; k++) {

                // we should be able to draw x's from any distribution; noise should be drawn from a normal distribution
                ContinuousDistribution xd = new LogisticDistribution();
                ContinuousDistribution nd = new NormalDistribution(0.0, 2.0);

                // generate a synthetic data set
                BivariateSample sample = new BivariateSample();
                for (int i = 0; i < 16; i++) {
                    double x = xd.GetRandomValue(rng);
                    double y = a0 + b0 * x + nd.GetRandomValue(rng);
                    sample.Add(x, y);
                }

                // do the regression
                LinearRegressionResult result = sample.LinearRegression();

                // test consistancy
                Assert.IsTrue(result.Intercept == result.Parameters["Intercept"].Estimate);
                Assert.IsTrue(result.Intercept.Value == result.Parameters.Best[result.Parameters.IndexOf("Intercept")]);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Intercept.Uncertainty, Math.Sqrt(result.Parameters.CovarianceOf("Intercept", "Intercept"))));
                Assert.IsTrue(result.Slope == result.Parameters["Slope"].Estimate);
                Assert.IsTrue(result.Slope.Value == result.Parameters.Best[result.Parameters.IndexOf("Slope")]);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Slope.Uncertainty, Math.Sqrt(result.Parameters.CovarianceOf("Slope","Slope"))));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(result.R.Statistic, sample.CorrelationCoefficient));

                // Named parameters
                Assert.IsTrue(result.Parameters[nameof(result.Intercept)].Estimate == result.Intercept);
                Assert.IsTrue(result.Parameters[nameof(result.Slope)].Estimate == result.Slope);
                Assert.IsTrue(result.Parameters.CovarianceOf(nameof(result.Intercept), nameof(result.Slope)) == result.Parameters.Covariance[0, 1]);

                // record best fit parameters
                double a = result.Parameters.Best[0];
                double b = result.Parameters.Best[1];
                pSample.Add(a, b);

                // record estimated covariances
                caa += result.Parameters.Covariance[0, 0];
                cbb += result.Parameters.Covariance[1, 1];
                cab += result.Parameters.Covariance[0, 1];

                UncertainValue yPredict = result.Predict(x0);
                ySample.Add(yPredict.Value);
                ySigma += yPredict.Uncertainty;

                double SST = 0.0;
                foreach (double y in sample.Y) SST += MoreMath.Sqr(y - sample.Y.Mean);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(SST, result.Anova.Total.SumOfSquares));

                double SSR = 0.0;
                foreach (double z in result.Residuals) SSR += z * z;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(SSR, result.Anova.Residual.SumOfSquares));

            }

            caa /= pSample.Count;
            cbb /= pSample.Count;
            cab /= pSample.Count;
            ySigma /= pSample.Count;

            // check that mean parameter estimates are what they should be: the underlying population parameters
            //Assert.IsTrue(pSample.X.PopulationMean.ConfidenceInterval(0.95).ClosedContains(a0));
            //Assert.IsTrue(pSample.Y.PopulationMean.ConfidenceInterval(0.95).ClosedContains(b0));

            //Console.WriteLine("{0} {1}", caa, pSample.X.PopulationVariance);
            //Console.WriteLine("{0} {1}", cbb, pSample.Y.PopulationVariance);

            // check that parameter covarainces are what they should be: the reported covariance estimates
            Assert.IsTrue(pSample.X.PopulationVariance.ConfidenceInterval(0.95).ClosedContains(caa));
            Assert.IsTrue(pSample.Y.PopulationVariance.ConfidenceInterval(0.95).ClosedContains(cbb));
            Assert.IsTrue(pSample.PopulationCovariance.ConfidenceInterval(0.95).ClosedContains(cab));

            // Check that the predicted ys conform to the model and the asserted uncertainty.
            Assert.IsTrue(ySample.PopulationMean.ConfidenceInterval(0.95).ClosedContains(a0 + x0 * b0));
            //Assert.IsTrue(ySample.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(ySigma));

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

            Sample rSample = new Sample();
            ContinuousDistribution rDistribution = null;

            Sample fSample = new Sample();
            ContinuousDistribution fDistribution = null;

            for (int i = 0; i < 127; i++) {

                BivariateSample sample = new BivariateSample();
                for (int j = 0; j < 7; j++) {
                    sample.Add(xd.GetRandomValue(rng), yd.GetRandomValue(rng));
                }
                LinearRegressionResult result = sample.LinearRegression();

                rSample.Add(result.R.Statistic);
                rDistribution = result.R.Distribution;

                fSample.Add(result.F.Statistic);
                fDistribution = result.F.Distribution;

                Assert.IsTrue(result.F.Statistic == result.Anova.Result.Statistic);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    result.R.Probability, result.F.Probability,
                    new EvaluationSettings() { RelativePrecision = 1.0E-14, AbsolutePrecision = 1.0E-16 }
                ));

            }

            Assert.IsTrue(rSample.KuiperTest(rDistribution).Probability > 0.05);
            Assert.IsTrue(fSample.KuiperTest(fDistribution).Probability > 0.05);

        }

        [TestMethod]
        public void BivariatePolynomialRegressionFit () {

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

        }


        [TestMethod]
        public void BivariatePolynomialRegressionCovariance () {

            // do a set of polynomial regression fits
            // make sure not only that the fit parameters are what they should be, but that their variances/covariances are as claimed

            Random rng = new Random(271828);

            // define logistic parameters
            double[] a = new double[] { 0.0, -1.0, 2.0, -3.0 };

            // keep track of sample of returned a and b fit parameters
            MultivariateSample A = new MultivariateSample(a.Length);

            // also keep track of returned covariance estimates
            // since these vary slightly from fit to fit, we will average them
            SymmetricMatrix C = new SymmetricMatrix(a.Length);

            // also keep track of test statistics
            Sample F = new Sample();

            // do 100 fits
            for (int k = 0; k < 100; k++) {

                // we should be able to draw x's from any distribution; noise should be drawn from a normal distribution
                ContinuousDistribution xd = new CauchyDistribution();
                ContinuousDistribution nd = new NormalDistribution(0.0, 4.0);

                // generate a synthetic data set
                BivariateSample s = new BivariateSample();
                for (int j = 0; j < 20; j++) {
                    double x = xd.GetRandomValue(rng);
                    double y = nd.GetRandomValue(rng);
                    for (int i = 0; i < a.Length; i++) {
                        y += a[i] * MoreMath.Pow(x, i);
                    }
                    s.Add(x, y);
                }

                // do the regression
                PolynomialRegressionResult r = s.PolynomialRegression(a.Length - 1);

                ColumnVector ps = r.Parameters.Best;
                //Console.WriteLine("{0} {1} {2}", ps[0], ps[1], ps[2]);

                // record best fit parameters
                A.Add(ps);

                // record estimated covariances
                C += r.Parameters.Covariance;

                // record the fit statistic
                F.Add(r.F.Statistic);
                //Console.WriteLine("F={0}", r.GoodnessOfFit.Statistic);

            }

            C = (1.0 / A.Count) * C; // allow matrix division by real numbers

            // check that mean parameter estimates are what they should be: the underlying population parameters
            for (int i = 0; i < A.Dimension; i++) {
                Console.WriteLine("{0} {1}", A.Column(i).PopulationMean, a[i]);
                Assert.IsTrue(A.Column(i).PopulationMean.ConfidenceInterval(0.95).ClosedContains(a[i]));
            }

            // check that parameter covarainces are what they should be: the reported covariance estimates
            for (int i = 0; i < A.Dimension; i++) {
                for (int j = i; j < A.Dimension; j++) {
                    Console.WriteLine("{0} {1} {2} {3}", i, j, C[i, j], A.TwoColumns(i, j).PopulationCovariance);
                    Assert.IsTrue(A.TwoColumns(i, j).PopulationCovariance.ConfidenceInterval(0.95).ClosedContains(C[i, j]));
                }
            }

            // check that F is distributed as it should be
            //Console.WriteLine(fs.KolmogorovSmirnovTest(new FisherDistribution(2, 48)).LeftProbability);

        }

        [TestMethod]
        public void BivariateLinearPolynomialRegressionAgreement () {

            // A degree-1 polynomial fit should give the same answer as a linear fit

            BivariateSample B = new BivariateSample();
            B.Add(0.0, 5.0);
            B.Add(3.0, 6.0);
            B.Add(1.0, 7.0);
            B.Add(4.0, 8.0);
            B.Add(2.0, 9.0);
            GeneralLinearRegressionResult PR = B.PolynomialRegression(1);
            GeneralLinearRegressionResult LR = B.LinearRegression();
            Assert.IsTrue(TestUtilities.IsNearlyEqual(PR.Parameters.Best, LR.Parameters.Best));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(PR.Parameters.Covariance, LR.Parameters.Covariance));

        }

        [TestMethod]
        public void BivariateNonlinearFit () {

            // Verify that we can fit a non-linear function,
            // that the estimated parameters do cluster around the true values,
            // and that the estimated parameter covariances do reflect the actually observed covariances

            double a = 2.7;
            double b = 3.1;

            ContinuousDistribution xDistribution = new ExponentialDistribution(2.0);
            ContinuousDistribution eDistribution = new NormalDistribution(0.0, 4.0);

            FrameTable parameters = new FrameTable(new ColumnDefinition<double>("a"), new ColumnDefinition<double>("b"));
            //MultivariateSample parameters = new MultivariateSample("a", "b");
            MultivariateSample covariances = new MultivariateSample(3);

            for (int i = 0; i < 64; i++) {

                BivariateSample sample = new BivariateSample();
                Random rng = new Random(i);
                for (int j = 0; j < 8; j++) {
                    double x = xDistribution.GetRandomValue(rng);
                    double y = a * Math.Pow(x, b) + eDistribution.GetRandomValue(rng);
                    sample.Add(x, y);
                }

                NonlinearRegressionResult fit = sample.NonlinearRegression(
                    (IReadOnlyList<double> p, double x) => p[0] * Math.Pow(x, p[1]),
                    new double[] { 1.0, 1.0 }
                );

                parameters.AddRow(fit.Parameters.Best);
                //parameters.Add(fit.Parameters.Best);
                covariances.Add(fit.Parameters.Covariance[0, 0], fit.Parameters.Covariance[1, 1], fit.Parameters.Covariance[0, 1]);

            }

            Assert.IsTrue(parameters["a"].As<double>().PopulationMean().ConfidenceInterval(0.99).ClosedContains(a));
            Assert.IsTrue(parameters["b"].As<double>().PopulationMean().ConfidenceInterval(0.99).ClosedContains(b));

            Assert.IsTrue(parameters["a"].As<double>().PopulationVariance().ConfidenceInterval(0.99).ClosedContains(covariances.Column(0).Mean));
            Assert.IsTrue(parameters["b"].As<double>().PopulationVariance().ConfidenceInterval(0.99).ClosedContains(covariances.Column(1).Mean));
            Assert.IsTrue(parameters["a"].As<double>().PopulationCovariance(parameters["b"].As<double>()).ConfidenceInterval(0.99).ClosedContains(covariances.Column(2).Mean));
            Assert.IsTrue(Bivariate.PopulationCovariance(parameters["a"].As<double>(), parameters["b"].As<double>()).ConfidenceInterval(0.99).ClosedContains(covariances.Column(2).Mean));
            //Assert.IsTrue(parameters.TwoColumns(0, 1).PopulationCovariance.ConfidenceInterval(0.99).ClosedContains(covariances.Column(2).Mean));

        }

    }
}
