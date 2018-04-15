using System;
using System.Collections.Generic;
using System.Linq;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Data;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;
using Meta.Numerics.Matrices;

namespace Test {

    [TestClass]
    public class MultivariateSampleTest {

        public MultivariateSample CreateMultivariateNormalSample (ColumnVector M, SymmetricMatrix C, int n) {

            int d = M.Dimension;

            MultivariateSample S = new MultivariateSample(d);

            SquareMatrix A = C.CholeskyDecomposition().SquareRootMatrix();

            Random rng = new Random(1);
            ContinuousDistribution normal = new NormalDistribution();


            for (int i = 0; i < n; i++) {

                // create a vector of normal deviates
                ColumnVector V = new ColumnVector(d);
                for (int j = 0; j < d; j++) {
                    double y = rng.NextDouble();
                    double z = normal.InverseLeftProbability(y);
                    V[j] = z;
                }

                // form the multivariate distributed vector
                ColumnVector X = M + A * V;

                // add it to the sample
                S.Add(X);


            }

            return (S);

        }

        [TestMethod]
        public void MultivariateManipulations () {

            MultivariateSample S = new MultivariateSample(3);

            Assert.IsTrue(S.Dimension == 3);

            Assert.IsTrue(S.Count == 0);

            S.Add(1.1, 1.2, 1.3);
            S.Add(2.1, 2.2, 2.3);

            Assert.IsTrue(S.Count == 2);

            // check that an entry is there, remove it, check that it is not there
            Assert.IsTrue(S.Contains(1.1, 1.2, 1.3));
            Assert.IsTrue(S.Remove(1.1, 1.2, 1.3));
            Assert.IsFalse(S.Contains(1.1, 1.2, 1.3));

            // clear it and check that the count went to zero
            S.Clear();
            Assert.IsTrue(S.Count == 0);

        }

        [TestMethod]
        public void MultivariateNormalSummaryStatistics () {

            ColumnVector V = new ColumnVector( new double[] { 1.0, 2.0} );
            SymmetricMatrix C = new SymmetricMatrix(2);
            C[0, 0] = 1.0;
            C[1, 1] = 2.0;
            C[0, 1] = 0.5;
            int N = 100;
            MultivariateSample S = CreateMultivariateNormalSample(V, C, 100);

            Assert.IsTrue(S.Count == N);

            // check the population means
            Assert.IsTrue(S.Column(0).PopulationMean.ConfidenceInterval(0.95).ClosedContains(1.0));
            Assert.IsTrue(S.Column(1).PopulationMean.ConfidenceInterval(0.95).ClosedContains(2.0));

            // check the population variances
            Assert.IsTrue(S.Column(0).PopulationVariance.ConfidenceInterval(0.95).ClosedContains(C[0, 0]));
            //Assert.IsTrue(S.PopulationCovariance(0, 1).ConfidenceInterval(0.95).ClosedContains(C[0, 1]));
            //Assert.IsTrue(S.PopulationCovariance(1, 0).ConfidenceInterval(0.95).ClosedContains(C[1, 0]));
            Assert.IsTrue(S.Column(1).PopulationVariance.ConfidenceInterval(0.95).ClosedContains(C[1, 1]));
            //Console.WriteLine(S.PopulationCovariance(0, 0));
            //Console.WriteLine(S.PopulationCovariance(1, 1));
            //Console.WriteLine(S.PopulationCovariance(0, 1));

            Console.WriteLine("--");
            // add tests of known higher moments for multivariate normal distribution
            // at the momement that is hard because we don't have uncertainty estimates for them
            //Console.WriteLine(S.Moment(0, 0));
            //Console.WriteLine(S.Mean(0));
            //Console.WriteLine(S.Moment(1, 0));
            //Console.WriteLine(S.Variance(0));
            //Console.WriteLine(S.MomentAboutMean(2, 0));
            //Console.WriteLine(S.MomentAboutMean(3, 0));
            //Console.WriteLine(S.MomentAboutMean(4, 0));

        }

        [TestMethod]
        public void BivariateNullAssociation () {

            Random rng = new Random(31415926);

            // Create a data structure to hold the results of Pearson, Spearman, and Kendall tests.
            FrameTable data = new FrameTable();
            data.AddColumn<double>("r");
            data.AddColumn<double>("ρ");
            data.AddColumn<double>("τ");

            // Create variables to hold the claimed distribution of each test statistic.
            ContinuousDistribution PRD = null;
            ContinuousDistribution SRD = null;
            ContinuousDistribution KTD = null;

            // Generate a large number of bivariate samples and conduct our three tests on each.
            ContinuousDistribution xDistribution = new LognormalDistribution();
            ContinuousDistribution yDistribution = new CauchyDistribution();
            for (int j = 0; j < 100; j++) {

                List<double> x = new List<double>();
                List<double> y = new List<double>();
                for (int i = 0; i < 100; i++) {
                    x.Add(xDistribution.GetRandomValue(rng));
                    y.Add(yDistribution.GetRandomValue(rng));
                }

                TestResult PR = Bivariate.PearsonRTest(x, y);
                TestResult SR = Bivariate.SpearmanRhoTest(x, y);
                TestResult KT = Bivariate.KendallTauTest(x, y);

                PRD = PR.Statistic.Distribution;
                SRD = SR.Statistic.Distribution;
                KTD = KT.Statistic.Distribution;

                data.AddRow(new Dictionary<string, object>() {
                    {"r", PR.Statistic.Value}, {"ρ", SR.Statistic.Value}, {"τ", KT.Statistic.Value}
                });
            }

            Assert.IsTrue(data["r"].As<double>().KolmogorovSmirnovTest(PRD).Probability > 0.05);
            Assert.IsTrue(data["ρ"].As<double>().KolmogorovSmirnovTest(SRD).Probability > 0.05);
            Assert.IsTrue(data["τ"].As<double>().KolmogorovSmirnovTest(KTD).Probability > 0.05);
        }
        
        [TestMethod]
        public void PairedStudentTTest () {

            BivariateSample s = new BivariateSample();
            s.Add(3, 5);
            s.Add(0, 1);
            s.Add(6, 5);
            s.Add(7, 7);
            s.Add(4, 10);
            s.Add(3, 9);
            s.Add(2, 7);
            s.Add(1, 11);
            s.Add(4, 8);
            TestResult r = s.PairedStudentTTest();
            // Maybe we should assert something here?
        }

        [TestMethod]
        public void MultivariateLinearRegressionSimple () {

            // define model y = a + b0 * x0 + b1 * x1 + noise
            double a = 1.0;
            double b0 = -2.0;
            double b1 = 3.0;
            ContinuousDistribution x0distribution = new CauchyDistribution(10.0, 5.0);
            ContinuousDistribution x1distribution = new UniformDistribution(Interval.FromEndpoints(-10.0, 20.0));
            ContinuousDistribution noise = new NormalDistribution(0.0, 10.0);

            // draw a sample from the model
            Random rng = new Random(1);
            MultivariateSample sample = new MultivariateSample("x0", "x1", "y");
            FrameTable table = new FrameTable();
            table.AddColumns<double>("x0", "x1", "y");

            for (int i = 0; i < 100; i++) {
                double x0 = x0distribution.GetRandomValue(rng);
                double x1 = x1distribution.GetRandomValue(rng);
                double eps = noise.GetRandomValue(rng);
                double y = a + b0 * x0 + b1 * x1 + eps;
                sample.Add(x0, x1, y);
                table.AddRow(x0, x1, y);
            }

            // do a linear regression fit on the model
            ParameterCollection oldResult = sample.LinearRegression(2).Parameters;
            MultiLinearRegressionResult newResult = table["y"].As<double>().MultiLinearRegression(
                table["x0"].As<double>(), table["x1"].As<double>()
            );

            // the result should have the appropriate dimension
            Assert.IsTrue(oldResult.Count == 3);
            Assert.IsTrue(newResult.Parameters.Count == 3);

            // The parameters should match the model
            Assert.IsTrue(oldResult[0].Estimate.ConfidenceInterval(0.90).ClosedContains(b0));
            Assert.IsTrue(oldResult[1].Estimate.ConfidenceInterval(0.90).ClosedContains(b1));
            Assert.IsTrue(oldResult[2].Estimate.ConfidenceInterval(0.90).ClosedContains(a));

            Assert.IsTrue(newResult.CoefficientOf(0).ConfidenceInterval(0.99).ClosedContains(b0));
            Assert.IsTrue(newResult.CoefficientOf("x1").ConfidenceInterval(0.99).ClosedContains(b1));
            Assert.IsTrue(newResult.Intercept.ConfidenceInterval(0.99).ClosedContains(a));

            // The residuals should be compatible with the model predictions
            for (int i = 0; i < table.Rows.Count; i++) {
                FrameRow row = table.Rows[i];
                double x0 = (double) row["x0"];
                double x1 = (double) row["x1"];
                double yp = newResult.Predict(x0, x1).Value;
                double y = (double) row["y"];
                Assert.IsTrue(TestUtilities.IsNearlyEqual(newResult.Residuals[i], y - yp));
            }

        }


        [TestMethod]
        public void MultivariateLinearRegressionVariances () {

            // define model y = a + b0 * x0 + b1 * x1 + noise
            double a = -3.0;
            double b0 = 2.0;
            double b1 = -1.0;
            ContinuousDistribution x0distribution = new LaplaceDistribution();
            ContinuousDistribution x1distribution = new CauchyDistribution();
            ContinuousDistribution eDistribution = new NormalDistribution(0.0, 4.0);

            FrameTable data = new FrameTable();
            data.AddColumns<double>("a", "da", "b0", "db0", "b1", "db1", "ab1Cov", "p", "dp");

            // draw a sample from the model
            Random rng = new Random(4);
            for (int j = 0; j < 64; j++) {
                List<double> x0s = new List<double>();
                List<double> x1s = new List<double>();
                List<double> ys = new List<double>();

                for (int i = 0; i < 16; i++) {
                    double x0 = x0distribution.GetRandomValue(rng);
                    double x1 = x1distribution.GetRandomValue(rng);
                    double e = eDistribution.GetRandomValue(rng);
                    double y = a + b0 * x0 + b1 * x1 + e;
                    x0s.Add(x0);
                    x1s.Add(x1);
                    ys.Add(y);
                }

                // do a linear regression fit on the model
                MultiLinearRegressionResult result = ys.MultiLinearRegression(
                    new Dictionary<string, IReadOnlyList<double>> {
                        {"x0", x0s }, {"x1", x1s }
                    }
                );
                UncertainValue pp = result.Predict(-5.0, 6.0);

                data.AddRow(
                    result.Intercept.Value, result.Intercept.Uncertainty,
                    result.CoefficientOf("x0").Value, result.CoefficientOf("x0").Uncertainty,
                    result.CoefficientOf("x1").Value, result.CoefficientOf("x1").Uncertainty,
                    result.Parameters.CovarianceOf("Intercept", "x1"),
                    pp.Value, pp.Uncertainty
                );

            }

            // The estimated parameters should agree with the model that generated the data.

            // The variances of the estimates should agree with the claimed variances
            Assert.IsTrue(data["a"].As<double>().PopulationStandardDeviation().ConfidenceInterval(0.99).ClosedContains(data["da"].As<double>().Mean()));
            Assert.IsTrue(data["b0"].As<double>().PopulationStandardDeviation().ConfidenceInterval(0.99).ClosedContains(data["db0"].As<double>().Mean()));
            Assert.IsTrue(data["b1"].As<double>().PopulationStandardDeviation().ConfidenceInterval(0.99).ClosedContains(data["db1"].As<double>().Mean()));
            Assert.IsTrue(data["a"].As<double>().PopulationCovariance(data["b1"].As<double>()).ConfidenceInterval(0.99).ClosedContains(data["ab1Cov"].As<double>().Mean()));
            Assert.IsTrue(data["p"].As<double>().PopulationStandardDeviation().ConfidenceInterval(0.99).ClosedContains(data["dp"].As<double>().Median()));

        }

        [TestMethod]
        public void MultivariateLinearRegressionBadInputTest () {

            // create a sample
            MultivariateSample sample = new MultivariateSample(3);
            sample.Add(1, 2, 3);
            sample.Add(2, 3, 4);
            
            // try to predict with too little data
            try {
                sample.LinearRegression(2);
                Assert.IsTrue(false);
            } catch (InvalidOperationException) {
                Assert.IsTrue(true);
            }

            // add enough data
            sample.Add(3, 4, 5);
            sample.Add(4, 5, 6);

            // try to predict a non-existent variable
            try {
                sample.LinearRegression(-1);
                Assert.IsTrue(false);
            } catch (ArgumentOutOfRangeException) {
                Assert.IsTrue(true);
            }

            try {
                sample.LinearRegression(3);
                Assert.IsTrue(false);
            } catch (ArgumentOutOfRangeException) {
                Assert.IsTrue(true);
            }

        }

        [TestMethod]
        public void MultivariateLinearRegressionAgreement2 () {

            // A multivariate linear regression with just one x-column should be the same as a bivariate linear regression.

            double intercept = 1.0;
            double slope = -2.0;
            ContinuousDistribution yErrDist = new NormalDistribution(0.0, 3.0);
            UniformDistribution xDist = new UniformDistribution(Interval.FromEndpoints(-2.0, 3.0));
            Random rng = new Random(1111111);

            MultivariateSample multi = new MultivariateSample("x", "y");
            for (int i = 0; i < 10; i++) {
                double x = xDist.GetRandomValue(rng);
                double y = intercept + slope * x + yErrDist.GetRandomValue(rng);
                multi.Add(x, y);
            }

            // Old multi linear regression code.
            MultiLinearRegressionResult result1 = multi.LinearRegression(1);

            // Simple linear regression code.
            LinearRegressionResult result2 = multi.TwoColumns(0, 1).LinearRegression();
            Assert.IsTrue(TestUtilities.IsNearlyEqual(result1.Parameters["Intercept"].Estimate, result2.Parameters["Intercept"].Estimate));

            // New multi linear regression code.
            MultiLinearRegressionResult result3 = multi.Column(1).ToList().MultiLinearRegression(multi.Column(0).ToList());
            Assert.IsTrue(TestUtilities.IsNearlyEqual(result1.Parameters["Intercept"].Estimate, result3.Parameters["Intercept"].Estimate));

        }


        [TestMethod]
        public void MultivariateLinearLogisticRegressionSimple () {

            // define model y = a + b0 * x0 + b1 * x1 + noise
            double a = 1.0;
            double b0 = -1.0 / 2.0;
            double b1 = 1.0 / 3.0;
            ContinuousDistribution x0distribution = new LaplaceDistribution();
            ContinuousDistribution x1distribution = new NormalDistribution();

            // draw a sample from the model
            Random rng = new Random(1);
            MultivariateSample old = new MultivariateSample("y", "x0", "x1");
            FrameTable table = new FrameTable();
            table.AddColumn<double>("x0");
            table.AddColumn<double>("x1");
            table.AddColumn<bool>("y");
            
            for (int i = 0; i < 100; i++) {
                double x0 = x0distribution.GetRandomValue(rng);
                double x1 = x1distribution.GetRandomValue(rng);
                double t = a + b0 * x0 + b1 * x1;
                double p = 1.0 / (1.0 + Math.Exp(-t));
                bool y = (rng.NextDouble() < p);
                old.Add(y ? 1.0 : 0.0, x0, x1);
                table.AddRow(x0, x1, y);
            }

            // do a linear regression fit on the model
            MultiLinearLogisticRegressionResult oldResult = old.LogisticLinearRegression(0);
            MultiLinearLogisticRegressionResult newResult = table["y"].As<bool>().MultiLinearLogisticRegression(
                table["x0"].As<double>(), table["x1"].As<double>()
            );

            // the result should have the appropriate dimension
            Assert.IsTrue(newResult.Parameters.Count == 3);

            // The parameters should match the model
            Assert.IsTrue(newResult.CoefficientOf(0).ConfidenceInterval(0.99).ClosedContains(b0));
            Assert.IsTrue(newResult.CoefficientOf("x1").ConfidenceInterval(0.99).ClosedContains(b1));
            Assert.IsTrue(newResult.Intercept.ConfidenceInterval(0.99).ClosedContains(a));

            // Our predictions should be better than chance.
            int correct = 0;
            for (int i = 0; i < table.Rows.Count; i++) {
                FrameRow row = table.Rows[i];
                double x0 = (double) row["x0"];
                double x1 = (double) row["x1"];
                double p = newResult.Predict(x0, x1).Value;
                bool y = (bool) row["y"];
                if ((y && p > 0.5) || (!y & p < 0.5)) correct++;
            }
            Assert.IsTrue(correct > 0.5 * table.Rows.Count);

        }

        [TestMethod]
        public void MultivariateLinearLogisticRegressionVariances () {

            // define model y = a + b0 * x0 + b1 * x1 + noise
            double a = -3.0;
            double b0 = 2.0;
            double b1 = 1.0;
            ContinuousDistribution x0distribution = new ExponentialDistribution();
            ContinuousDistribution x1distribution = new LognormalDistribution();

            FrameTable data = new FrameTable();
            data.AddColumns<double>("a", "da", "b0", "db0", "b1", "db1", "p", "dp" );

            // draw a sample from the model
            Random rng = new Random(2);
            for (int j = 0; j < 32; j++) {
                List<double> x0s = new List<double>();
                List<double> x1s = new List<double>();
                List<bool> ys = new List<bool>();

                FrameTable table = new FrameTable();
                table.AddColumn<double>("x0");
                table.AddColumn<double>("x1");
                table.AddColumn<bool>("y");
                
                for (int i = 0; i < 32; i++) {
                    double x0 = x0distribution.GetRandomValue(rng);
                    double x1 = x1distribution.GetRandomValue(rng);
                    double t = a + b0 * x0 + b1 * x1;
                    double p = 1.0 / (1.0 + Math.Exp(-t));
                    bool y = (rng.NextDouble() < p);
                    x0s.Add(x0);
                    x1s.Add(x1);
                    ys.Add(y);
                }

                // do a linear regression fit on the model
                MultiLinearLogisticRegressionResult result = ys.MultiLinearLogisticRegression(
                    new Dictionary<string, IReadOnlyList<double>> {
                        {"x0", x0s }, {"x1", x1s }
                    }
                );
                UncertainValue pp = result.Predict(0.0, 1.0);

                data.AddRow(
                    result.Intercept.Value, result.Intercept.Uncertainty,
                    result.CoefficientOf("x0").Value, result.CoefficientOf("x0").Uncertainty,
                    result.CoefficientOf("x1").Value, result.CoefficientOf("x1").Uncertainty,
                    pp.Value, pp.Uncertainty
                );

            }

            // The estimated parameters should agree with the model that generated the data.

            // The variances of the estimates should agree with the claimed variances
            Assert.IsTrue(data["a"].As<double>().PopulationStandardDeviation().ConfidenceInterval(0.99).ClosedContains(data["da"].As<double>().Mean()));
            Assert.IsTrue(data["b0"].As<double>().PopulationStandardDeviation().ConfidenceInterval(0.99).ClosedContains(data["db0"].As<double>().Mean()));
            Assert.IsTrue(data["b1"].As<double>().PopulationStandardDeviation().ConfidenceInterval(0.99).ClosedContains(data["db1"].As<double>().Mean()));
            Assert.IsTrue(data["p"].As<double>().PopulationStandardDeviation().ConfidenceInterval(0.99).ClosedContains(data["dp"].As<double>().Mean()));

        }

        [TestMethod]
        public void MultivariateMoments () {

            // create a random sample
            MultivariateSample M = new MultivariateSample(3);
            ContinuousDistribution d0 = new NormalDistribution();
            ContinuousDistribution d1 = new ExponentialDistribution();
            ContinuousDistribution d2 = new UniformDistribution();
            Random rng = new Random(1);
            int n = 10;
            for (int i = 0; i < n; i++) {
                M.Add(d0.GetRandomValue(rng), d1.GetRandomValue(rng), d2.GetRandomValue(rng));
            }

            // test that moments agree
            for (int i = 0; i < 3; i++) {
                int[] p = new int[3];
                p[i] = 1;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(M.Column(i).Mean, M.RawMoment(p)));
                p[i] = 2;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(M.Column(i).Variance, M.CentralMoment(p)));
                for (int j = 0; j < i; j++) {
                    int[] q = new int[3];
                    q[i] = 1;
                    q[j] = 1;
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(M.TwoColumns(i, j).Covariance, M.CentralMoment(q)));
                }
            }


        }

        [TestMethod]
        public void MultivariateLinearRegressionNullDistribution () {

            int d = 4;

            Random rng = new Random(1);
            NormalDistribution n = new NormalDistribution();

            Sample fs = new Sample();

            for (int i = 0; i < 64; i++) {
                MultivariateSample ms = new MultivariateSample(d);
                for (int j = 0; j < 8; j++) {
                    double[] x = new double[d];
                    for (int k = 0; k < d; k++) {
                        x[k] = n.GetRandomValue(rng);
                    }
                    ms.Add(x);
                }
                GeneralLinearRegressionResult r = ms.LinearRegression(0);
                fs.Add(r.F.Statistic.Value);
            }

            // conduct a KS test to check that F follows the expected distribution
            TestResult ks = fs.KolmogorovSmirnovTest(new FisherDistribution(3, 4));
            Assert.IsTrue(ks.Probability > 0.05);

        }

        [TestMethod]
        public void MultivariateLinearRegressionAgreement () {

            Random rng = new Random(1);

            MultivariateSample SA = new MultivariateSample(2);
            for (int i = 0; i < 10; i++) {
                SA.Add(rng.NextDouble(), rng.NextDouble());
            }
            GeneralLinearRegressionResult RA = SA.LinearRegression(0);
            ColumnVector PA = RA.Parameters.ValuesVector;
            SymmetricMatrix CA = RA.Parameters.CovarianceMatrix;

            MultivariateSample SB = SA.Columns(1, 0);
            GeneralLinearRegressionResult RB = SB.LinearRegression(1);
            ColumnVector PB = RB.Parameters.ValuesVector;
            SymmetricMatrix CB = RB.Parameters.CovarianceMatrix;

            Assert.IsTrue(TestUtilities.IsNearlyEqual(PA[0], PB[1])); Assert.IsTrue(TestUtilities.IsNearlyEqual(PA[1], PB[0]));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(CA[0, 0], CB[1, 1])); Assert.IsTrue(TestUtilities.IsNearlyEqual(CA[0, 1], CB[1, 0])); Assert.IsTrue(TestUtilities.IsNearlyEqual(CA[1, 1], CB[0, 0]));

            BivariateSample SC = SA.TwoColumns(1, 0);
            GeneralLinearRegressionResult RC = SC.LinearRegression();
            ColumnVector PC = RC.Parameters.ValuesVector;
            SymmetricMatrix CC = RC.Parameters.CovarianceMatrix;

            Assert.IsTrue(TestUtilities.IsNearlyEqual(PA, PC));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(CA, CC));

        }

        private double GetTotalVariance (MultivariateSample sample) {
            double total = 0.0;
            for (int i = 0; i < sample.Dimension; i++) {
                total += sample.Column(i).Variance;
            }
            return (total);
        }

        [TestMethod]
        public void PrincipalComponentAnalysis () {

            int D = 3;
            int N = 10;

            // construct a sample
            Random rng = new Random(1);
            MultivariateSample sample = new MultivariateSample(D);
            for (int i = 0; i < N; i++) {
                double x = 1.0 * rng.NextDouble() - 1.0;
                double y = 4.0 * rng.NextDouble() - 2.0;
                double z = 9.0 * rng.NextDouble() - 3.0;
                sample.Add(x, y, z);
            }

            // get its column means
            RowVector mu = new RowVector(D);
            for (int i = 0; i < D; i++) {
                mu[i] = sample.Column(i).Mean;
            }

            // get total variance
            double tVariance = GetTotalVariance(sample);
            Console.WriteLine(tVariance);

            // do a principal component analysis
            PrincipalComponentAnalysis pca = sample.PrincipalComponentAnalysis();
            Assert.IsTrue(pca.Dimension == sample.Dimension);
            Assert.IsTrue(pca.Count == sample.Count);

            // check that the PCs behave as expected
            for (int i = 0; i < pca.Dimension; i++) {
                PrincipalComponent pc = pca.Components[i];
                Assert.IsTrue(pc.Index == i);
                Assert.IsTrue(pc.Analysis == pca);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(pc.Weight * pc.NormalizedVector, pc.ScaledVector()));
                Assert.IsTrue((0.0 <= pc.VarianceFraction) && (pc.VarianceFraction <= 1.0));
                if (i == 0) {
                    Assert.IsTrue(pc.VarianceFraction == pc.CumulativeVarianceFraction);
                } else {
                    PrincipalComponent ppc = pca.Components[i - 1];
                    Assert.IsTrue(pc.VarianceFraction <= ppc.VarianceFraction);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(ppc.CumulativeVarianceFraction + pc.VarianceFraction, pc.CumulativeVarianceFraction));
                }
            }

            // express the sample in terms of principal components
            MultivariateSample csample = pca.TransformedSample();

            // check that the explained variances are as claimed
            for (int rD = 1; rD <= D; rD++) {
                MultivariateSample rSample = new MultivariateSample(D);
                foreach (double[] cEntry in csample) {
                    RowVector x = mu.Copy();
                    for (int i = 0; i < rD; i++) {
                        PrincipalComponent pc = pca.Components[i];
                        x += (cEntry[i] * pc.Weight) * pc.NormalizedVector;
                    }
                    rSample.Add(x);
                }
                double rVariance = GetTotalVariance(rSample);
                Console.WriteLine("{0} {1}", rD, rVariance);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(rVariance / tVariance, pca.Components[rD-1].CumulativeVarianceFraction));
            }

        }

        [TestMethod]
        public void MeansClustering2 () {

            ColumnVector[] centers = new ColumnVector[] {
                new ColumnVector(0.0, 0.0, 0.0),
                new ColumnVector(2.0, 0.0, 0.0),
                new ColumnVector(0.0, 2.0, 0.0),
                new ColumnVector(0.0, 0.0, 2.0)
            };

            FrameTable table = new FrameTable();
            string alphabet = "abcdefghijklmnopqrstuvwxyz";
            for(int j = 0; j < 3; j++) {
                table.AddColumn<double>(alphabet[j].ToString());
            }
            
             
            List<int> inputAssignments = new List<int>();
            List<ColumnVector> inputVectors = new List<ColumnVector>();
            Random rng = new Random(2);
            ContinuousDistribution dist = new NormalDistribution(0.0, 1.0);
            for (int i = 0; i < 100; i++) {
                int inputAssignment = rng.Next(0, centers.Length);
                inputAssignments.Add(inputAssignment);
                ColumnVector inputVector = centers[inputAssignment].Copy();
                for (int k = 0; k < inputVector.Dimension; k++) {
                    inputVector[k] += dist.GetRandomValue(rng);
                }
                inputVectors.Add(inputVector);
                table.AddRow<double>(inputVector);
            }

            MeansClusteringResult result = table.AsColumns<double>().MeansClustering(centers.Length);

            //MultivariateSample s = new MultivariateSample(3);
            //foreach (ColumnVector v in inputVectors) { s.Add(v); }
            //MeansClusteringResult result = s.MeansClustering(centers.Length);

            List<int> outputAssignments = new List<int>();
            for (int i = 0; i < inputVectors.Count; i++) {
                int assignment = result.Classify(inputVectors[i]);
                outputAssignments.Add(assignment);
            }

            // Map the output centroids to the original centroids
            Dictionary<int, int> map = new Dictionary<int, int>();
            for (int outputIndex = 0; outputIndex < result.Count; outputIndex++) {
                ColumnVector centroid = result.Centroid(outputIndex);
                int mappedInputIndex = -1;
                double mappedInputDistance = Double.MaxValue;
                for (int inputIndex = 0; inputIndex < centers.Length; inputIndex++) {
                    double distance = (centroid - centers[inputIndex]).Norm();
                    if (distance < mappedInputDistance) {
                        mappedInputIndex = inputIndex;
                        mappedInputDistance = distance;
                    }
                }
                Assert.IsTrue(mappedInputIndex >= 0);
                Assert.IsTrue(mappedInputDistance < 1.0);
                map.Add(outputIndex, mappedInputIndex);
            }

            int correctCount = 0;
            for (int i = 0; i < outputAssignments.Count; i++) {
                if (map[outputAssignments[i]] == inputAssignments[i]) correctCount++;
            }
            Assert.IsTrue(correctCount >= 0.50 * outputAssignments.Count);

        }

        [TestMethod]
        public void MeansClustering () {

            // Re-create the mouse test

            double[] x = new double[3];
            double[] y = new double[3];
            double[] s = new double[3];

            x[0] = 0.25;
            y[0] = 0.75;
            s[0] = 0.1;

            x[1] = 0.75;
            y[1] = 0.75;
            s[1] = 0.1;

            x[2] = 0.5;
            y[2] = 0.5;
            s[2] = 0.2;

            MultivariateSample points = new MultivariateSample(2);
            Random rng = new Random(1);
            NormalDistribution d = new NormalDistribution();
            for (int i = 0; i < 100; i++) {
                int k = rng.Next(3);
                points.Add(x[k] + s[k] * d.GetRandomValue(rng), y[k] + s[k] * d.GetRandomValue(rng));
            }

            MeansClusteringResult result = points.MeansClustering(3);

            Assert.IsTrue(result.Count == 3);
            Assert.IsTrue(result.Dimension == 2);

        }


    }

}
