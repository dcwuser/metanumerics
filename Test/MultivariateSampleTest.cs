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
using Meta.Numerics.Matrices;

namespace Test {

    [TestClass]
    public class MultivariateSampleTest {

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
            FrameTable table = new FrameTable();
            table.AddColumns<double>("x0", "x1", "y");

            for (int i = 0; i < 100; i++) {
                double x0 = x0distribution.GetRandomValue(rng);
                double x1 = x1distribution.GetRandomValue(rng);
                double eps = noise.GetRandomValue(rng);
                double y = a + b0 * x0 + b1 * x1 + eps;
                table.AddRow(x0, x1, y);
            }

            // do a linear regression fit on the model
            MultiLinearRegressionResult newResult = table["y"].As<double>().MultiLinearRegression(
                table["x0"].As<double>(), table["x1"].As<double>()
            );

            // the result should have the appropriate dimension
            Assert.IsTrue(newResult.Parameters.Count == 3);

            // The parameters should match the model
            Assert.IsTrue(newResult.CoefficientOf(0).ConfidenceInterval(0.99).ClosedContains(b0));
            Assert.IsTrue(newResult.CoefficientOf("x1").ConfidenceInterval(0.99).ClosedContains(b1));
            Assert.IsTrue(newResult.Intercept.ConfidenceInterval(0.99).ClosedContains(a));

            // The residuals should be compatible with the model predictions
            double ssr = 0.0;
            for (int i = 0; i < table.Rows.Count; i++) {
                FrameRow row = table.Rows[i];
                double x0 = (double) row["x0"];
                double x1 = (double) row["x1"];
                double yp = newResult.Predict(x0, x1).Value;
                double y = (double) row["y"];
                double z = y - yp;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(newResult.Residuals[i], z));
                ssr += z * z;
            }
            Assert.IsTrue(TestUtilities.IsNearlyEqual(newResult.SumOfSquaredResiduals, ssr));

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
        public void MultivariateLinearRegressionAgreement2 () {

            // A multivariate linear regression with just one x-column should be the same as a bivariate linear regression.

            double intercept = 1.0;
            double slope = -2.0;
            ContinuousDistribution yErrDist = new NormalDistribution(0.0, 3.0);
            UniformDistribution xDist = new UniformDistribution(Interval.FromEndpoints(-2.0, 3.0));
            Random rng = new Random(1111111);

            int n = 10;
            double[] x = new double[n];
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                x[i] = xDist.GetRandomValue(rng);
                y[i] = intercept + slope * x[i] + yErrDist.GetRandomValue(rng);
            }

            // Simple linear regression code.
            LinearRegressionResult result2 = y.LinearRegression(x);

            // New multi linear regression code.
            MultiLinearRegressionResult result3 = y.MultiLinearRegression(x);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(result2.Parameters.ValuesVector, result3.Parameters.ValuesVector));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(result2.Parameters.CovarianceMatrix, result3.Parameters.CovarianceMatrix));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(result2.RSquared, result3.RSquared));

        }

        [TestMethod]
        public void MultivariateRegressionRSquaredDistribution() {

            // Collect r^2 values from multivariate linear regressions.

            double cz = 1.0;
            double cx = 0.0;
            double cy = 0.0;

            Random rng = new Random(1001110000);
            ContinuousDistribution xDistribution = new UniformDistribution(Interval.FromEndpoints(-4.0, 8.0));
            ContinuousDistribution yDistribution = new UniformDistribution(Interval.FromEndpoints(-8.0, 4.0));
            ContinuousDistribution eDistribution = new NormalDistribution();

            List<double> r2Sample = new List<double>();

            for (int i = 0; i < 500; i++) {

                List<double> xs = new List<double>();
                List<double> ys = new List<double>();
                List<double> zs = new List<double>();
                for (int k = 0; k < 12; k++) {
                    double x = xDistribution.GetRandomValue(rng);
                    xs.Add(x);
                    double y = yDistribution.GetRandomValue(rng);
                    ys.Add(y);
                    double z = cx * x + cy * y + cz + eDistribution.GetRandomValue(rng);
                    zs.Add(z);
                }
                MultiLinearRegressionResult fit = zs.MultiLinearRegression(xs, ys);
                r2Sample.Add(fit.RSquared);
            }

            // r^2 values should be distributed as expected.
            ContinuousDistribution r2Distribution = new BetaDistribution((3 - 1) / 2.0, (12 - 3) / 2.0);

            TestResult ks = r2Sample.KolmogorovSmirnovTest(r2Distribution);
            Assert.IsTrue(ks.Probability > 0.05);
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
                table.AddRow(x0, x1, y);
            }

            // do a linear regression fit on the model
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
        public void MultivariateLinearRegressionNullDistribution () {

            Random rng = new Random(1);
            NormalDistribution err = new NormalDistribution();
            double[] y = new double[12];
            double[][] x = new double[4][];
            for (int k = 0; k < x.Length; k++) x[k] = new double[y.Length];

            List<double> F = new List<double>();
            ContinuousDistribution FDistribution = null;

            for (int i = 0; i < 64; i++) {
                for (int j = 0; j < y.Length; j++) {
                    y[j] = 6.0 + err.GetRandomValue(rng);
                    for (int k = 0; k < x.Length; k++) {
                        x[k][j] = rng.NextDouble();
                    }
                }
                GeneralLinearRegressionResult r = y.MultiLinearRegression(x);
                F.Add(r.F.Statistic.Value);
                FDistribution = r.F.Statistic.Distribution;
            }

            // conduct a KS test to check that F follows the expected distribution
            TestResult D = F.KolmogorovSmirnovTest(FDistribution);
            Assert.IsTrue(D.Probability > 0.05);

        }

        private double GetTotalVariance (MultivariateSample sample) {
            double total = 0.0;
            for (int i = 0; i < sample.Dimension; i++) {
                total += sample.Column(i).Variance;
            }
            return total;
        }

        private double GetTotalVariance(IReadOnlyList<IReadOnlyList<double>> sample) {
            double total = 0.0;
            for (int i = 0; i < sample.Count; i++) {
                total += sample[i].Variance();
            }
            return total;
        }

        [TestMethod]
        public void PrincipalComponentAnalysis () {

            int D = 3;
            int N = 10;

            // construct a sample
            Random rng = new Random(1);
            double[][] sample = new double[D][];
            for (int i = 0; i < D; i++) sample[i] = new double[N];
            for (int i = 0; i < N; i++) {
                sample[0][i] = 1.0 * rng.NextDouble() - 1.0;
                sample[1][i] = 4.0 * rng.NextDouble() - 2.0;
                sample[2][i] = 9.0 * rng.NextDouble() - 3.0;
            }

            // get its column means
            RowVector mu = new RowVector(D);
            for (int i = 0; i < D; i++) {
                mu[i] = sample[i].Mean();
            }

            // get total variance
            double tVariance = GetTotalVariance(sample);

            // do a principal component analysis
            PrincipalComponentAnalysis pca = sample.PrincipalComponentAnalysis();
            Assert.IsTrue(pca.Dimension == D);
            Assert.IsTrue(pca.Count == N);

            // check that the PCs behave as expected
            Assert.IsTrue(pca.Components.Count == pca.Dimension);
            for (int i = 0; i < pca.Dimension; i++) {
                PrincipalComponent pc = pca.Components[i];
                Assert.IsTrue(pc.Index == i);
                Assert.IsTrue(pc.Analysis == pca);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(pc.Weight * pc.NormalizedVector, pc.ScaledVector()));
                Assert.IsTrue(pca.MinimumDimension(pc.CumulativeVarianceFraction) == i + 1);
            }

            // Check enumerator, and verify that variance fractions behave as expected.
            int count = 0;
            double cumulative = 0.0;
            double previous = Double.PositiveInfinity;
            foreach (PrincipalComponent pc in pca.Components) {
                Assert.IsTrue(pc.Index == count);
                count++;
                Assert.IsTrue((0.0 <= pc.VarianceFraction) && (pc.VarianceFraction <= 1.0));
                Assert.IsTrue(pc.VarianceFraction <= previous);
                previous = pc.VarianceFraction;
                cumulative += pc.VarianceFraction;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(cumulative, pc.CumulativeVarianceFraction));
            }
            Assert.IsTrue(count == pca.Components.Count);

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
