using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Functions;
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
            Distribution normal = new NormalDistribution();


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

            Random rng = new Random(314159265);

            // Create sample sets for our three test statisics
            Sample PS = new Sample();
            Sample SS = new Sample();
            Sample KS = new Sample();

            // variables to hold the claimed distribution of teach test statistic
            Distribution PD = null;
            Distribution SD = null;
            Distribution KD = null;

            // generate a large number of bivariate samples and conduct our three tests on each

            for (int j = 0; j < 100; j++) {

                BivariateSample S = new BivariateSample();

                // sample size should be large so that asymptotic assumptions are justified
                for (int i = 0; i < 100; i++) {
                    double x = rng.NextDouble();
                    double y = rng.NextDouble();
                    S.Add(x, y);
                }

                TestResult PR = S.PearsonRTest();
                PS.Add(PR.Statistic);
                PD = PR.Distribution;
                TestResult SR = S.SpearmanRhoTest();
                SS.Add(SR.Statistic);
                SD = SR.Distribution;
                TestResult KR = S.KendallTauTest();
                KS.Add(KR.Statistic);
                KD = KR.Distribution;

            }

            // do KS to test whether the samples follow the claimed distributions
            //Console.WriteLine(PS.KolmogorovSmirnovTest(PD).LeftProbability);
            //Console.WriteLine(SS.KolmogorovSmirnovTest(SD).LeftProbability);
            //Console.WriteLine(KS.KolmogorovSmirnovTest(KD).LeftProbability);
            Assert.IsTrue(PS.KolmogorovSmirnovTest(PD).LeftProbability < 0.95);
            Assert.IsTrue(SS.KolmogorovSmirnovTest(SD).LeftProbability < 0.95);
            Assert.IsTrue(KS.KolmogorovSmirnovTest(KD).LeftProbability < 0.95);

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
            Console.WriteLine(s.Count);
            TestResult r = s.PairedStudentTTest();
            Console.WriteLine(r.Statistic);
            Console.WriteLine(r.LeftProbability);
        }

        [TestMethod]
        public void MultivariateLinearRegressionTest () {

            // define model y = a + b0 * x0 + b1 * x1 + noise
            double a = 1.0;
            double b0 = -2.0;
            double b1 = 3.0;
            Distribution noise = new NormalDistribution(0.0, 10.0);

            // draw a sample from the model
            Random rng = new Random(1);
            MultivariateSample sample = new MultivariateSample(3);
            for (int i = 0; i < 100; i++) {
                double x0 = -10.0 + 20.0 * rng.NextDouble();
                double x1 = -10.0 + 20.0 * rng.NextDouble();
                double eps = noise.InverseLeftProbability(rng.NextDouble());
                double y = a + b0 * x0 + b1 * x1 + eps;
                sample.Add(x0, x1, y);
            }

            // do a linear regression fit on the model
            FitResult result = sample.LinearRegression(2);

            // the result should have the appropriate dimension
            Assert.IsTrue(result.Dimension == 3);

            // the result should be significant
            Console.WriteLine("{0} {1}", result.GoodnessOfFit.Statistic, result.GoodnessOfFit.LeftProbability);
            Assert.IsTrue(result.GoodnessOfFit.LeftProbability > 0.95);

            // the parameters should match the model
            Console.WriteLine(result.Parameter(0));
            Assert.IsTrue(result.Parameter(0).ConfidenceInterval(0.90).ClosedContains(b0));
            Console.WriteLine(result.Parameter(1));
            Assert.IsTrue(result.Parameter(1).ConfidenceInterval(0.90).ClosedContains(b1));
            Console.WriteLine(result.Parameter(2));
            Assert.IsTrue(result.Parameter(2).ConfidenceInterval(0.90).ClosedContains(a));

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

        public void OldMultivariateLinearRegressionTest () {

            MultivariateSample sample = new MultivariateSample(3);
            
            sample.Add(98322, 81449, 269465);
            sample.Add(65060, 31749, 121900);
            sample.Add(36052, 14631, 37004);
            sample.Add(31829, 27732, 91400);
            sample.Add(7101, 9693, 54900);
            sample.Add(41294, 4268, 16160);
            sample.Add(16614, 4697, 21500);
            sample.Add(3449, 4233, 9306);
            sample.Add(3386, 5293, 38300);
            sample.Add(6242, 2039, 13369);
            sample.Add(14036, 7893, 29901);
            sample.Add(2636, 3345, 10930);
            sample.Add(869, 1135, 5100);
            sample.Add(452, 727, 7653);
            
            /*
            sample.Add(41.9, 29.1, 251.3);
            sample.Add(43.4, 29.3, 251.3);
            sample.Add(43.9, 29.5, 248.3);
            sample.Add(44.5, 29.7, 267.5);
            sample.Add(47.3, 29.9, 273.0);
            sample.Add(47.5, 30.3, 276.5);
            sample.Add(47.9, 30.5, 270.3);
            sample.Add(50.2, 30.7, 274.9);
            sample.Add(52.8, 30.8, 285.0);
            sample.Add(53.2, 30.9, 290.0);
            sample.Add(56.7, 31.5, 297.0);
            sample.Add(57.0, 31.7, 302.5);
            sample.Add(63.5, 31.9, 304.5);
            sample.Add(65.3, 32.0, 309.3);
            sample.Add(71.1, 32.1, 321.7);
            sample.Add(77.0, 32.5, 330.7);
            sample.Add(77.8, 32.9, 349.0);
            */

            Console.WriteLine(sample.Count);

            //sample.LinearRegression(0);
            sample.LinearRegression(0);

        }

        [TestMethod]
        public void MultivariateMoments () {

            // create a random sample
            MultivariateSample M = new MultivariateSample(3);
            Distribution d0 = new NormalDistribution();
            Distribution d1 = new ExponentialDistribution();
            Distribution d2 = new UniformDistribution();
            Random rng = new Random(1);
            int n = 10;
            for (int i = 0; i < n; i++) {
                M.Add(d0.GetRandomValue(rng), d1.GetRandomValue(rng), d2.GetRandomValue(rng));
            }

            // test that moments agree
            for (int i = 0; i < 3; i++) {
                int[] p = new int[3];
                p[i] = 1;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(M.Column(i).Mean, M.Moment(p)));
                p[i] = 2;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(M.Column(i).Variance, M.MomentAboutMean(p)));
                for (int j = 0; j < i; j++) {
                    int[] q = new int[3];
                    q[i] = 1;
                    q[j] = 1;
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(M.TwoColumns(i, j).Covariance, M.MomentAboutMean(q)));
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
                FitResult r = ms.LinearRegression(0);
                fs.Add(r.GoodnessOfFit.Statistic);
            }

            // conduct a KS test to check that F follows the expected distribution
            TestResult ks = fs.KolmogorovSmirnovTest(new FisherDistribution(3, 4));
            Assert.IsTrue(ks.LeftProbability < 0.95);

        }

        [TestMethod]
        public void MultivariateLinearRegressionAgreement () {

            Random rng = new Random(1);

            MultivariateSample SA = new MultivariateSample(2);
            for (int i = 0; i < 10; i++) {
                SA.Add(rng.NextDouble(), rng.NextDouble());
            }
            FitResult RA = SA.LinearRegression(0);
            double[] PA = RA.Parameters();
            SymmetricMatrix CA = RA.CovarianceMatrix();

            MultivariateSample SB = SA.Columns(1, 0);
            FitResult RB = SB.LinearRegression(1);
            double[] PB = RB.Parameters();
            SymmetricMatrix CB = RB.CovarianceMatrix();

            Assert.IsTrue(TestUtilities.IsNearlyEqual(PA[0], PB[1])); Assert.IsTrue(TestUtilities.IsNearlyEqual(PA[1], PB[0]));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(CA[0, 0], CB[1, 1])); Assert.IsTrue(TestUtilities.IsNearlyEqual(CA[0, 1], CB[1, 0])); Assert.IsTrue(TestUtilities.IsNearlyEqual(CA[1, 1], CB[0, 0]));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(RA.GoodnessOfFit.Statistic, RB.GoodnessOfFit.Statistic));

            BivariateSample SC = SA.TwoColumns(1, 0);
            FitResult RC = SC.LinearRegression();
            double[] PC = RC.Parameters();
            SymmetricMatrix CC = RC.CovarianceMatrix();

            Assert.IsTrue(TestUtilities.IsNearlyEqual(PA, PC));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(CA, CC));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(RA.GoodnessOfFit.Statistic, RC.GoodnessOfFit.Statistic));

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
                PrincipalComponent pc = pca.Component(i);
                Assert.IsTrue(pc.Index == i);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(pc.Weight * pc.NormalizedVector(), pc.ScaledVector()));
                Assert.IsTrue((0.0 <= pc.VarianceFraction) && (pc.VarianceFraction <= 1.0));
                if (i == 0) {
                    Assert.IsTrue(pc.VarianceFraction == pc.CumulativeVarianceFraction);
                } else {
                    PrincipalComponent ppc = pca.Component(i - 1);
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
                        PrincipalComponent pc = pca.Component(i);
                        x += (cEntry[i] * pc.Weight) * pc.NormalizedVector();
                    }
                    rSample.Add(x);
                }
                double rVariance = GetTotalVariance(rSample);
                Console.WriteLine("{0} {1}", rD, rVariance);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(rVariance / tVariance, pca.Component(rD-1).CumulativeVarianceFraction));
            }

        }

    }
}
