using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Statistics;
using Meta.Numerics.Matrices;

namespace Test
{


    [TestClass()]
    public class MultivariateSampleTest {


        private TestContext testContextInstance;

        public TestContext TestContext {
            get {
                return testContextInstance;
            }
            set {
                testContextInstance = value;
            }
        }

        #region Additional test attributes
        // 
        //You can use the following additional attributes as you write your tests:
        //
        //Use ClassInitialize to run code before running the first test in the class
        //[ClassInitialize()]
        //public static void MyClassInitialize(TestContext testContext)
        //{
        //}
        //
        //Use ClassCleanup to run code after all tests in a class have run
        //[ClassCleanup()]
        //public static void MyClassCleanup()
        //{
        //}
        //
        //Use TestInitialize to run code before running each test
        //[TestInitialize()]
        //public void MyTestInitialize()
        //{
        //}
        //
        //Use TestCleanup to run code after each test has run
        //[TestCleanup()]
        //public void MyTestCleanup()
        //{
        //}
        //
        #endregion

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
        public void MultivariateManipulationsTest () {

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
        public void MultivariateNormalTest () {

            ColumnVector V = new ColumnVector( new double[] { 1.0, 2.0} );
            SymmetricMatrix C = new SymmetricMatrix(2);
            C[0, 0] = 1.0;
            C[1, 1] = 2.0;
            C[0, 1] = 0.5;
            int N = 100;
            MultivariateSample S = CreateMultivariateNormalSample(V, C, 100);

            Assert.IsTrue(S.Count == N);

            // check the population means
            Assert.IsTrue(S.PopulationMean(0).ConfidenceInterval(0.95).ClosedContains(1.0));
            Assert.IsTrue(S.PopulationMean(1).ConfidenceInterval(0.95).ClosedContains(2.0));

            // check the population variances
            Assert.IsTrue(S.PopulationCovariance(0, 0).ConfidenceInterval(0.95).ClosedContains(C[0, 0]));
            Assert.IsTrue(S.PopulationCovariance(0, 1).ConfidenceInterval(0.95).ClosedContains(C[0, 1]));
            Assert.IsTrue(S.PopulationCovariance(1, 0).ConfidenceInterval(0.95).ClosedContains(C[1, 0]));
            Assert.IsTrue(S.PopulationCovariance(1, 1).ConfidenceInterval(0.95).ClosedContains(C[1, 1]));
            //Console.WriteLine(S.PopulationCovariance(0, 0));
            //Console.WriteLine(S.PopulationCovariance(1, 1));
            //Console.WriteLine(S.PopulationCovariance(0, 1));

            Console.WriteLine("--");
            // add tests of known higher moments for multivariate normal distribution
            // at the momement that is hard because we don't have uncertainty estimates for them
            Console.WriteLine(S.Moment(0, 0));
            Console.WriteLine(S.Mean(0));
            Console.WriteLine(S.Moment(1, 0));
            Console.WriteLine(S.Variance(0));
            Console.WriteLine(S.MomentAboutMean(2, 0));
            Console.WriteLine(S.MomentAboutMean(3, 0));
            Console.WriteLine(S.MomentAboutMean(4, 0));

        }

        [TestMethod]
        public void MultivariateNullAssociationTest () {

            Random rng = new Random(314159265);

            // Create sample sets for our three test statisics
            Sample PS = new Sample();
            Sample SS = new Sample();
            Sample KS = new Sample();

            // variables to hold the claimed distribution of teach test statistic
            Distribution PD = null;
            Distribution SD = null;
            Distribution KD = null;

            // generate a large number of multivariate samples and conduct our three tests on each

            for (int j = 0; j < 100; j++) {

                MultivariateSample S = new MultivariateSample(2);

                for (int i = 0; i < 100; i++) {
                    double x = rng.NextDouble();
                    double y = rng.NextDouble();
                    S.Add(x, y);
                }

                TestResult PR = S.PearsonR(0, 1);
                PS.Add(PR.Statistic);
                PD = PR.Distribution;
                TestResult SR = S.SpearmanRho(0, 1);
                SS.Add(SR.Statistic);
                SD = SR.Distribution;
                TestResult KR = S.KendallTau(0, 1);
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
        public void MultivariateSampleAgreementTest () {


            // check agreement of multivariate and univariate samples

            MultivariateSample ms = new MultivariateSample(1);
            Sample us = new Sample();

            // add some values

            Distribution dist = new WeibullDistribution(1.0, 2.0);
            Random rng = new Random(1);
            for (int i = 0; i < 100; i++) {
                double x = dist.InverseLeftProbability(rng.NextDouble());
                ms.Add(x);
                us.Add(x);
            }

            Console.WriteLine("count");
            Assert.IsTrue(ms.Count == us.Count);

            Console.WriteLine("descriptive");
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ms.Mean(0), us.Mean));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ms.StandardDeviation(0), us.StandardDeviation));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ms.Variance(0), us.Variance));

            Console.WriteLine("raw moments");
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ms.Moment(0), us.Moment(0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ms.Moment(1), us.Moment(1)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ms.Moment(2), us.Moment(2)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ms.Moment(3), us.Moment(3)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ms.Moment(4), us.Moment(4)));

            Console.WriteLine("central moments");
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ms.MomentAboutMean(0), us.MomentAboutMean(0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ms.MomentAboutMean(1), us.MomentAboutMean(1)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ms.MomentAboutMean(2), us.MomentAboutMean(2)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ms.MomentAboutMean(3), us.MomentAboutMean(3)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ms.MomentAboutMean(4), us.MomentAboutMean(4)));

            Console.WriteLine("population means");
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ms.PopulationMean(0).Value, us.PopulationMean.Value));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ms.PopulationMean(0).Uncertainty, us.PopulationMean.Uncertainty));

            Console.WriteLine("population variance");
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ms.PopulationCovariance(0, 0).Value, us.PopulationVariance.Value));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ms.PopulationCovariance(0, 0).Uncertainty, us.PopulationVariance.Uncertainty));

        }

    }
}
