using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;


namespace Test {
    
    [TestClass()]
    public class DiscreteDistributionTest {


        private TestContext testContextInstance;

        /// <summary>
        ///Gets or sets the test context which provides
        ///information about and functionality for the current test run.
        ///</summary>
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


        private DiscreteDistribution[] distributions = new DiscreteDistribution[] {
            new BernoulliDistribution(0.1),
            new BinomialDistribution(0.2, 30), new BinomialDistribution(0.4, 5),
            new PoissonDistribution(4.5),
            new DiscreteUniformDistribution(5, 11),
            new GeometricDistribution(0.6)
        };

        [TestMethod]
        public void DiscreteDistributionMomentSpecialCases () {

            foreach (DiscreteDistribution distribution in distributions) {
                Assert.IsTrue(distribution.Moment(0) == 1.0);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(distribution.Moment(1), distribution.Mean));
                Assert.IsTrue(distribution.MomentAboutMean(0) == 1.0);
                Assert.IsTrue(distribution.MomentAboutMean(1) == 0.0);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(distribution.MomentAboutMean(2), distribution.Variance));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(Math.Sqrt(distribution.Variance), distribution.StandardDeviation));
            }

        }

        [TestMethod]
        public void DiscreteDistributionMomentSums () {
            foreach (DiscreteDistribution distribution in distributions) {
                Console.WriteLine(distribution.GetType().Name);
                // C2 = M2 - M1^2
                double M1 = distribution.Moment(1);
                double M2 = distribution.Moment(2);
                double C2 = distribution.MomentAboutMean(2);
                if (!Double.IsInfinity(C2)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(C2 + M1 * M1, M2));
                }
                // C3 = M3 - 3 M2 M1 + 2 M1^3
                double M3 = distribution.Moment(3);
                double C3 = distribution.MomentAboutMean(3);
                if (!Double.IsInfinity(C3)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(C3 + 3.0 * M2 * M1, M3 + 2.0 * M1 * M1 * M1));
                }
                // C4 = M4 - 4 M3 M1 + 6 M2 M1^2 - 3 M1^4
            }
        }

        [TestMethod]
        public void DiscreteDistributionUnitarity () {
            foreach (DiscreteDistribution distribution in distributions) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    distribution.ExpectationValue(delegate (int k) { return (1.0); }), 1.0
                ));
            }
        }

        [TestMethod]
        public void DiscreteDistributionMean () {
            foreach (DiscreteDistribution distribution in distributions) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    distribution.ExpectationValue(delegate(int k) { return (k); }), distribution.Mean
                ));
            }
        }

        [TestMethod]
        public void DiscreteDistributionVariance () {
            foreach (DiscreteDistribution distribution in distributions) {
                double m = distribution.Mean;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    distribution.ExpectationValue(delegate(int x) { return (Math.Pow(x-m, 2)); }), distribution.Variance
                ));
            }
        }

        [TestMethod]
        public void DiscreteDistributionSkewness () {
            foreach (DiscreteDistribution distribution in distributions) {
                Console.WriteLine(distribution.GetType().FullName);
                //Console.WriteLine(distribution.MomentAboutMean(3));
                //Console.WriteLine(distribution.MomentAboutMean(2));
                //Console.WriteLine(distribution.Skewness);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    distribution.Skewness, distribution.MomentAboutMean(3) / Math.Pow(distribution.MomentAboutMean(2), 3.0 / 2.0)
                ));
            }
        }

        [TestMethod]
        public void DiscreteDistributionProbabilityAxioms () {

            foreach (DiscreteDistribution distribution in distributions) {

                // some of these values will be outside the support, but that's fine, our results should still be consistent with probability axioms
                foreach (int k in TestUtilities.GenerateUniformIntegerValues(-10, +100, 6)) {

                    double DP = distribution.ProbabilityMass(k);
                    Assert.IsTrue(DP >= 0.0); Assert.IsTrue(DP <= 1.0);

                    double P = distribution.LeftProbability(k);
                    double Q = distribution.RightProbability(k);

                    Assert.IsTrue(P >= 0.0); Assert.IsTrue(P <= 1.0);
                    Assert.IsTrue(Q >= 0.0); Assert.IsTrue(Q <= 1.0);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(P + Q, 1.0));

                }

            }

        }

        [TestMethod]
        public void DiscreteDistributionInverseCDF () {

            Random rng = new Random(1);
            for (int i = 0; i < 10; i++) {

                double P = rng.NextDouble();

                foreach (DiscreteDistribution distribution in distributions) {
                    int x = distribution.InverseLeftProbability(P);
                    Console.WriteLine("{0} {1} {2} {3}", distribution.GetType().Name, P, x, distribution.LeftProbability(x));
                    Assert.IsTrue(distribution.LeftProbability(x - 1) < P);
                    Assert.IsTrue(P <= distribution.LeftProbability(x));
                }


            }

        }

        [TestMethod]
        public void DiscreteContinuousAgreement () {

            DiscreteDistribution dd = new BinomialDistribution(0.6, 7);
            Distribution cd = new DiscreteAsContinuousDistribution(dd);

            Assert.IsTrue(cd.Mean == dd.Mean);
            Assert.IsTrue(cd.StandardDeviation == dd.StandardDeviation);
            Assert.IsTrue(cd.Variance == dd.Variance);
            Assert.IsTrue(cd.Skewness == dd.Skewness);
            Assert.IsTrue(cd.Moment(5) == dd.Moment(5));
            Assert.IsTrue(cd.MomentAboutMean(5) == dd.MomentAboutMean(5));

            // this should cause an interval conversion
            Assert.IsTrue(cd.Support == dd.Support);

            Assert.IsTrue(cd.LeftProbability(4.5) == dd.LeftProbability(4));
            Assert.IsTrue(cd.RightProbability(4.5) == dd.RightProbability(4));

        }

        [TestMethod]
        public void OutsideDiscreteDistributionSupport () {
            foreach (DiscreteDistribution distribution in distributions) {
                DiscreteInterval support = distribution.Support;
                if (support.LeftEndpoint > Int32.MinValue) {
                    Assert.IsTrue(distribution.ProbabilityMass(support.LeftEndpoint - 1) == 0.0);
                    Assert.IsTrue(distribution.LeftProbability(support.LeftEndpoint - 1) == 0.0);
                    Assert.IsTrue(distribution.RightProbability(support.LeftEndpoint - 1) == 1.0);
                }
                if (support.RightEndpoint < Int32.MaxValue) {
                    Assert.IsTrue(distribution.ProbabilityMass(support.RightEndpoint + 1) == 0.0);
                    Assert.IsTrue(distribution.LeftProbability(support.RightEndpoint + 1) == 1.0);
                    Assert.IsTrue(distribution.RightProbability(support.RightEndpoint + 1) == 0.0);
                }
            }
        }

        [TestMethod]
        public void PoissonBug () {

            PoissonDistribution pd = new PoissonDistribution(0.5);
            double x = pd.InverseLeftProbability(0.7716);
            Console.WriteLine(x);

        }

        [TestMethod]
        public void DiscreteRandomValueDistribution () {

            PoissonDistribution distribution = new PoissonDistribution(3.0);
            Random rng = new Random(1);

            //foreach (DiscreteDistribution distribution in distributions) {

                Console.WriteLine(distribution.GetType().FullName);

                Histogram h = new Histogram(10);
            
                for (int i = 0; i < 100; i++) {
                    h.Add(distribution.GetRandomValue(rng));
                }
              
                for (int i = 0; i < 10; i++) {
                    Console.WriteLine(h[i].Counts);
                }

                TestResult test = h.ChiSquaredTest(distribution);
                Console.WriteLine(test.Statistic);

            


            //}


        }

    }
}
