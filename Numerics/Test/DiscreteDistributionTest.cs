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
            new PoissonDistribution(4.5), new PoissonDistribution(400.0),
            new DiscreteUniformDistribution(5, 11),
            new GeometricDistribution(0.6),
            new NegativeBinomialDistribution(7.8, 0.4)
        };

        public static DiscreteDistribution[] GetDistributions () {
            return (new DiscreteDistribution[] {
                new BernoulliDistribution(0.1),
                new BinomialDistribution(0.2, 30), new BinomialDistribution(0.4, 5),
                new PoissonDistribution(4.5), new PoissonDistribution(400.0),
                new DiscreteUniformDistribution(5, 11),
                new GeometricDistribution(0.6),
                new NegativeBinomialDistribution(7.8, 0.4)
            });
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

                    Console.WriteLine("{0} {1}", distribution.GetType().Name, k);

                    double DP = distribution.ProbabilityMass(k);
                    Assert.IsTrue(DP >= 0.0); Assert.IsTrue(DP <= 1.0);

                    double P = distribution.LeftInclusiveProbability(k);
                    double Q = distribution.RightExclusiveProbability(k);
                    Console.WriteLine("{0} {1} {2}", P, Q, P + Q);

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
                    int k = distribution.InverseLeftProbability(P);
                    Console.WriteLine("{0} P={1} K={2} P(k<K)={3} P(k<=K)={4}", distribution.GetType().Name, P, k, distribution.LeftExclusiveProbability(k), distribution.LeftInclusiveProbability(k));
                    Console.WriteLine("    {0} {1}", distribution.LeftExclusiveProbability(k) < P, P <= distribution.LeftInclusiveProbability(k));
                    Assert.IsTrue(distribution.LeftExclusiveProbability(k) < P);
                    Assert.IsTrue(P <= distribution.LeftInclusiveProbability(k));
                }


            }

        }

        /*
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
            Assert.IsTrue(TestUtilities.IsNearlyEqual(cd.Support.LeftEndpoint, dd.Minimum));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(cd.Support.RightEndpoint, dd.Maximum));

            //Assert.IsTrue(cd.LeftProbability(4.5) == dd.LeftProbability(4));
            //Assert.IsTrue(cd.RightProbability(4.5) == dd.RightProbability(4));

            // Switch LeftProbablity for discrete distributions to be exclusive.
            // This is already the case for internal distributions used for exact null distributions, but not for public discrete distributions.

        }
        */

        [TestMethod]
        public void OutsideDiscreteDistributionSupport () {
            foreach (DiscreteDistribution distribution in distributions) {
                if (distribution.Minimum > Int32.MinValue) {
                    Assert.IsTrue(distribution.ProbabilityMass(distribution.Minimum - 1) == 0.0);
                    Assert.IsTrue(distribution.LeftInclusiveProbability(distribution.Minimum - 1) == 0.0);
                    Assert.IsTrue(distribution.RightExclusiveProbability(distribution.Minimum - 1) == 1.0);
                }
                if (distribution.Maximum < Int32.MaxValue) {
                    Assert.IsTrue(distribution.ProbabilityMass(distribution.Maximum + 1) == 0.0);
                    Assert.IsTrue(distribution.LeftInclusiveProbability(distribution.Maximum + 1) == 1.0);
                    Assert.IsTrue(distribution.RightExclusiveProbability(distribution.Maximum + 1) == 0.0);
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
        public void DiscreteDistributionRandomValues () {

            foreach (DiscreteDistribution distribution in distributions) {

                int max;
                if (distribution.Maximum < 128) {
                    max = distribution.Maximum + 1;
                } else {
                    max = (int) Math.Round(distribution.Mean + 2.0 * distribution.StandardDeviation);
                }
                Console.WriteLine("{0} {1}", distribution.GetType().Name, max);

                Histogram h = new Histogram(max);

                Random rng = new Random(314159265);
                for (int i = 0; i < 1024; i++) h.Add(distribution.GetRandomValue(rng));
                TestResult result = h.ChiSquaredTest(distribution);
                Console.WriteLine("{0} {1}", result.Statistic, result.RightProbability);
                Assert.IsTrue(result.RightProbability > 0.05);

            }

        }

        [TestMethod]
        [Ignore]
        public void DiscreteDistributionBase () {

            DiscreteDistribution D = new DiscreteTestDistribution();

            double P = 0.0;
            double M1 = 0.0;
            double M2 = 0.0;
            for (int k = 1; k <= 3; k++) {
                P += D.ProbabilityMass(k);
                M1 += k * D.ProbabilityMass(k);
                M2 += k * k * D.ProbabilityMass(k);
            }
            Assert.IsTrue(TestUtilities.IsNearlyEqual(P, 1.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(D.Moment(0), 1.0));

            double C2 = M2 - M1 * M1;
            Assert.IsTrue(TestUtilities.IsNearlyEqual(M1, D.Mean));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(M1, D.Moment(1)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(M2, D.Moment(2)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(C2, D.Variance));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(C2, D.MomentAboutMean(2)));

            Assert.IsTrue(D.InverseLeftProbability(D.LeftInclusiveProbability(2)) == 2);

        }

    }

    // a minimal implementation to test base methods on abstract DiscreteDistribution class

    public class DiscreteTestDistribution : DiscreteDistribution {

#if PAST
        public override DiscreteInterval Support {
            get { return DiscreteInterval.FromEndpoints(1, 3); }
        }
#endif

        public override int Minimum {
            get { return (1); }
        }

        public override int Maximum {
            get { return (3); }
        }

        public override double ProbabilityMass (int k) {
            switch (k) {
                case 1:
                    return (1.0 / 6.0);
                case 2:
                    return (2.0 / 6.0);
                case 3:
                    return (3.0 / 6.0);
                default:
                    return (0.0);
            }
        }

    }

}
