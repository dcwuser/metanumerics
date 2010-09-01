using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Statistics;


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
            new BinomialDistribution(0.2, 30)
        };

        [TestMethod]
        public void DiscreteDistributionMomentSpecialCases () {

            foreach (DiscreteDistribution distribution in distributions) {
                Assert.IsTrue(distribution.Moment(0) == 1.0);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(distribution.Moment(1), distribution.Mean));
                Assert.IsTrue(distribution.MomentAboutMean(0) == 1.0);
                Assert.IsTrue(distribution.MomentAboutMean(1) == 0.0);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(distribution.MomentAboutMean(2), distribution.Variance));                
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
                    distribution.ExpectationValue(delegate (double k) { return (1.0); }), 1.0
                ));
            }
        }

        [TestMethod]
        public void DiscreteDistributionMean () {
            foreach (DiscreteDistribution distribution in distributions) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    distribution.ExpectationValue(delegate(double k) { return (k); }), distribution.Mean
                ));
            }
        }

        [TestMethod]
        public void DiscreteDistributionVariance () {
            foreach (DiscreteDistribution distribution in distributions) {
                double m = distribution.Mean;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    distribution.ExpectationValue(delegate(double x) { return (Math.Pow(x-m, 2)); }), distribution.Variance
                ));
            }
        }

    }
}
