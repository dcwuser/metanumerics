
using Meta.Numerics;
using Meta.Numerics.Statistics;

using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Test {
    
    
    /// <summary>
    ///This is a test class for BinaryContingencyTableTest and is intended
    ///to contain all BinaryContingencyTableTest Unit Tests
    ///</summary>
    [TestClass()]
    public class BinaryContingencyTableTest {


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

        private static BinaryContingencyTable CreateExperiment (double p, double q0, double q1, int N) {

            BinaryContingencyTable e = new BinaryContingencyTable();

            Random rng = new Random(1);
            for (int i = 0; i < N; i++) {

                int r, c;
                if (rng.NextDouble() < p) {
                    r = 0;
                    if (rng.NextDouble() < q0) {
                        c = 0;
                    } else {
                        c = 1;
                    }
                } else {
                    r = 1;
                    if (rng.NextDouble() < q1) {
                        c = 0;
                    } else {
                        c = 1;
                    }
                }
                e[r, c] += 1;

            }


            return (e);

        }


        [TestMethod]
        public void BinaryContingencyNullTest () {

            BinaryContingencyTable e0 = CreateExperiment(0.25, 0.33, 0.33, 200);

            Assert.IsTrue(e0.Total == 200);
            Assert.IsTrue(e0.RowTotal(0) + e0.RowTotal(1) == e0.Total);
            Assert.IsTrue(e0.ColumnTotal(0) + e0.ColumnTotal(1) == e0.Total);

            UncertainValue lnr = e0.LogOddsRatio;
            Assert.IsTrue(lnr.ConfidenceInterval(0.95).ClosedContains(0.0));

            UncertainValue r = e0.OddsRatio;
            Assert.IsTrue(r.ConfidenceInterval(0.95).ClosedContains(1.0));

            TestResult f = e0.FisherExactTest();
            Assert.IsTrue(f.RightProbability < 0.95, f.RightProbability.ToString());

        }

        [TestMethod]
        [ExpectedException(typeof(DimensionMismatchException))]
        public void BinaryContingencyInvalidConstructionTest () {

            int[,] M = new int[2,3];
            BinaryContingencyTable t = new BinaryContingencyTable(M);
        }

        [TestMethod]
        public void BinaryContingencyTest () {

            BinaryContingencyTable e1 = CreateExperiment(0.50, 0.50, 0.75, 200);

            Assert.IsTrue(e1.RowTotal(0) + e1.RowTotal(1) == e1.Total);
            Assert.IsTrue(e1.ColumnTotal(0) + e1.ColumnTotal(1) == e1.Total);

            UncertainValue lnr = e1.LogOddsRatio;
            Assert.IsFalse(lnr.ConfidenceInterval(0.95).ClosedContains(0.0));

            UncertainValue r = e1.OddsRatio;
            Assert.IsFalse(r.ConfidenceInterval(0.95).ClosedContains(1.0));

            TestResult p = e1.PearsonChiSquaredTest();
            Assert.IsTrue(p.LeftProbability > 0.95, p.RightProbability.ToString());

            TestResult f = e1.FisherExactTest();
            Assert.IsTrue(f.RightProbability > 0.95);

        }

    }
}
