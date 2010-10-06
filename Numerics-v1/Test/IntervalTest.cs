using System;
using Meta.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Test {
    
    
    /// <summary>
    ///This is a test class for IntervalTest and is intended
    ///to contain all IntervalTest Unit Tests
    ///</summary>
    [TestClass()]
    public class IntervalTest {


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


        private double a = -3.1;
        private double b = -2.0;
        private double c = 2.7;
        private double d = 5.1;

        [TestMethod]
        public void IntervalZeroWidthTest () {
            Interval aa = Interval.FromEndpoints(a, a);
            Assert.AreEqual<double>(aa.Width, 0.0);
            Assert.AreEqual<double>(aa.LeftEndpoint, aa.RightEndpoint);
            Assert.AreEqual<double>(aa.Midpoint, a);
        }

        [TestMethod()]
        public void IntervalWidthTest () {
            Interval ab = Interval.FromEndpoints(a, b);
            Assert.AreEqual<double>(ab.Width, Math.Abs(a - b));
        }

        /// <summary>
        ///A test for To
        ///</summary>
        [TestMethod()]
        public void ToTest () {
            Interval ab = Interval.FromEndpoints(a, b);
            Assert.AreEqual<double>(ab.RightEndpoint, b);
        }

        /// <summary>
        ///A test for Midpoint
        ///</summary>
        [TestMethod()]
        public void MidpointTest () {
            Interval ab = Interval.FromEndpoints(a, b);
            Assert.AreEqual(ab.Midpoint, (a + b) / 2.0);
        }

        /// <summary>
        ///A test for From
        ///</summary>
        [TestMethod()]
        public void FromTest () {
            Interval ab = Interval.FromEndpoints(a, b);
            Assert.AreEqual<double>(ab.LeftEndpoint, a);
        }

        /// <summary>
        ///A test for OpenContains
        ///</summary>
        [TestMethod()]
        public void IntervalOpenContainsTest () {
            Interval ac = Interval.FromEndpoints(a, c);
            Assert.IsFalse(ac.OpenContains(a));
            Assert.IsFalse(ac.OpenContains(c));
            Assert.IsTrue(ac.OpenContains(b));
            Assert.IsFalse(ac.OpenContains(d));
        }

        /// <summary>
        ///A test for FromMidpointAndWidth
        ///</summary>
        [TestMethod()]
        public void FromMidpointAndWidthTest () {
        }

        /// <summary>
        ///A test for FromEndpoints
        ///</summary>
        [TestMethod()]
        public void IntervalFromEndpointsTest () {
            Interval ab = Interval.FromEndpoints(a, b);
            Interval ba = Interval.FromEndpoints(b, a);
            Assert.IsTrue(ab.LeftEndpoint == ba.LeftEndpoint);
            Assert.IsTrue(ab.RightEndpoint == ba.RightEndpoint);
            Assert.IsTrue(ab.Width == ba.Width);
        }

        /// <summary>
        ///A test for FromEndpointAndWidth
        ///</summary>
        [TestMethod()]
        public void FromEndpointAndWidthTest () {
        }

        /// <summary>
        ///A test for ClosedContains
        ///</summary>
        [TestMethod()]
        public void IntervalClosedContainsTest () {
            Interval ac = Interval.FromEndpoints(a, c);
            Assert.IsTrue(ac.ClosedContains(a));
            Assert.IsTrue(ac.ClosedContains(c));
            Assert.IsTrue(ac.ClosedContains(b));
            Assert.IsFalse(ac.ClosedContains(d));
        }

    }
}
