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
        public void IntervalZeroWidth () {
            Interval aa = Interval.FromEndpoints(a, a);
            Assert.IsTrue(aa.Width == 0.0);
            Assert.IsTrue(aa.LeftEndpoint == aa.RightEndpoint);
            Assert.IsTrue(aa.Midpoint == a);
        }

        [TestMethod()]
        public void IntervalWidth () {
            Interval ab = Interval.FromEndpoints(a, b);
            Assert.IsTrue(ab.Width == Math.Abs(a - b));
        }

        [TestMethod]
        public void IntervalEndpoints () {
            Interval ab = Interval.FromEndpoints(a, b);
            Assert.IsTrue(ab.LeftEndpoint == a);
            Assert.IsTrue(ab.RightEndpoint == b);
        }

        [TestMethod()]
        public void ToTest () {
            Interval ab = Interval.FromEndpoints(a, b);
            Assert.AreEqual<double>(ab.RightEndpoint, b);
        }

        [TestMethod]
        public void IntervalMidpoint () {
            Interval ab = Interval.FromEndpoints(a, b);
            Assert.IsTrue(ab.Midpoint == (a + b) / 2.0);
        }

        /// <summary>
        ///A test for From
        ///</summary>
        [TestMethod()]
        public void FromTest () {
            Interval ab = Interval.FromEndpoints(a, b);
            Assert.AreEqual<double>(ab.LeftEndpoint, a);
        }

        [TestMethod]
        public void IntervalOpenContains () {
            Interval ac = Interval.FromEndpoints(a, c);
            Assert.IsFalse(ac.OpenContains(a));
            Assert.IsFalse(ac.OpenContains(c));
            Assert.IsTrue(ac.OpenContains(b));
            Assert.IsFalse(ac.OpenContains(d));
        }

        [TestMethod]
        public void IntervalFromMidpointAndWidth () {
            double m = 1.0;
            double w = 4.0;
            Interval ab = Interval.FromMidpointAndWidth(m, w);
            Assert.IsTrue(ab.Width == w);
            Assert.IsTrue(ab.Midpoint == m);
            Assert.IsTrue(ab.LeftEndpoint == m - w / 2.0);
            Assert.IsTrue(ab.RightEndpoint == m + w / 2.0);
        }

        [TestMethod]
        public void IntervalFromEndpoints () {
            Interval ab = Interval.FromEndpoints(a, b);
            Interval ba = Interval.FromEndpoints(b, a);
            Assert.IsTrue(ab.LeftEndpoint == ba.LeftEndpoint);
            Assert.IsTrue(ab.RightEndpoint == ba.RightEndpoint);
            Assert.IsTrue(ab.Width == ba.Width);
            Assert.IsTrue(ab == ba);
            Assert.IsFalse(ab != ba);
            Assert.IsTrue(ab.GetHashCode() == ba.GetHashCode());
        }

        [TestMethod]
        public void IntervalClosedContains () {
            Interval ac = Interval.FromEndpoints(a, c);
            Assert.IsTrue(ac.ClosedContains(a));
            Assert.IsTrue(ac.ClosedContains(c));
            Assert.IsTrue(ac.ClosedContains(b));
            Assert.IsFalse(ac.ClosedContains(d));
        }

    }
}
