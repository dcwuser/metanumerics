using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Statistics;


namespace Test {
    
    
    /// <summary>
    ///This is a test class for UncertainMathTest and is intended
    ///to contain all UncertainMathTest Unit Tests
    ///</summary>
    [TestClass()]
    public class UncertainMathTest {


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


        UncertainValue a = new UncertainValue(2.0, 0.5);

        [TestMethod]
        public void UncertainMathFunctionValuesTest () {
            // We need to move to IsNearlyEqual because exact equality fails in release mode, even though all-digit printed values are the same
            // I suspect this has to do with the optomized comparison including bits beyond those that strictly belong to the value and are random
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Cos(a).Value, Math.Cos(a.Value)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Sin(a).Value, Math.Sin(a.Value)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Tan(a).Value, Math.Tan(a.Value)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Cosh(a).Value, Math.Cosh(a.Value)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Sinh(a).Value, Math.Sinh(a.Value)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Tanh(a).Value, Math.Tanh(a.Value)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Log(a).Value, Math.Log(a.Value)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Exp(a).Value, Math.Exp(a.Value)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Sqrt(a).Value, Math.Sqrt(a.Value)));
        }

        [TestMethod]
        public void UncertainMathSqrtTest () {
            UncertainValue sa = UncertainMath.Sqrt(a);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(sa.RelativeUncertainty, a.RelativeUncertainty / 2.0), String.Format("a={0} sa={1}", a, sa));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Pow(sa, 2.0).Value, a.Value));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Pow(sa, 2.0).Uncertainty, a.Uncertainty));
        }

        [TestMethod]
        public void UncertainMathLogTest () {
            UncertainValue la = UncertainMath.Log(a);
            Assert.IsTrue(la.Uncertainty == a.RelativeUncertainty);
            Assert.IsTrue(UncertainMath.Exp(la) == a);
        }

    }
}
