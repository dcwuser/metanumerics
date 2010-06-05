using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;


namespace Test {
    
    
    /// <summary>
    ///This is a test class for FunctionMathTest and is intended
    ///to contain all FunctionMathTest Unit Tests
    ///</summary>
    [TestClass()]
    public class DifferentiateTest {

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


        [TestMethod()]
        public void DifferentiationTest () {

            Function<double, double> f = delegate(double x) {
                return (AdvancedMath.Gamma(x));
                //return (Math.Sin(x));
            };
            double dfdx = FunctionMath.Differentiate(f, 1.0);
            Console.WriteLine(dfdx);

        }
    }
}
