using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Statistics;

namespace Test {

    internal class TestDerivative {

        public TestDerivative (Function<double, double> function, Function<double, double> derivative) {
            this.Function = function;
            this.Derivative = derivative;
        }

        public Function<double, double> Function { get; internal set; }

        public Function<double, double> Derivative { get; internal set; }

    }

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


        private TestDerivative[] derivatives = new TestDerivative[] {

            new TestDerivative( delegate (double x) { return(x*x); }, delegate (double x) { return(2.0 * x); } ),

            new TestDerivative( delegate (double x) { return(x); }, delegate (double x) { return(1.0); } ),

            new TestDerivative( delegate (double x) { return(1.0); }, delegate (double x) { return(0.0); } ),

            // (0,
            //new TestDerivative(delegate (double x) { return(1.0 / x); }, delegate (double x) {return(-1.0 / (x*x)); } ),

            // [0,
            //new TestDerivative( delegate (double x) { return(Math.Sqrt(x)); } , delegate (double x) { return(1.0 / 2.0 / Math.Sqrt(x)); } ),

            new TestDerivative( delegate (double x) { return(Math.Sin(x)); }, delegate (double x) { return(Math.Cos(x)); } ),

            //new TestDerivative( delegate (double x) { return(Math.Cos(x)); }, delegate (double x) { return(-Math.Sin(x)); } ),

            //new TestDerivative( delegate (double x) { return(Math.Exp(x)); }, delegate (double x) { return(Math.Exp(x)); } ), 

            // (0, 
            new TestDerivative( delegate (double x) { return(Math.Log(x)); }, delegate (double x) { return(1.0 / x); } ),

            // (0,
            //new TestDerivative (delegate (double x) { return(Math.Pow(x, x)); }, delegate (double x) { return(Math.Pow(x, x) * ( Math.Log(x) + 1.0)); } ), 

            new TestDerivative (delegate (double x) { return(Math.Atan(x)); }, delegate (double x) { return(1.0 / ( 1.0 + x * x)); } )

        };

        [TestMethod()]
        public void DifferentiationTest () {

            int i = 0;
            foreach (TestDerivative test in derivatives) {
                i++;

                foreach (double x in TestUtilities.GenerateRealValues(0.1, 100.0, 5)) {

                    Console.WriteLine("{0} {1}", i, x);

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        FunctionMath.Differentiate(test.Function, x),
                        test.Derivative(x),
                        Math.Pow(2, -42)
                    ));

                }

            }


        }
    }
}
