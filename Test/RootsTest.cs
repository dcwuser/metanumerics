using Meta.Numerics.Functions;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Meta.Numerics;

namespace Test
{
    
    
    /// <summary>
    ///This is a test class for RootsTest and is intended
    ///to contain all RootsTest Unit Tests
    ///</summary>
    [TestClass()]
    public class RootsTest {


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

        [TestMethod]
        public void FindRootOfEi () {
            Function<double, double> f = new Function<double, double>(AdvancedMath.IntegralEi);
            double x = FunctionMath.FindZero(f, 0.5); // fails for x=1, tries to cross y-axis boundary
            Assert.IsTrue(TestUtilities.IsNearlyEqual(x, 0.37250741078136663446));
        }

        [TestMethod]
        public void FindRootOfPsi () {
            Function<double, double> f = new Function<double, double>(AdvancedMath.Psi);
            double x = FunctionMath.FindZero(f, 1.5);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(x, 1.46163214496836234126));
        }

        [TestMethod]
        public void FindRootOfJ0 () {
            Function<double, double> f = delegate(double x) {
                return (AdvancedMath.BesselJ(0, x));
            };
            double y = FunctionMath.FindZero(f, Interval.FromEndpoints(2.0, 4.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(y, 2.40482555769577276862)); 
        }

    }
}
