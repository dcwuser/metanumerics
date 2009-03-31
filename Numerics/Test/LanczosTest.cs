using Meta.Numerics.Functions;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace Test
{


    [TestClass()]
    public class LanczosTest {


        private TestContext testContextInstance;

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


        /*
        [TestMethod()]
        public void LanczosCoefficientsTest () {
            int n = 14;
            double g = 5.2421875;
            decimal[] fs = Lanczos.LanczosCoefficients(n, g);
            for (int i=0; i<fs.Length; i++) {
                Console.WriteLine("{0} {1:g20}", i, fs[i]);
            }
        }
        */
    }
}
