using System;
using Meta.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Test
{
    
    
    /// <summary>
    ///This is a test class for ComplexTest and is intended
    ///to contain all ComplexTest Unit Tests
    ///</summary>
    [TestClass()]
    public class ComplexTest {


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


        private Complex a = new Complex(1.5, -2.2);
        private Complex b = new Complex(-7.4, 3.8);
        private Complex c = new Complex(5.9, 6.0);

        /// <summary>
        ///A test for Re
        ///</summary>
        [TestMethod()]
        public void ReTest () {
            Complex z = new Complex(1.5, -2.2);
            Assert.AreEqual<double>(z.Re, 1.5);
        }

        /// <summary>
        ///A test for Im
        ///</summary>
        [TestMethod()]
        public void ImTest () {
            Complex z = new Complex(1.5, -2.2);
            Assert.AreEqual<double>(z.Im, -2.2);
        }

        /// <summary>
        ///A test for Conjugate
        ///</summary>
        [TestMethod()]
        public void ComplexConjugateTest () {

            Assert.AreEqual<double>(a.Re, a.Conjugate.Re);
            Assert.AreEqual<double>(a.Im, -a.Conjugate.Im);
            Assert.IsTrue((a * a.Conjugate).Im == 0.0);
        }

        /// <summary>
        ///A test for op_UnaryNegation
        ///</summary>
        [TestMethod()]
        public void ComplexNegationTest () {
            Assert.IsTrue((a + (-a)) == 0);
        }

        /// <summary>
        ///A test for op_Subtraction
        ///</summary>
        [TestMethod()]
        public void ComplexSubtractionTest () {
            Assert.IsTrue((a - a) == 0);
            Assert.IsTrue((0 - a) == (-a));
            Assert.IsTrue((a - b) == -(b - a));
        }


        [TestMethod()]
        public void ComplexMultiplyTest () {
            Assert.IsTrue(a * b == b * a);
            Assert.IsTrue(a * (b + c) == (a * b + a * c));
        }

        [TestMethod]
        public void ComplexExplicitCastTest () {
            Complex a = new Complex(1.0, 0.0);
            double b = (double) a;
            Assert.IsTrue(b == a.Re);
        }

        [TestMethod]
        [ExpectedException(typeof(InvalidCastException))]
        public void ComplexExplicitCastFailureTest () {
            Complex a = new Complex(1.0,1.0);
            double b = (double) a;
        }

        [TestMethod]
        public void ComplexEqualityTest () {
            Assert.IsTrue(a == a);
            Assert.IsFalse(a == b);
        }

        public void ComplexInequalityTest () {
            Assert.IsFalse(a != a);
            Assert.IsTrue(a != b);
        }


        /// <summary>
        ///A test for op_Addition
        ///</summary>
        [TestMethod()]
        public void ComplexAdditionTest () {
            Assert.IsTrue(a + b == b + a);
            Assert.IsTrue(a + (-b) == a - b);
        }


        [TestMethod()]
        public void EqualsTest () {
            Assert.IsTrue(Complex.Equals(a, a));
            Assert.IsFalse(Complex.Equals(a, b));
        }

    }
}
