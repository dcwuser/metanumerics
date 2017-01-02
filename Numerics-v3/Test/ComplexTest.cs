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

        [TestMethod]
        public void ComplexComponents () {
            Complex z = new Complex(-1.0, 2.0);
            Assert.IsTrue(z.Re == -1.0);
            Assert.IsTrue(z.Im == 2.0);
        }

        [TestMethod]
        public void ComplexEquality () {

            Complex ac = a;

            // Equality operator
            Assert.IsTrue(a == ac);
            Assert.IsTrue(ac == a);
            Assert.IsFalse(a == b);
            Assert.IsFalse(b == a);

            // Inequality operator
            Assert.IsFalse(a != ac);
            Assert.IsFalse(ac != a);
            Assert.IsTrue(a != b);
            Assert.IsTrue(b != a);

            // Equals method
            Assert.IsTrue(a.Equals(ac));
            Assert.IsFalse(a.Equals(b));
            Assert.IsFalse(a.Equals(new object()));
            Assert.IsFalse(a.Equals(null));

            // Hash
            Assert.IsTrue(a.GetHashCode() == ac.GetHashCode());

        }

        [TestMethod()]
        public void ComplexConjugation () {
            Assert.IsTrue(a.Conjugate.Re == a.Re);
            Assert.IsTrue(a.Conjugate.Im == -a.Im);
            Assert.IsTrue((a * a.Conjugate).Im == 0.0);
        }

        [TestMethod]
        public void ComplexNegation () {
            Assert.IsTrue((a + (-a)) == 0);
        }

        [TestMethod]
        public void ComplexSubtraction () {
            Assert.IsTrue((a - a) == 0);
            Assert.IsTrue((0 - a) == (-a));
            Assert.IsTrue((a - b) == -(b - a));
        }


        [TestMethod]
        public void ComplexMultiplication () {
            Assert.IsTrue(a * b == b * a);
            Assert.IsTrue(a * (b + c) == (a * b + a * c));
        }

        [TestMethod]
        public void ComplexExplicitCast () {
            Complex a = new Complex(1.0, 0.0);
            double b = (double) a;
            Assert.IsTrue(b == a.Re);
        }

        [TestMethod]
        [ExpectedException(typeof(InvalidCastException))]
        public void ComplexExplicitCastFailure () {
            Complex a = new Complex(1.0,1.0);
            double b = (double) a;
        }

        [TestMethod]
        public void ComplexAddition () {
            Assert.IsTrue(a + b == b + a);
            Assert.IsTrue(a + (-b) == a - b);
        }


        [TestMethod()]
        public void ComplexEquals () {
            Assert.IsTrue(Complex.Equals(a, a));
            Assert.IsFalse(Complex.Equals(a, b));
        }

    }
}
