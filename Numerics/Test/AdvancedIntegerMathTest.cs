using System;
using Meta.Numerics.Functions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Test {
    
    
    /// <summary>
    ///This is a test class for AdvancedIntegerMathTest and is intended
    ///to contain all AdvancedIntegerMathTest Unit Tests
    ///</summary>
    [TestClass()]
    public class AdvancedIntegerMathTest {


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

        // Factorial

        private int nmax = 30;

        [TestMethod]
        public void FactorialSpecialCasesTest () {
            Assert.IsTrue(AdvancedIntegerMath.Factorial(0) == 1);
            Assert.IsTrue(AdvancedIntegerMath.Factorial(1) == 1);
            Assert.IsTrue(AdvancedIntegerMath.Factorial(2) == 2);
            Assert.IsTrue(AdvancedIntegerMath.Factorial(3) == 6);
        }

        [TestMethod()]
        public void FactorialRecurrenceTest () {
            for (int n = 1; n < nmax; n++) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedIntegerMath.Factorial(n), n * AdvancedIntegerMath.Factorial(n - 1)));
            }
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void FactorialNegativeArgumentTest () {
            AdvancedIntegerMath.Factorial(-1);
        }

        // Binomial coefficients

        [TestMethod]
        public void BinomialCoefficientRecurrenceTest () {
            // this is the Pascal triangle recurrence
            for (int n = 1; n < nmax; n++) {
                for (int m = 1; m < n; m++) {
                    //Assert.AreEqual<long>(AdvancedIntegerMath.BinomialCoefficient(n + 1, m), AdvancedIntegerMath.BinomialCoefficient(n, m) + AdvancedIntegerMath.BinomialCoefficient(n, m - 1));
                    Assert.AreEqual<long>(AdvancedIntegerMath.BinomialCoefficient(n + 1, m + 1), AdvancedIntegerMath.BinomialCoefficient(n, m) + AdvancedIntegerMath.BinomialCoefficient(n, m+1));

                }
            }
        }

        [TestMethod]
        public void BinomialCoefficientInequalityTest () {
            for (int n = 1; n < nmax; n++) {
                for (int m = 1; m < n; m++) {
                    double lower = Math.Pow(1.0 * n / m, m);
                    double upper = Math.Pow(Math.E * n / m, m);
                    double value = AdvancedIntegerMath.BinomialCoefficient(n, m);
                    Assert.IsTrue((lower <= value) && (value <= upper));
                }
            }
        }

        [TestMethod]
        public void BinomialCoefficientSpecialCasesTest () {
            for (int n = 0; n < nmax; n++) {
                Assert.IsTrue(AdvancedIntegerMath.BinomialCoefficient(n, 0) == 1);
                Assert.IsTrue(AdvancedIntegerMath.BinomialCoefficient(n, n) == 1);
            }
        }

        [TestMethod]
        public void BinomialCoefficientSumTest () {
            for (int n = 0; n < nmax; n++) {
                long sum = 0;
                long product = 1;
                for (int m = 0; m <= n; m++) {
                    sum += AdvancedIntegerMath.BinomialCoefficient(n, m);
                    if (m > 0) product *= 2;
                }
                Console.WriteLine("n={0}, sum={0}, product={1}", n, sum, product);
                Assert.IsTrue(sum == product, String.Format("n={0}, sum={0}, product={1}", n, sum, product));
            }
        }

        [TestMethod]
        public void BinomialCoefficientSumOfSquaresTest () {
            for (int n = 0; n < nmax; n++) {
                Console.WriteLine(n);
                long sum = 0;
                for (int m = 0; m <= n; m++) {
                    Console.WriteLine("  {0} {1}", m, AdvancedIntegerMath.BinomialCoefficient(n, m));
                    long value = AdvancedIntegerMath.BinomialCoefficient(n, m);
                    sum += value * value;
                }
                Console.WriteLine(AdvancedIntegerMath.BinomialCoefficient(2 * n, n));
                Assert.IsTrue(sum == AdvancedIntegerMath.BinomialCoefficient(2 * n, n));
            }
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void BinomialCoefficientInvalidArgumentTest1 () {
            AdvancedIntegerMath.BinomialCoefficient(-1, -1);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void BinomialCoefficientInvalidArgumentTest2 () {
            AdvancedIntegerMath.BinomialCoefficient(2, -1);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void BinomialCoefficientInvalidArgumentTest3 () {
            AdvancedIntegerMath.BinomialCoefficient(2, 3);
        }

        // GCD and LCM

        [TestMethod]
        public void GcdLcmTest () {
            long a = 3457832408;
            long b = 56789309233;
            Assert.AreEqual<long>(AdvancedIntegerMath.GCF(a, b) * AdvancedIntegerMath.LCM(a, b), a * b);
        }

        [TestMethod]
        public void IntegerParticianTest () {

            IntegerPartitionEnumerator e = new IntegerPartitionEnumerator(120);
            int P = 0;
            while (e.MoveNext()) {
                P++;
                /*
                int[] p = e.Current;
                for (int i = 0; i < p.Length; i++) {
                    Console.Write("{0} ", p[i]);
                }
                Console.WriteLine();
                */
            }
            Console.WriteLine(P);

        }

        [TestMethod]
        public void PowModTest () {
            Console.WriteLine(AdvancedIntegerMath.PowMod(130000000, 670, 59));
        }

    }
}
