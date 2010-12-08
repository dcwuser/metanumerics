using System;
using System.Collections;

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
        public void FactorialSpecialCases () {
            Assert.IsTrue(AdvancedIntegerMath.Factorial(0) == 1);
            Assert.IsTrue(AdvancedIntegerMath.Factorial(1) == 1);
            Assert.IsTrue(AdvancedIntegerMath.Factorial(2) == 2);
            Assert.IsTrue(AdvancedIntegerMath.Factorial(3) == 6);
            Assert.IsTrue(AdvancedIntegerMath.Factorial(4) == 24);
        }

        [TestMethod()]
        public void FactorialRecurrence () {

            foreach (int n in TestUtilities.GenerateIntegerValues(1, 30, 5)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedIntegerMath.Factorial(n), n * AdvancedIntegerMath.Factorial(n - 1)));
            }
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void FactorialNegativeArgumentTest () {
            AdvancedIntegerMath.Factorial(-1);
        }

        [TestMethod]
        public void DoubleFactorialSpecialCases () {
            Assert.IsTrue(AdvancedIntegerMath.DoubleFactorial(0) == 1);
            Assert.IsTrue(AdvancedIntegerMath.DoubleFactorial(1) == 1);
            Assert.IsTrue(AdvancedIntegerMath.DoubleFactorial(2) == 2);
            Assert.IsTrue(AdvancedIntegerMath.DoubleFactorial(3) == 3);
            Assert.IsTrue(AdvancedIntegerMath.DoubleFactorial(4) == 8);
            Assert.IsTrue(AdvancedIntegerMath.DoubleFactorial(5) == 15);
        }

        [TestMethod]
        public void FactorialDoubleFactorialRelationship () {

            // n! = n!! (n-1)!!
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedIntegerMath.LogFactorial(n),
                    AdvancedIntegerMath.LogDoubleFactorial(n) + AdvancedIntegerMath.LogDoubleFactorial(n - 1)
                ));
            }
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
        public void BinomialCoefficientSpecialCases () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 30, 5)) {
                Assert.IsTrue(AdvancedIntegerMath.BinomialCoefficient(n, 0) == 1);
                Assert.IsTrue(AdvancedIntegerMath.BinomialCoefficient(n, 1) == n);
                Assert.IsTrue(AdvancedIntegerMath.BinomialCoefficient(n, n - 1) == n);
                Assert.IsTrue(AdvancedIntegerMath.BinomialCoefficient(n, n) == 1);
            }
        }

        [TestMethod]
        public void BinomialCoefficientSum () {

            foreach (int n in TestUtilities.GenerateIntegerValues(1, 30, 5)) {
                long sum = 0;
                long product = 1;
                for (int m = 0; m <= n; m++) {
                    sum += AdvancedIntegerMath.BinomialCoefficient(n, m);
                    if (m > 0) product *= 2;
                }
                Console.WriteLine("n={0}, sum={0}, product={1}", n, sum, product);
                Assert.IsTrue(sum == product);
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
            Assert.IsTrue(AdvancedIntegerMath.GCF(a, b) * AdvancedIntegerMath.LCM(a, b) == a * b);
        }

        [TestMethod]
        public void IntegerPartitionCounts () {

            // these counts are from Table 21.5 of Abromowitz & Stegun
            Assert.IsTrue(CountValues(AdvancedIntegerMath.Partitions(1)) == 1);
            Assert.IsTrue(CountValues(AdvancedIntegerMath.Partitions(2)) == 2);
            Assert.IsTrue(CountValues(AdvancedIntegerMath.Partitions(3)) == 3);
            Assert.IsTrue(CountValues(AdvancedIntegerMath.Partitions(4)) == 5);
            Assert.IsTrue(CountValues(AdvancedIntegerMath.Partitions(5)) == 7);
            Assert.IsTrue(CountValues(AdvancedIntegerMath.Partitions(6)) == 11);
            Assert.IsTrue(CountValues(AdvancedIntegerMath.Partitions(7)) == 15);
            Assert.IsTrue(CountValues(AdvancedIntegerMath.Partitions(8)) == 22);

        }

        private int CountValues (IEnumerable enumeration) {
            int count = 0;
            foreach (object value in enumeration) count++;
            return(count);
        }

        [TestMethod]
        public void IntegerPartitionSums () {

            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {

                foreach (int[] partition in AdvancedIntegerMath.Partitions(n)) {
                    int s = 0;
                    foreach (int i in partition) s += i;
                    Assert.IsTrue(s == n);
                }

            }

        }

        [TestMethod]
        public void PowModTest () {
            Console.WriteLine(AdvancedIntegerMath.PowMod(130000000, 670, 59));
            Console.WriteLine(AdvancedIntegerMath.PowMod(2, 360360, 779167));
            // some examples from http://modular.math.washington.edu/edu/2007/spring/ent/ent-html/node81.html
            Assert.IsTrue(AdvancedIntegerMath.PowMod(2, 360360, 779167) == 584876);
        }

    }
}
