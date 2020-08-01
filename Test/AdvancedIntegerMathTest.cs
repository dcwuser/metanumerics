using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics.Distributions;

using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Test {
    
    [TestClass]
    public class AdvancedIntegerMathTest {

        // Factorial

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

        [TestMethod]
        public void HarmonicSpecialCases () {

            Assert.IsTrue(AdvancedIntegerMath.HarmonicNumber(0) == 0.0);
            Assert.IsTrue(AdvancedIntegerMath.HarmonicNumber(1) == 1.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedIntegerMath.HarmonicNumber(2), 3.0 / 2.0));

        }

        [TestMethod]
        public void HarmonicPsiAgreement () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedIntegerMath.HarmonicNumber(n - 1),
                    AdvancedMath.Psi(n) + AdvancedMath.EulerGamma
                ));
            }
        }

        [TestMethod]
        public void PrimeSpecialCases () {

            // small numbers
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(1));
            Assert.IsTrue(AdvancedIntegerMath.IsPrime(2));
            Assert.IsTrue(AdvancedIntegerMath.IsPrime(3));
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(4));
            Assert.IsTrue(AdvancedIntegerMath.IsPrime(5));
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(6));
            Assert.IsTrue(AdvancedIntegerMath.IsPrime(7));
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(8));
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(9));
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(10));

            // Mersene prime candidates
            Assert.IsTrue(AdvancedIntegerMath.IsPrime(3)); // 2^2 - 1
            Assert.IsTrue(AdvancedIntegerMath.IsPrime(7)); // 2^3 - 1
            Assert.IsTrue(AdvancedIntegerMath.IsPrime(31)); // 2^5 - 1
            Assert.IsTrue(AdvancedIntegerMath.IsPrime(127)); // 2^7 - 1
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(2047)); // 2^11 - 1
            Assert.IsTrue(AdvancedIntegerMath.IsPrime(8191)); // 2^13 - 1
            Assert.IsTrue(AdvancedIntegerMath.IsPrime(131071)); // 2^17 - 1
            Assert.IsTrue(AdvancedIntegerMath.IsPrime(524287)); // 2^19 - 1
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(8388607)); // 2^23 - 1
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(536870911)); // 2^29 - 1
            Assert.IsTrue(AdvancedIntegerMath.IsPrime(2147483647)); // 2^31 - 1

            // Pseudoprimes
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(2047));
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(3277));
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(4033));
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(121));
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(703));
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(1891));
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(781));
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(1541));

            // Fermat primes
            Assert.IsTrue(AdvancedIntegerMath.IsPrime(5)); // 2^2 + 1
            Assert.IsTrue(AdvancedIntegerMath.IsPrime(17)); // 2^4 + 1
            Assert.IsTrue(AdvancedIntegerMath.IsPrime(257)); // 2^8 + 1
            Assert.IsTrue(AdvancedIntegerMath.IsPrime(65537)); // 2^16 + 1

            // Euler's nonprime
            Assert.IsFalse(AdvancedIntegerMath.IsPrime(1000009)); // 1000009 = 293 * 3413

        }

        [TestMethod]
        public void PrimeCounts () {
            Assert.IsTrue(CountPrimes(10) == 4);
            Assert.IsTrue(CountPrimes(100) == 25);
            Assert.IsTrue(CountPrimes(1000) == 168);
            Assert.IsTrue(CountPrimes(10000) == 1229);
        }

        private int CountPrimes (int n) {
            int count = 0;
            for (int k = 1; k <= n; k++) {
                if (AdvancedIntegerMath.IsPrime(k)) count++;
            }
            return (count);
        }

        // Binomial coefficients

        [TestMethod]
        public void BinomialCoefficientRecurrence () {
            // this is the Pascal triangle recurrence
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (int m in TestUtilities.GenerateUniformIntegerValues(0, n-1, 5)) {
                    Console.WriteLine("n = {0}, m = {1}", n, m);
                    Console.WriteLine(AdvancedIntegerMath.BinomialCoefficient(n + 1, m + 1));
                    Console.WriteLine(AdvancedIntegerMath.BinomialCoefficient(n, m + 1));
                    Console.WriteLine(AdvancedIntegerMath.BinomialCoefficient(n, m));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedIntegerMath.BinomialCoefficient(n + 1, m + 1),
                        AdvancedIntegerMath.BinomialCoefficient(n, m) + AdvancedIntegerMath.BinomialCoefficient(n, m + 1)
                    ));
                }
            }
        }

        [TestMethod]
        public void BinomialCoefficientInequality () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 6)) {
                foreach (int m in TestUtilities.GenerateUniformIntegerValues(0, n, 4)) {
                    double lower = Math.Pow(1.0 * n / m, m);
                    double upper = Math.Pow(Math.E * n / m, m);
                    double value = AdvancedIntegerMath.BinomialCoefficient(n, m);
                    Assert.IsTrue((lower <= value) && (value <= upper));
                }
            }
        }

        [TestMethod]
        public void BinomialCoefficientSpecialCases () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 6)) {
                Assert.IsTrue(AdvancedIntegerMath.BinomialCoefficient(n, 0) == 1);
                Assert.IsTrue(AdvancedIntegerMath.BinomialCoefficient(n, 1) == n);
                Assert.IsTrue(AdvancedIntegerMath.BinomialCoefficient(n, n - 1) == n);
                Assert.IsTrue(AdvancedIntegerMath.BinomialCoefficient(n, n) == 1);
            }
        }

        [TestMethod]
        public void BinomialCoefficientSums () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 8)) {
                double S0 = 0.0;
                double S1 = 0.0;
                double S2 = 0.0;
                int k = 0;
                foreach (double B in AdvancedIntegerMath.BinomialCoefficients(n)) {
                    S0 += B;
                    S1 += k * B;
                    S2 += k * k * B;
                    k += 1;
                }
                Assert.IsTrue(TestUtilities.IsNearlyEqual(S0, MoreMath.Pow(2, n)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(S1, MoreMath.Pow(2, n - 1) * n));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(S2, MoreMath.Pow(2, n - 2) * n * (n + 1)));
            }
        }

        [TestMethod]
        public void BinomialCoefficientAgreement () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                int m = -1;
                IEnumerator B = AdvancedIntegerMath.BinomialCoefficients(n).GetEnumerator();
                while (B.MoveNext()) {
                    m++;
                    Assert.IsTrue(TestUtilities.IsNearlyEqual((double)B.Current, AdvancedIntegerMath.BinomialCoefficient(n, m)));
                }

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
        public void PowModTest () {
            // an example from http://en.wikipedia.org/wiki/Modular_exponentiation
            Assert.IsTrue(AdvancedIntegerMath.PowMod(4, 13, 497) == 445);
            // from the calculator http://www.math.umn.edu/~garrett/crypto/a01/FastPow.html
            Assert.IsTrue(AdvancedIntegerMath.PowMod(130000000, 670, 59) == 3);
            // some examples from http://modular.math.washington.edu/edu/2007/spring/ent/ent-html/node81.html
            Assert.IsTrue(AdvancedIntegerMath.PowMod(2, 60, 779167) == 710981);
            Assert.IsTrue(AdvancedIntegerMath.PowMod(2, 360360, 779167) == 584877);
        }

        [TestMethod]
        public void FibonacciSpecialCases () {
            Assert.IsTrue(AdvancedIntegerMath.FibonacciNumber(0) == 0.0);
            Assert.IsTrue(AdvancedIntegerMath.FibonacciNumber(1) == 1.0);
            Assert.IsTrue(AdvancedIntegerMath.FibonacciNumber(2) == 1.0);
            Assert.IsTrue(AdvancedIntegerMath.FibonacciNumber(3) == 2.0);
        }

        [TestMethod]
        public void FibonacciRecurrence () {
            foreach (int n in TestUtilities.GenerateIntegerValues(2, 1000, 8)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedIntegerMath.FibonacciNumber(n),
                    AdvancedIntegerMath.FibonacciNumber(n - 1) + AdvancedIntegerMath.FibonacciNumber(n - 2)
                ));
            }
        }

        [TestMethod]
        public void BellNumberAsPoissonMoment () {

            // The Bell numbers are the moments of the Poisson distribution
            // with \mu = 1.

            UnivariateDistribution d = new PoissonDistribution(1.0);

            foreach (int r in TestUtilities.GenerateIntegerValues(1, 50, 4)) {

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedIntegerMath.BellNumber(r),
                    d.RawMoment(r)
                ));
            }

        }

        [TestMethod]
        public void BellNumberRecurrence () {

            // B_{n+1} = \sum_{k=0}^{n} {n \choose k} B_k

            foreach (int r in TestUtilities.GenerateIntegerValues(50, 100, 2)) {
                double s = 0.0;
                IEnumerator<double> b = AdvancedIntegerMath.BinomialCoefficients(r).GetEnumerator();
                for (int k = 0; k <= r; k++) {
                    b.MoveNext();
                    s += b.Current * AdvancedIntegerMath.BellNumber(k);
                }
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s, AdvancedIntegerMath.BellNumber(r + 1)));
            }

        }

        [TestMethod]
        public void StirlingNumber1SpecialCases () {

            Assert.IsTrue(AdvancedIntegerMath.StirlingNumber1(0, 0) == 1);
            Assert.IsTrue(AdvancedIntegerMath.StirlingNumber1(1, 0) == 0);
            Assert.IsTrue(AdvancedIntegerMath.StirlingNumber1(1, 1) == 1);
            Assert.IsTrue(AdvancedIntegerMath.StirlingNumber1(2, 0) == 0);
            Assert.IsTrue(AdvancedIntegerMath.StirlingNumber1(2, 1) == 1);
            Assert.IsTrue(AdvancedIntegerMath.StirlingNumber1(2, 2) == 1);

            foreach (int m in TestUtilities.GenerateIntegerValues(2, 100, 4)) {

                Assert.IsTrue(AdvancedIntegerMath.StirlingNumber1(m, 0) == 0);
                Assert.IsTrue(AdvancedIntegerMath.StirlingNumber1(m, 1) == AdvancedIntegerMath.Factorial(m - 1));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedIntegerMath.StirlingNumber1(m, 2),
                    AdvancedIntegerMath.Factorial(m - 1) * AdvancedIntegerMath.HarmonicNumber(m - 1)
                ));

                Assert.IsTrue(AdvancedIntegerMath.StirlingNumber1(m, m) == 1);
                Assert.IsTrue(AdvancedIntegerMath.StirlingNumber1(m, m - 1) == AdvancedIntegerMath.BinomialCoefficient(m, 2));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedIntegerMath.StirlingNumber1(m, m - 2),
                    (3 * m - 1) * AdvancedIntegerMath.BinomialCoefficient(m, 3) / 4
                ));

            }

        }

        [TestMethod]
        public void StirlingNumber2SpecialCases () {

            Assert.IsTrue(AdvancedIntegerMath.StirlingNumber2(0, 0) == 1);
            Assert.IsTrue(AdvancedIntegerMath.StirlingNumber2(1, 0) == 0);
            Assert.IsTrue(AdvancedIntegerMath.StirlingNumber2(1, 1) == 1);
            Assert.IsTrue(AdvancedIntegerMath.StirlingNumber2(2, 0) == 0);
            Assert.IsTrue(AdvancedIntegerMath.StirlingNumber2(2, 1) == 1);
            Assert.IsTrue(AdvancedIntegerMath.StirlingNumber2(2, 2) == 1);
            Assert.IsTrue(AdvancedIntegerMath.StirlingNumber2(3, 0) == 0);
            Assert.IsTrue(AdvancedIntegerMath.StirlingNumber2(3, 1) == 1);
            Assert.IsTrue(AdvancedIntegerMath.StirlingNumber2(3, 2) == 3);
            Assert.IsTrue(AdvancedIntegerMath.StirlingNumber2(3, 3) == 1);

            foreach (int m in TestUtilities.GenerateIntegerValues(2, 100, 4)) {

                Assert.IsTrue(AdvancedIntegerMath.StirlingNumber2(m, 0) == 0);
                Assert.IsTrue(AdvancedIntegerMath.StirlingNumber2(m, 1) == 1);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedIntegerMath.StirlingNumber2(m, 2),
                    MoreMath.Pow(2, m - 1) - 1
                ));

                Assert.IsTrue(AdvancedIntegerMath.StirlingNumber2(m, m) == 1);
                Assert.IsTrue(AdvancedIntegerMath.StirlingNumber2(m, m - 1) == AdvancedIntegerMath.BinomialCoefficient(m, 2));

            }

        }

        [TestMethod]
        public void StirlingNumberInequality () {

            // Rows are log-concave

            foreach (int n in TestUtilities.GenerateIntegerValues(10, 100, 4)) {
                foreach (int k in TestUtilities.GenerateRealValues(2, n-1, 4)) {

                    Assert.IsTrue(
                        MoreMath.Sqr(AdvancedIntegerMath.StirlingNumber1(n, k)) >=
                        AdvancedIntegerMath.StirlingNumber1(n, k - 1) *
                        AdvancedIntegerMath.StirlingNumber1(n, k + 1)
                    );

                    Assert.IsTrue(
                        MoreMath.Sqr(AdvancedIntegerMath.StirlingNumber2(n, k)) >=
                        AdvancedIntegerMath.StirlingNumber2(n, k - 1) *
                        AdvancedIntegerMath.StirlingNumber2(n, k + 1)
                    );

                }
            }

        }

        [TestMethod]
        public void StirlingNumbers1RowSum () {

            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {

                double sum = 0.0;
                foreach (double s in AdvancedIntegerMath.StirlingNumbers1(n)) {
                    sum += s;
                }
                Assert.IsTrue(TestUtilities.IsNearlyEqual(sum, AdvancedIntegerMath.Factorial(n)));

            }

        }

        [TestMethod]
        public void StirlingNumbers2RowSum () {

            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {

                double sum = 0.0;
                foreach (double s in AdvancedIntegerMath.StirlingNumbers2(n)) {
                    sum += s;
                }
                Assert.IsTrue(TestUtilities.IsNearlyEqual(sum, AdvancedIntegerMath.BellNumber(n)));

            }

        }

        [TestMethod]
        public void StirlingNumbers2ColumnSum () {

            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (int k in TestUtilities.GenerateIntegerValues(1, n, 4)) {
                    
                    double sum = 0.0;
                    for (int m = k; m <= n; m++) {
                        sum +=
                            AdvancedIntegerMath.BinomialCoefficient(n, m) *
                            AdvancedIntegerMath.StirlingNumber2(m, k);
                    }

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(sum, AdvancedIntegerMath.StirlingNumber2(n + 1, k + 1)));

                }
            }

        }

        [TestMethod]
        public void StirlingNumber2Parity () {

            // Sterling numbers of the 2nd kind provably have the same parity
            // as a particular related binomial coefficient.

            // This is actually quite a difficult test because it tests the least
            // significant bits of the numbers. We can't operate on values
            // with more than 52 bits, since for those numbers we will loose
            // the least significant bits.

            foreach (int n in TestUtilities.GenerateIntegerValues(1, 50, 4)) {
                foreach (int k in TestUtilities.GenerateIntegerValues(1, n, 4)) {

                    long B = (long) AdvancedIntegerMath.BinomialCoefficient(n - k / 2 - 1, n - k);
                    long S = (long) AdvancedIntegerMath.StirlingNumber2(n, k);

                    Assert.IsTrue(B % 2 == S % 2);

                }
            }

        }

        [TestMethod]
        public void StirlingNumberMatrixInverse () {

            int n = 8;

            SquareMatrix S1 = new SquareMatrix(n);
            for (int i = 0; i < n; i++) {
                double[] s = AdvancedIntegerMath.StirlingNumbers1(i);
                for (int j = 0; j < s.Length; j++) {
                    if ((i - j) % 2 == 0) {
                        S1[i, j] = s[j];
                    } else {
                        S1[i, j] = -s[j];
                    }
                }
            }

            SquareMatrix S2 = new SquareMatrix(n);
            for (int i = 0; i < n; i++) {
                double[] s = AdvancedIntegerMath.StirlingNumbers2(i);
                for (int j = 0; j < s.Length; j++) {
                    S2[i, j] = s[j];
                }
            }

            SquareMatrix S12 = S1 * S2;

            SquareMatrix I = new SquareMatrix(n);
            for (int i = 0; i < n; i++) I[i, i] = 1.0;

            Assert.IsTrue(TestUtilities.IsNearlyEqual(S12, I));

        }

        // http://www.cs.columbia.edu/~cs4205/files/CM5.pdf
        // contains many lesser-known relationships among Stirling numbers

        [TestMethod]
        public void BernoulliNumbers () {
            Assert.IsTrue(AdvancedIntegerMath.BernoulliNumber(0) == 1.0);
            Assert.IsTrue(AdvancedIntegerMath.BernoulliNumber(1) == -1.0 / 2.0);
            Assert.IsTrue(AdvancedIntegerMath.BernoulliNumber(2) == 1.0 / 6.0);
            Assert.IsTrue(AdvancedIntegerMath.BernoulliNumber(3) == 0.0);
        }
    }

}
