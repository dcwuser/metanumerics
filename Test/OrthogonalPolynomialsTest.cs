using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;


namespace Test {



    [TestClass()]
    public class OrthogonalPolynomialsTest {


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


        // Orthogonal Polynomials

        private int[] orders = new int[] { 3, 30, 300 };
        private double[] arguments = new double[] { 0.11, 1.1, 11.1, 111.1 };


        // Hermite

        [TestMethod]
        public void HermiteSpecialCaseTest () {
            foreach (double x in arguments) {
                Assert.IsTrue(OrthogonalPolynomials.HermiteH(0, x) == 1.0);
                Assert.IsTrue(OrthogonalPolynomials.HermiteH(1, x) == 2.0 * x);
            }
            foreach (int n in orders) {
                if (n > 100) continue;
                if (n % 2 == 0) {
                    int m = n / 2;
                    int s = 1;
                    if (m % 2 != 0) s = -s;
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.HermiteH(n, 0.0), s * AdvancedIntegerMath.Factorial(n) / AdvancedIntegerMath.Factorial(m)));
                } else {
                    Assert.IsTrue(OrthogonalPolynomials.HermiteH(n, 0.0) == 0.0);
                }
            }
        }

        [TestMethod]
        public void HermiteReflectionTest () {
            foreach (int n in orders) {
                foreach (double x in arguments) {
                    if (n % 2 == 0) {
                        Assert.IsTrue(OrthogonalPolynomials.HermiteH(n, -x) == OrthogonalPolynomials.HermiteH(n, x));
                    } else {
                        Assert.IsTrue(OrthogonalPolynomials.HermiteH(n, -x) == -OrthogonalPolynomials.HermiteH(n, x));
                    }
                }
            }
        }

        [TestMethod]
        public void HermiteRecurrenceTest () {
            foreach (int n in orders) {
                if (n < 100) {
                    foreach (double x in arguments) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.HermiteH(n + 1, x), 2.0 * x * OrthogonalPolynomials.HermiteH(n, x) - 2.0 * n * OrthogonalPolynomials.HermiteH(n - 1, x)), String.Format("n={0}, x={1}, H={2}", n, x, OrthogonalPolynomials.HermiteH(n, x)));
                    }
                }
            }
        }

        [TestMethod]
        public void HermiteAdditionTest () {
            foreach (int n in orders) {
                if (n > 100) continue;
                foreach (double x in arguments) {
                    foreach (double y in arguments) {
                        double value = OrthogonalPolynomials.HermiteH(n, x + y);
                        double sum = 0.0;
                        for (int k = 0; k <= n; k++) {
                            sum += AdvancedIntegerMath.BinomialCoefficient(n, k) * OrthogonalPolynomials.HermiteH(k, x) * Math.Pow(2 * y, n-k);
                        }
                        Console.WriteLine("n={0}, x={1}, y={2}, H(x+y)={3}, sum={4}", n, x, y, value, sum);
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(value, sum));
                    }
                }
            }
        }

        [TestMethod]
        public void HermiteSumTest () {
            // accuracy of equality falls at high n because of cancelations among large terms
            for (int n=0; n<12; n++) {
                double sum = 0.0;
                for (int k = 0; k <= n; k++) {
                    double term = AdvancedIntegerMath.BinomialCoefficient(n, k) * OrthogonalPolynomials.HermiteH(n, k);
                    Console.WriteLine("  k={0}, term={1}", k, term);
                    if ((n - k) % 2 == 0) {
                        sum += term;
                    } else {
                        sum -= term;
                    }
                }
                Console.WriteLine("n={0}, sum={1}", n, sum);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(sum, Math.Pow(2, n) * AdvancedIntegerMath.Factorial(n)));
            }
        }

        [TestMethod]
        public void HermiteOrthonormalityTest () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 30, 3)) {
                foreach (int m in TestUtilities.GenerateIntegerValues(1, 30, 3)) {
                    Function<double, double> f = delegate(double x) {
                        return (Math.Exp(-x * x) * OrthogonalPolynomials.HermiteH(n, x) * OrthogonalPolynomials.HermiteH(m, x));
                    };
                    Interval r = Interval.FromEndpoints(Double.NegativeInfinity, Double.PositiveInfinity);
                    double I = FunctionMath.Integrate(f, r);
                    double N = Math.Sqrt(Math.PI) * Math.Pow(2.0, n) * AdvancedIntegerMath.Factorial(n);
                    Console.WriteLine("{0} {1} {2} {3}", n, m, I, N);
                    if (n == m) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(I, N));
                    } else {
                        Assert.IsTrue(Math.Abs(I) <= TestUtilities.TargetPrecision);
                    }
                }
            }
        }

        // test HermiteHe orthonormality, too

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void HermiteInvalidOrderTest () {
            OrthogonalPolynomials.HermiteH(-1, 1.0);
        }

        // Laguerre

        [TestMethod]
        public void LaguerreSpecialCaseTest () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                Assert.IsTrue(OrthogonalPolynomials.LaguerreL(n, 0.0) == 1.0);
            }
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 5)) {
                Assert.IsTrue(OrthogonalPolynomials.LaguerreL(0, x) == 1.0);
            }
        }

        [TestMethod]
        public void LaguerreRecurrenceTest () {
            foreach (int n in TestUtilities.GenerateRealValues(1.0, 1.0E2, 5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E4,10)) {
                    double LP = OrthogonalPolynomials.LaguerreL(n + 1, x);
                    double L = OrthogonalPolynomials.LaguerreL(n, x);
                    double LM = OrthogonalPolynomials.LaguerreL(n - 1, x);
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual((n + 1) * LP, n * LM, (2 * n + 1 - x) * L));
                }
            }
        }

        [TestMethod]
        public void LaguerreInequalityTest () {
            foreach (int n in TestUtilities.GenerateRealValues(1.0, 1.0E2, 5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 5)) {
                    Assert.IsTrue(OrthogonalPolynomials.LaguerreL(n,x) <= Math.Exp(x/2.0));
                }
            }
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void LaguerreInvalidOrderTest () {
            OrthogonalPolynomials.LaguerreL(-1, 1.0);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void LaguerreInvalidArgumentTest () {
            OrthogonalPolynomials.LaguerreL(1, -0.1);
        }

        [TestMethod]
        public void LaguerreOrthonormalityTest () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 3)) {
                foreach (int m in TestUtilities.GenerateIntegerValues(1, 100, 3)) {
                    Function<double, double> f = delegate(double x) {
                        return (Math.Exp(-x) * OrthogonalPolynomials.LaguerreL(n, x) * OrthogonalPolynomials.LaguerreL(m, x));
                    };
                    Interval r = Interval.FromEndpoints(0.0, Double.PositiveInfinity);
                    double I = FunctionMath.Integrate(f, r);
                    Console.WriteLine("{0} {1} {2}", n, m, I);
                    if (n == m) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(I, 1.0));
                    } else {
                        Assert.IsTrue(Math.Abs(I) < TestUtilities.TargetPrecision);
                    }
                }
            }
        }



        // Chebyshev

        private double[] cArguments = new double[] { 0.00001, 0.02, 0.34, 0.56, 0.78, 0.99999 };

        [TestMethod]
        public void ChebyshevSpecialCaseTest () {
            foreach (double x in cArguments) {
                Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(0, x) == 1.0);
                Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(1, x) == x);
            }
            foreach (int n in orders) {
                if (n % 2 == 0) {
                    if (n % 4 == 0) {
                        Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(n, 0.0) == 1.0);
                    } else {
                        Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(n, 0.0) == -1.0);
                    }
                } else {
                    Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(n, 0.0) == 0.0);
                }
                Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(n, 1.0) == 1.0);
            }
        }

        [TestMethod]
        public void ChebyshevInequalityTest () {
            foreach (int n in orders) {
                foreach (double x in cArguments) {
                    Assert.IsTrue(Math.Abs(OrthogonalPolynomials.ChebyshevT(n, x)) <= 1.0);
                }
            }
        }

        [TestMethod]
        public void ChebyshevRecurrenceTest () {
            foreach (int n in orders) {
                foreach (double x in cArguments) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ChebyshevT(n + 1, x), 2.0 * x * OrthogonalPolynomials.ChebyshevT(n, x) - OrthogonalPolynomials.ChebyshevT(n - 1, x)));
                }
            }
        }

        [TestMethod]
        public void ChebyshevReflectionTest () {
            foreach (int n in orders) {
                foreach (double x in cArguments) {
                    if (n % 2 == 0) {
                        Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(n, -x) == OrthogonalPolynomials.ChebyshevT(n, x));
                    } else {
                        Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(n, -x) == -OrthogonalPolynomials.ChebyshevT(n, x));
                    }
                }
            }
        }

        [TestMethod]
        public void ChebyshevDoublingTest () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1,100,5)) {
                foreach (double x in cArguments) {
                    Console.WriteLine("n={0}, x={1}, T2={2}, T1={3}", n, x, OrthogonalPolynomials.ChebyshevT(2 * n, x), OrthogonalPolynomials.ChebyshevT(n, 1.0 - 2 * x * x));
                    if (n % 2 == 0) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ChebyshevT(2 * n, x), OrthogonalPolynomials.ChebyshevT(n, 1.0 - 2 * x * x)));
                    } else {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ChebyshevT(2 * n, x), -OrthogonalPolynomials.ChebyshevT(n, 1.0 - 2 * x * x)));
                    }
                    //Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ChebyshevT(2 * n, Math.Sqrt((x + 1.0) / 2.0)), OrthogonalPolynomials.ChebyshevT(n, x)));
                }
            }
        }

        [TestMethod]
        public void ChebyshevMultiplicationTest () {
            foreach (int n in orders) {
                foreach (int m in orders) {
                    if (n * m > 500) continue; // we don't test very high order
                    foreach (double x in cArguments) {
                        Console.WriteLine("n={0}, m={1}, x={2}", n, m, x);
                        double Tm = OrthogonalPolynomials.ChebyshevT(m, x);
                        Console.WriteLine("T_m(x)= {0}", Tm);
                        double Tn = OrthogonalPolynomials.ChebyshevT(n, Tm);
                        Console.WriteLine("T_n(T_m) = {0}", Tn);
                        double Tnm = OrthogonalPolynomials.ChebyshevT(n * m, x);
                        Console.WriteLine("T_nm(x) = {0}", Tnm);
                        double d = Math.Abs(Tn - Tnm) / Math.Abs(Tnm);
                        Console.WriteLine("relative d = {0}", d);
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ChebyshevT(n, OrthogonalPolynomials.ChebyshevT(m, x)), OrthogonalPolynomials.ChebyshevT(n * m, x)));
                    }
                }
            }
        }

        [TestMethod]
        public void ChebyshevCosineTest () {
            // test T_n(cos t) = cos(n t)
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double t in GenerateRandomAngles(-Math.PI, Math.PI, 5)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ChebyshevT(n, Math.Cos(t)), Math.Cos(n * t)));
                }
            }
        }


        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void ChebyshevInvalidArgumentTest () {
            OrthogonalPolynomials.ChebyshevT(2, -1.1);
        }

        // orthonormality test for chebyshev fails because endpoint singularity of weight causes
        // the numerical integral not to converge


        // Legendre

        [TestMethod]
        public void LegendreSpecialCaseTest () {
            foreach (double x in cArguments) {
                Assert.IsTrue(OrthogonalPolynomials.LegendreP(0, x) == 1.0);
                Assert.IsTrue(OrthogonalPolynomials.LegendreP(1, x) == x);
            }
            foreach (int n in orders) {
                Assert.IsTrue(OrthogonalPolynomials.LegendreP(n, 1.0) == 1.0);
            }
        }

        [TestMethod]
        public void LegendreInequalityTest () {
            foreach (int n in orders) {
                foreach (double x in cArguments) {
                    Assert.IsTrue(Math.Abs(OrthogonalPolynomials.LegendreP(n, x)) <= 1.0);
                }
            }
        }

        [TestMethod]
        public void LegendreRecurrenceTest () {
            foreach (int n in orders) {
                foreach (double x in cArguments) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual((n + 1) * OrthogonalPolynomials.LegendreP(n + 1, x), (2 * n + 1) * x * OrthogonalPolynomials.LegendreP(n, x) - n * OrthogonalPolynomials.LegendreP(n - 1, x)));
                }
            }
        }

        [TestMethod]
        public void LegendreReflectionTest () {
            foreach (int n in orders) {
                foreach (double x in cArguments) {
                    if (n % 2 == 0) {
                        Assert.IsTrue(OrthogonalPolynomials.LegendreP(n, -x) == OrthogonalPolynomials.LegendreP(n, x));
                    } else {
                        Assert.IsTrue(OrthogonalPolynomials.LegendreP(n, -x) == -OrthogonalPolynomials.LegendreP(n, x));
                    }
                }
            }
        }

        [TestMethod]
        public void LegendreNegativeOrderTest () {
            foreach (int n in orders) {
                foreach (double x in cArguments) {
                    Assert.IsTrue(OrthogonalPolynomials.LegendreP(-n, x) == OrthogonalPolynomials.LegendreP(n - 1, x));
                }
            }

        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void LegendreInvalidArgumentTest () {
            OrthogonalPolynomials.LegendreP(2, -1.1);
        }

        [TestMethod]
        public void LegendreOrthonormalityTest () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 3)) {
                foreach (int m in TestUtilities.GenerateIntegerValues(1, 100, 3)) {
                    Function<double, double> f = delegate(double x) {
                        return (OrthogonalPolynomials.LegendreP(n, x) * OrthogonalPolynomials.LegendreP(m, x));
                    };
                    Interval r = Interval.FromEndpoints(-1.0, 1.0);
                    double I = FunctionMath.Integrate(f, r);
                    Console.WriteLine("{0} {1} {2}", n, m, I);
                    if (n == m) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(I, 2.0 / (2 * n + 1)));
                    } else {
                        Assert.IsTrue(Math.Abs(I) < TestUtilities.TargetPrecision);
                    }
                }
            }
        }

        // spherical harmonics

        private double[] GenerateRandomAngles (double lower, double upper, int n) {
            double[] x = new double[n];
            Random rng = new Random(1);
            for (int i = 0; i < n; i++) {
                x[i] = lower + rng.NextDouble() * (upper - lower);
            }
            return (x);
        }

        [TestMethod]
        public void SphericalHarmonicSpecialCaseTest () {
            // ell=0, m=0
            foreach (double theta in GenerateRandomAngles(-Math.PI, Math.PI, 5)) {
                foreach (double phi in GenerateRandomAngles(0.0, 2.0 * Math.PI, 5)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.SphericalHarmonic(0, 0, theta, phi), 1.0 / Math.Sqrt(4.0 * Math.PI)));
                }
            }
            // theta=0, m=0
            foreach (int ell in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double phi in GenerateRandomAngles(0.0, 2.0 * Math.PI, 5)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.SphericalHarmonic(ell, 0, 0.0, phi), Math.Sqrt((2 * ell + 1) / (4.0 * Math.PI))));
                }
            }

        }

        [TestMethod]
        public void SphericalHarmonicLegendreTest () {
            foreach (int ell in TestUtilities.GenerateIntegerValues(1, 100, 6)) {
                foreach (double theta in GenerateRandomAngles(-Math.PI, Math.PI, 5)) {
                    // this loop shouldn't matter since the answer is phi-independent
                    foreach (double phi in GenerateRandomAngles(0.0, 2.0 * Math.PI, 3)) {
                        Complex Y = AdvancedMath.SphericalHarmonic(ell, 0, theta, phi);
                        double LP = Math.Sqrt((2 * ell + 1) / (4.0 * Math.PI)) * OrthogonalPolynomials.LegendreP(ell, Math.Cos(theta));
                        Console.WriteLine("l={0}, m=0, t={1}, p={2}, Y={3}, LP={4}", ell, theta, phi, Y, LP);
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(Y, LP));
                    }
                }
            }
        }

        [TestMethod]
        public void SphericalHarmonicConjugationTest () {
            foreach (int ell in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double theta in GenerateRandomAngles(-Math.PI, Math.PI, 5)) {
                    foreach (double phi in GenerateRandomAngles(0.0, 2.0 * Math.PI, 5)) {
                        for (int m = 0; m <= ell; m++) {
                            int s = 1 - 2 * (m % 2);
                            Complex YPM = AdvancedMath.SphericalHarmonic(ell, m, theta, phi);
                            Complex YNM = s * AdvancedMath.SphericalHarmonic(ell, -m, theta, phi).Conjugate;
                            Console.WriteLine("l={0} m={1} t={2} p={3} Y1={4} Y2={5}", ell, m, theta, phi, YPM, YNM);
                            Assert.IsTrue(YPM == YNM);
                        }
                    }
                }
            }
        }

        [TestMethod]
        public void SphericalHarmonicNormalizationTest () {
            foreach (int ell in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double theta in GenerateRandomAngles(-Math.PI, Math.PI, 5)) {
                    foreach (double phi in GenerateRandomAngles(0.0, 2.0 * Math.PI, 5)) {
                        double sum = 0.0;
                        for (int m = -ell; m <= ell; m++) {
                            Complex Y = AdvancedMath.SphericalHarmonic(ell, m, theta, phi);
                            sum += Y.Re * Y.Re + Y.Im * Y.Im;
                        }
                        Console.WriteLine("l={0}, t={1}, p={2}, sum={3}", ell, theta, phi, sum);
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(sum, (2 * ell + 1) / (4.0 * Math.PI)));
                    }
                }
            }
        }

        [TestMethod]
        public void SphericalHarmonicAdditionTest () {
            foreach (int ell in TestUtilities.GenerateIntegerValues(1, 100, 6)) {
                foreach (double theta in GenerateRandomAngles(-Math.PI, Math.PI, 4)) {
                    foreach (double phi in GenerateRandomAngles(0, 2.0 * Math.PI, 3)) {
                        Complex sum = 0.0;
                        for (int m = -ell; m <= ell; m++) {
                            sum += AdvancedMath.SphericalHarmonic(ell, m, theta, phi) *
                                AdvancedMath.SphericalHarmonic(ell, m, 0.0, phi).Conjugate;
                        }
                        sum = 4.0 * Math.PI / (2 * ell + 1) * sum;
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(sum, OrthogonalPolynomials.LegendreP(ell, Math.Cos(theta))));
                    }
                }
            }
        }

    }
}
