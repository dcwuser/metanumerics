using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using Meta.Numerics.Statistics.Distributions;


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



        // For all orthogonal polynomials, we should test
        //                     Legendre  Chebyshev  Hermine  Laguerre
        //   Orthonormality
        //   Recurrsion
        //   Special cases
        //   Symmetry
        //   Inequality
        // and in some cases other relations are known.

        // Orthogonal Polynomials

        // In order to ensure that very small and very large arguments are tested, or in the
        // bounded case arguments very close to 1, we specify arguments to use
        private int[] orders = new int[] { 3, 4, 33, 44, 333, 444 };
        private double[] aArguments = new double[] { 0.00011, 0.11, 1.1, 11.1, 111.1, 1111.1 };
        private double[] bArguments = new double[] { 0.00001, 0.02, 0.34, 0.56, 0.78, 0.99999 };

        // Hermite

        [TestMethod]
        public void HermiteSpecialCases () {
            foreach (double x in aArguments) {
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
        public void HermiteReflection () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 10)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0, 1000.0, 10)) {
                    Console.WriteLine("n={0} x={1}", n, x);
                    if (n % 2 == 0) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.HermiteH(n, -x), OrthogonalPolynomials.HermiteH(n, x)));
                    } else {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.HermiteH(n, -x), -OrthogonalPolynomials.HermiteH(n, x)));
                    }
                }
            }
        }

        [TestMethod]
        public void HermiteRecurrence () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 10)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0, 1000.0, 10)) {
                    Console.WriteLine("n={0} x={1}", n, x);
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        OrthogonalPolynomials.HermiteH(n + 1, x), 2.0 * n * OrthogonalPolynomials.HermiteH(n - 1, x),
                        2.0 * x * OrthogonalPolynomials.HermiteH(n, x)
                    ));
                }
            }
        }

        [TestMethod]
        public void HermiteAddition () {
            foreach (int n in orders) {
                if (n > 100) continue;
                foreach (double x in aArguments) {
                    foreach (double y in aArguments) {
                        double value = OrthogonalPolynomials.HermiteH(n, x + y);
                        double[] terms = new double[n + 1]; double sum = 0.0;
                        for (int k = 0; k <= n; k++) {
                            terms[k] = AdvancedIntegerMath.BinomialCoefficient(n, k) * OrthogonalPolynomials.HermiteH(k, x) * Math.Pow(2 * y, n-k);
                            sum += terms[k];
                        }
                        Console.WriteLine("n={0}, x={1}, y={2}, H(x+y)={3}, sum={4}", n, x, y, value, sum);
                        Assert.IsTrue(TestUtilities.IsSumNearlyEqual(terms, value));
                    }
                }
            }
        }

        [TestMethod]
        public void HermiteSum () {
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
        public void HermiteOrthonormality () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 30, 3)) {
                foreach (int m in TestUtilities.GenerateIntegerValues(1, 30, 3)) {
                    Func<double, double> f = delegate(double x) {
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
        public void HermiteInvalidOrder () {
            OrthogonalPolynomials.HermiteH(-1, 1.0);
        }


        [TestMethod]
        public void HermiteNormalExpectation () {

            // Wikipedia article on Hermite polynomials
            // X ~ N(\mu, 1) => E(He_n(X)) = \mu^n

            NormalDistribution d = new NormalDistribution(2.0, 1.0);

            foreach(int n in TestUtilities.GenerateIntegerValues(1, 32, 7)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    d.ExpectationValue(x => OrthogonalPolynomials.HermiteHe(n, x)),
                    MoreMath.Pow(d.Mean, n)
                ));       
            }

        }

        // Laguerre

        [TestMethod]
        public void LaguerreSpecialCases () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                Assert.IsTrue(OrthogonalPolynomials.LaguerreL(n, 0.0) == 1.0);
            }
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 5)) {
                Assert.IsTrue(OrthogonalPolynomials.LaguerreL(0, x) == 1.0);
            }
        }

        [TestMethod]
        public void LaguerreRecurrence () {
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
        public void LaguerreInequality () {
            foreach (int n in TestUtilities.GenerateRealValues(1.0, 1.0E2, 5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 5)) {
                    Assert.IsTrue(OrthogonalPolynomials.LaguerreL(n,x) <= Math.Exp(x/2.0));
                }
            }
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void LaguerreInvalidOrder () {
            OrthogonalPolynomials.LaguerreL(-1, 1.0);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void LaguerreInvalidArgument () {
            OrthogonalPolynomials.LaguerreL(1, -0.1);
        }

        [TestMethod]
        public void LaguerreOrthonormality () {
            // to avoid highly oscilatory integrals, don't use very high orders
            int[] orders = TestUtilities.GenerateIntegerValues(1, 30, 3);
            for (int i = 0; i < orders.Length; i++) {
                int n = orders[i];
                for (int j = 0; j <= i; j++) {
                    int m = orders[j];

                    Func<double, double> f = delegate(double x) {
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
            /*
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
            */
        }

        // associated Laguerre

        [TestMethod]
        public void AssociatedLaguerreSpecialCases () {
            foreach (double a in TestUtilities.GenerateRealValues(0.01, 100.0, 5)) {
                foreach (double x in TestUtilities.GenerateRealValues(0.01, 100.0, 5)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        OrthogonalPolynomials.LaguerreL(0, a, x), 1.0
                    ));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        OrthogonalPolynomials.LaguerreL(1, a, x), 1.0 + a - x
                    ));
                }
            }
        }

        [TestMethod]
        public void AssociatedLaguerreAlphaRecurrence () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double a in TestUtilities.GenerateRealValues(0.01, 100.0, 5)) {
                    foreach (double x in TestUtilities.GenerateRealValues(0.01, 100.0, 5)) {
                        Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                            OrthogonalPolynomials.LaguerreL(n, a - 1.0, x), OrthogonalPolynomials.LaguerreL(n - 1, a, x),
                            OrthogonalPolynomials.LaguerreL(n, a, x)
                        ));
                    }
                }
            }
        }

        [TestMethod]
        public void AssociatedLaguerreSum () {
            foreach (int n in TestUtilities.GenerateRealValues(1, 100, 5)) {
                foreach (double a in TestUtilities.GenerateRealValues(0.1, 100.0, 5)) {
                    foreach (double x in TestUtilities.GenerateRealValues(0.01, 1000.0, 5)) {

                        Console.WriteLine("n={0}, a={1}, x={2}", n, a, x);

                        List<double> L = new List<double>(n + 1);
                        for (int k = 0; k <= n; k++) {
                            L.Add(OrthogonalPolynomials.LaguerreL(k, a, x));
                        }

                        Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                            L, OrthogonalPolynomials.LaguerreL(n, a + 1.0, x)
                        ));

                    }

                }
            }
        }

        [TestMethod]
        public void AssociatedLaguerreOrthonormality () {
            // don't let orders get too big, or (1) the Gamma function will overflow and (2) our integral will become highly oscilatory
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 10, 3)) {
                foreach (int m in TestUtilities.GenerateIntegerValues(1, 10, 3)) {
                    foreach (double a in TestUtilities.GenerateRealValues(0.1, 10.0, 5)) {

                        //int n = 2;
                        //int m = 4;
                        //double a = 3.5;

                        Console.WriteLine("n={0} m={1} a={2}", n, m, a);

                        // evaluate the orthonormal integral
                        Func<double, double> f = delegate(double x) {
                            return (Math.Pow(x, a) * Math.Exp(-x) *
                                OrthogonalPolynomials.LaguerreL(m, a, x) *
                                OrthogonalPolynomials.LaguerreL(n, a, x)
                            );
                        };
                        Interval r = Interval.FromEndpoints(0.0, Double.PositiveInfinity);

                        // need to loosen default evaluation settings in order to get convergence in some of these cases
                        // seems to have most convergence problems for large a
                        EvaluationSettings e = new EvaluationSettings();
                        e.AbsolutePrecision = TestUtilities.TargetPrecision;
                        e.RelativePrecision = TestUtilities.TargetPrecision;

                        double I = FunctionMath.Integrate(f, r, e).Value;
                        Console.WriteLine(I);
                        
                        // test for orthonormality
                        if (n == m) {
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                                I, AdvancedMath.Gamma(n + a + 1) / AdvancedIntegerMath.Factorial(n)
                            ));
                        } else {
                            Assert.IsTrue(Math.Abs(I) < TestUtilities.TargetPrecision);
                        }

                    }
                }
            }
        }


        // Chebyshev


        [TestMethod]
        public void ChebyshevSpecialCases () {
            foreach (double x in bArguments) {
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
        public void ChebyshevInequality () {
            foreach (int n in orders) {
                foreach (double x in bArguments) {
                    Assert.IsTrue(Math.Abs(OrthogonalPolynomials.ChebyshevT(n, x)) <= 1.0);
                }
            }
        }

        [TestMethod]
        public void ChebyshevRecurrence () {
            foreach (int n in orders) {
                foreach (double x in bArguments) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ChebyshevT(n + 1, x), 2.0 * x * OrthogonalPolynomials.ChebyshevT(n, x) - OrthogonalPolynomials.ChebyshevT(n - 1, x)));
                }
            }
        }

        [TestMethod]
        public void ChebyshevReflection () {
            foreach (int n in orders) {
                foreach (double x in bArguments) {
                    if (n % 2 == 0) {
                        Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(n, -x) == OrthogonalPolynomials.ChebyshevT(n, x));
                    } else {
                        Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(n, -x) == -OrthogonalPolynomials.ChebyshevT(n, x));
                    }
                }
            }
        }

        [TestMethod]
        public void ChebyshevDoubling () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1,100,5)) {
                foreach (double x in bArguments) {
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
        public void ChebyshevMultiplication () {

            // Verify T(n, T(m, x)) = T(n m, x)

            // Close to x = \pm 1, T_n(x) oscilates strongly for large n. This means T_n(x)
            // swings quickly between -1 and +1, so a very change in argument can
            // have a very large effect on value. In this region, we can't test this relationship
            // because even an ulp error in T(m, x) could change T(n, T(m, x)) by a significant
            // fraction. 

            // Since zeros are given by
            //   x_k = \cos\left(\frac{(2k+1)\pi}{2n}\right)
            // Distance between zeros is
            //   dx = \sqrt{1-x^2} \pi / n

            foreach (int n in orders) {
                foreach (int m in orders) {
                    foreach (double x in bArguments) {

                        double dx = Math.Sqrt(1.0 - x * x) * Math.PI / (n * m);
                        if (dx < 1.0E-2) continue;

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
        public void ChebyshevCosine () {
            // test T_n(cos t) = cos(n t)
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double t in GenerateRandomAngles(-Math.PI, Math.PI, 5)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ChebyshevT(n, Math.Cos(t)), Math.Cos(n * t)));
                }
            }
        }


        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void ChebyshevInvalidArgument () {
            OrthogonalPolynomials.ChebyshevT(2, -1.1);
        }

        // orthonormality test for chebyshev fails because endpoint singularity of weight causes
        // the numerical integral not to converge


        // Legendre

        [TestMethod]
        public void LegendreSpecialCases () {
            foreach (double x in bArguments) {
                Assert.IsTrue(OrthogonalPolynomials.LegendreP(0, x) == 1.0);
                Assert.IsTrue(OrthogonalPolynomials.LegendreP(1, x) == x);
            }
            foreach (int n in orders) {
                Assert.IsTrue(OrthogonalPolynomials.LegendreP(n, 1.0) == 1.0);
            }
        }

        [TestMethod]
        public void LegendreInequality () {
            foreach (int n in orders) {
                foreach (double x in bArguments) {
                    Assert.IsTrue(Math.Abs(OrthogonalPolynomials.LegendreP(n, x)) <= 1.0);
                }
            }
        }

        [TestMethod]
        public void LegendreTuronInequality () {
            foreach (int n in orders) {
                foreach (double x in bArguments) {
                    Assert.IsTrue(
                        MoreMath.Pow(OrthogonalPolynomials.LegendreP(n, x), 2) >=
                        OrthogonalPolynomials.LegendreP(n - 1, x) * OrthogonalPolynomials.LegendreP(n + 1, x)
                    );
                }
            }
        }

        [TestMethod]
        public void LegendreRecurrence () {
            foreach (int n in orders) {
                foreach (double x in bArguments) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual((n + 1) * OrthogonalPolynomials.LegendreP(n + 1, x), (2 * n + 1) * x * OrthogonalPolynomials.LegendreP(n, x) - n * OrthogonalPolynomials.LegendreP(n - 1, x)));
                }
            }
        }

        [TestMethod]
        public void LegendreReflection () {
            foreach (int n in orders) {
                foreach (double x in bArguments) {
                    if (n % 2 == 0) {
                        Assert.IsTrue(OrthogonalPolynomials.LegendreP(n, -x) == OrthogonalPolynomials.LegendreP(n, x));
                    } else {
                        Assert.IsTrue(OrthogonalPolynomials.LegendreP(n, -x) == -OrthogonalPolynomials.LegendreP(n, x));
                    }
                }
            }
        }

        [TestMethod]
        public void LegendreNegativeOrder () {
            foreach (int n in orders) {
                foreach (double x in bArguments) {
                    Assert.IsTrue(OrthogonalPolynomials.LegendreP(-n, x) == OrthogonalPolynomials.LegendreP(n - 1, x));
                }
            }

        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void LegendreInvalidArgument () {
            OrthogonalPolynomials.LegendreP(2, -1.1);
        }

        [TestMethod]
        public void LegendreOrthonormality () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 3)) {
                foreach (int m in TestUtilities.GenerateIntegerValues(1, 100, 3)) {
                    Func<double, double> f = delegate(double x) {
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
        public void SphericalHarmonicSpecialCases () {
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
        public void SphericalHarmonicLegendre () {
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
        public void SphericalHarmonicConjugation () {
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
        public void SphericalHarmonicNormalization () {
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
        public void SphericalHarmonicAddition () {
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

        [TestMethod]
        public void SphericalHarmonicLowOrders () {
            foreach (double theta in GenerateRandomAngles(-Math.PI, Math.PI, 4)) {
                foreach (double phi in GenerateRandomAngles(0, 2.0 * Math.PI, 4)) {

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.SphericalHarmonic(0, 0, theta, phi), 1.0 / Math.Sqrt(4.0 * Math.PI)));
                    //Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.SphericalHarmonic(1, -1, theta, phi), Math.Sqrt(3.0 / 8.0 / Math.PI) * Math.Sin(theta) * ComplexMath.Exp(-ComplexMath.I * phi)));
                    //Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.SphericalHarmonic(1, 1, theta, phi), -Math.Sqrt(3.0 / 8.0 / Math.PI) * Math.Sin(theta) * ComplexMath.Exp(ComplexMath.I * phi)));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.SphericalHarmonic(1, 0, theta, phi), Math.Sqrt(3.0 / 4.0 / Math.PI) * Math.Cos(theta)));

                }
            }
        }

        [TestMethod]
        public void ZernikeSpecialCases () {

            foreach (double x in TestUtilities.GenerateUniformRealValues(0.0, 1.0, 4)) {
                Console.WriteLine(x);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ZernikeR(0, 0, x), 1.0));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ZernikeR(1, 1, x), x));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ZernikeR(2, 0, x), 2.0 * x * x - 1.0));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ZernikeR(2, 2, x), x * x));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ZernikeR(3, 1, x), 3.0 * x * x * x - 2.0 * x));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ZernikeR(3, 3, x), x * x * x));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ZernikeR(4, 0, x), 6.0 * x * x * x * x - 6.0 * x * x + 1.0));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ZernikeR(4, 2, x), 4.0 * x * x * x * x - 3.0 * x * x));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ZernikeR(4, 4, x), x * x * x * x));
            }

        }

        [TestMethod]
        public void ZernikeOrthonormality () {

            foreach (int na in TestUtilities.GenerateIntegerValues(1, 30, 8)) {
                foreach (int nb in TestUtilities.GenerateIntegerValues(1, 30, 6)) {
                    foreach (int m in TestUtilities.GenerateIntegerValues(1, 30, 4)) {
                        // skip trivial cases
                        if ((m > na) || (m > nb)) continue;
                        if ((m % 2) != (na % 2)) continue;
                        if ((m % 2) != (nb % 2)) continue;

                        Func<double,double> f = delegate (double x) {
                            return (OrthogonalPolynomials.ZernikeR(na, m, x) * OrthogonalPolynomials.ZernikeR(nb, m, x) * x);
                        };

                        double I = FunctionMath.Integrate(f, Interval.FromEndpoints(0.0, 1.0));

                        Console.WriteLine("({0},{1}) ({2},{3}) {4}", na, m, nb, m, I);

                        if (na == nb) {
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(I, 0.5 / (na + 1)));
                        } else {
                            Assert.IsTrue(Math.Abs(I) < TestUtilities.TargetPrecision);
                        }

                    }

                }
            }

        }

        
        [TestMethod]
        public void ZernikeBessel () {

            foreach (int n in TestUtilities.GenerateIntegerValues(1, 30, 6)) {
                foreach (int m in TestUtilities.GenerateUniformIntegerValues(0, n, 6)) {
                    //int m = n;
                    if ((m % 2) != (n % 2)) continue; // skip trivial cases
                    foreach (double x in TestUtilities.GenerateRealValues(1.0, 10.0, 4)) {

                        Console.WriteLine("{0} {1} {2}", n, m, x);

                        Func<double, double> f = delegate(double rho) {
                            return (
                                OrthogonalPolynomials.ZernikeR(n, m, rho) *
                                AdvancedMath.BesselJ(m, x * rho) * rho
                             );
                        };

                        double I = FunctionMath.Integrate(f, Interval.FromEndpoints(0.0, 1.0));

                        double J = AdvancedMath.BesselJ(n + 1, x) / x;
                        if ((n - m) / 2 % 2 != 0) J = -J;

                        Console.WriteLine("  {0} {1}", I, J);

                        // near-zero integrals result from delicate cancelations, and thus
                        // cannot be expected to achieve the relative target precision;
                        // the can achieve the absolute target precisions, so translate
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(I, J, Math.Abs(TestUtilities.TargetPrecision / J)));

                    }
                }
            }

        }

        [TestMethod]
        public void AssociatedLegendreLowOrders () {

            foreach (double x in TestUtilities.GenerateRealValues(0.01, 1.0, 10)) {

                Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.LegendreP(0, 0, x), 1.0));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.LegendreP(1, 0, x), x));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.LegendreP(1, 1, x), -Math.Sqrt(1.0 - x * x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.LegendreP(2, 0, x), (3.0 * x * x - 1.0) / 2.0));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.LegendreP(2, 1, x), -3.0 * x * Math.Sqrt(1.0 - x * x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.LegendreP(2, 2, x), 3.0 * (1.0 - x * x)));

            }

        }

        [TestMethod]
        public void AssociatedLegendreAgreement () {
            foreach (int l in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double x in TestUtilities.GenerateRealValues(0.01, 1.0, 5)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        OrthogonalPolynomials.LegendreP(l, 0, x), OrthogonalPolynomials.LegendreP(l, x)
                    ));
                }
            }
        }

        [TestMethod]
        public void AssociatedLegendreRecurrenceL () {

            foreach (int l in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (int m in TestUtilities.GenerateUniformIntegerValues(0, l-1, 4)) {
                    foreach (double x in TestUtilities.GenerateRealValues(0.01, 1.0, 4)) {
                        Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                            (l-m+1) * OrthogonalPolynomials.LegendreP(l+1, m, x),
                            (l+m) * OrthogonalPolynomials.LegendreP(l-1, m, x),
                            (2*l+1) * x * OrthogonalPolynomials.LegendreP(l, m, x)
                        ));
                    }
                }
            }

        }

        [TestMethod]
        public void AssociatedLegendreRecurrenceM () {

            foreach (int l in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (int m in TestUtilities.GenerateUniformIntegerValues(0, l - 1, 4)) {
                    foreach (double x in TestUtilities.GenerateRealValues(0.1, 1.0, 4)) {
                        Console.WriteLine("l={0} m={1} x={2}", l, m, x);
                        Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                            OrthogonalPolynomials.LegendreP(l, m + 1, x),
                            (l + m) * (l - m + 1 ) * OrthogonalPolynomials.LegendreP(l, m - 1, x),
                            - 2 * m * x / Math.Sqrt(1.0 - x * x) * OrthogonalPolynomials.LegendreP(l, m, x)
                        ));
                    }
                }
            }

        }

        [TestMethod]
        public void AssociatedLegendreOrthonormalityL () {

            // don't let l and m get too big, or numerical integration will fail due to heavily oscilatory behavior

            int[] ells = TestUtilities.GenerateIntegerValues(1, 10, 4);
            for (int ki = 0; ki < ells.Length; ki++) {
                int k = ells[ki];
                for (int li = 0; li <= ki; li++) {
                    int l = ells[li];
                    foreach(int m in TestUtilities.GenerateUniformIntegerValues(0, Math.Min(k,l), 4)) {

                        Func<double, double> f = delegate(double x) {
                            return (OrthogonalPolynomials.LegendreP(k, m, x) * OrthogonalPolynomials.LegendreP(l, m, x));
                        };
                        double I = FunctionMath.Integrate(f, Interval.FromEndpoints(-1.0, 1.0));

                        Console.WriteLine("k={0} l={1} m={2} I={3}", k, l, m, I);

                        if (k == l) {
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                                I, 2.0 / (2 * k + 1) * Math.Exp(AdvancedIntegerMath.LogFactorial(l+m) - AdvancedIntegerMath.LogFactorial(l-m))
                            ));
                        } else {
                            Assert.IsTrue(Math.Abs(I) < TestUtilities.TargetPrecision);
                        }
                    }
                }
            }

        }

    }
}
