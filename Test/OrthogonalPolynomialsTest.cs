using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using Meta.Numerics.Statistics.Distributions;
using System.Runtime.Remoting.Metadata.W3cXsd2001;
using System.Runtime.CompilerServices;

namespace Test {



    [TestClass]
    public class OrthogonalPolynomialsTest {

        // For all orthogonal polynomials, we should test
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
        public void HermiteLowOrders() {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 8)) {
                Assert.IsTrue(OrthogonalPolynomials.HermiteH(0, x) == 1.0);
                Assert.IsTrue(OrthogonalPolynomials.HermiteH(1, x) == 2.0 * x);
                Assert.IsTrue(OrthogonalPolynomials.HermiteHe(0, x) == 1.0);
                Assert.IsTrue(OrthogonalPolynomials.HermiteHe(1, x) == x);
            }
        }

        [TestMethod]
        public void HermiteSpecialCases () {
            // Values at 0; A&S 22.4.8
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 8)) {
                if (n % 2 == 0) {
                    int m = n / 2;
                    int s = (m % 2 == 0) ? 1 : -1;
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.HermiteH(n, 0.0), s * AdvancedIntegerMath.Factorial(n) / AdvancedIntegerMath.Factorial(m)));
                } else {
                    Assert.IsTrue(OrthogonalPolynomials.HermiteH(n, 0.0) == 0.0);
                }
            }
        }

        [TestMethod]
        public void HermiteReflection () {
            // A&S 22.4.8
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 10)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0, 1000.0, 10)) {
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
            // A&S 22.7.13
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 10)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0, 1000.0, 10)) {
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        OrthogonalPolynomials.HermiteH(n + 1, x), 2.0 * n * OrthogonalPolynomials.HermiteH(n - 1, x),
                        2.0 * x * OrthogonalPolynomials.HermiteH(n, x)
                    ));
                }
            }
        }

        [TestMethod]
        public void HermiteTuranInequality () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E3, 8)) {
                    Assert.IsTrue(
                        MoreMath.Sqr(OrthogonalPolynomials.HermiteH(n, x)) >= OrthogonalPolynomials.HermiteH(n + 1, x) * OrthogonalPolynomials.HermiteH(n - 1, x)
                    );
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
            // Don't let order get too high, or integral will become too oscilatory to handle to high precision numerically.
            int[] orders = TestUtilities.GenerateIntegerValues(1, 12, 4);
            foreach (int n in orders) {
                foreach (int m in orders) {
                    double I = FunctionMath.Integrate(
                        (double x) => Math.Exp(-x * x / 2.0) * OrthogonalPolynomials.HermiteHe(n, x) * OrthogonalPolynomials.HermiteHe(m, x),
                        Double.NegativeInfinity, Double.PositiveInfinity, new IntegrationSettings() { EvaluationBudget = 10000 }
                    );
                    if (n == m) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(I, Math.Sqrt(2.0 * Math.PI) * AdvancedIntegerMath.Factorial(n)));
                    } else {
                        Assert.IsTrue(Math.Abs(I) <= TestUtilities.TargetPrecision);
                    }
                }
            }
        }

        [TestMethod]
        public void HermiteStatisticianPhysicistAgreement () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 8)) {
                    TestUtilities.IsNearlyEqual(
                        OrthogonalPolynomials.HermiteH(n, x),
                        Math.Pow(2.0, n / 2.0) * OrthogonalPolynomials.HermiteHe(n, Math.Sqrt(2.0) * x)
                    );
                }
            }

        }

        [TestMethod]
        public void HermiteInvalidOrder () {
            Assert.ThrowsException<ArgumentOutOfRangeException>(() => OrthogonalPolynomials.HermiteH(-1, 1.0));
            Assert.ThrowsException<ArgumentOutOfRangeException>(() => OrthogonalPolynomials.HermiteHe(-1, 1.0));
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
        public void LaguerreLowOrders () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 8)) {
                Assert.IsTrue(OrthogonalPolynomials.LaguerreL(0, x) == 1.0);
                Assert.IsTrue(OrthogonalPolynomials.LaguerreL(1, x) == 1.0 - x);
            }
        }

        [TestMethod]
        public void LaguerreSpecialCases () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                Assert.IsTrue(OrthogonalPolynomials.LaguerreL(n, 0.0) == 1.0);
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
            // A&S 22.14.12
            foreach (int n in TestUtilities.GenerateRealValues(1.0, 100, 5)) {
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
        public void AssociatedLaguerreLowOrders () {
            foreach (double a in TestUtilities.GenerateUniformRealValues(-1.0, 10.0, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 4)) {
                    Assert.IsTrue(OrthogonalPolynomials.LaguerreL(0, a, x) == 1.0);
                    Assert.IsTrue(OrthogonalPolynomials.LaguerreL(1, a, x) == 1.0 + a - x);
                }
            }
        }

        [TestMethod]
        public void AssociatedLaguerreSpecialCases () {
            // A&S 22.4.7 gives value at 0
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (double a in TestUtilities.GenerateUniformRealValues(-1.0, 10.0, 4)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        OrthogonalPolynomials.LaguerreL(n, a, 0.0),
                        AdvancedMath.Pochhammer(a + 1.0, n) / AdvancedIntegerMath.Factorial(n)
                    ));
                    // This is \binomial{n + \alpha}{n}. Maybe provide a double-valued version?
                }
            }
        }

        [TestMethod]
        public void AssociatedLaguerreTriangleRecurrence () {
            // Connects (n, a - 1) (n, a) (n - 1, a)
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
        public void AssociatedLaguerreAlphaRecurrence () {
            // \sum_{k = 0}^{n} L^{(\alpha)}_k = L^{(\alpha+1)}_n
            foreach (int n in TestUtilities.GenerateRealValues(1, 100, 5)) {
                foreach (double a in TestUtilities.GenerateRealValues(0.1, 100.0, 5)) {
                    foreach (double x in TestUtilities.GenerateRealValues(0.01, 1000.0, 5)) {

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
        public void AssociatedLaguerreInequality () {
            // A&S 22.14.13
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (double a in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 4)) {
                    foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 4)) {
                        Assert.IsTrue(Math.Abs(OrthogonalPolynomials.LaguerreL(n, a, x)) <= AdvancedMath.Pochhammer(a + 1.0, n) / AdvancedIntegerMath.Factorial(n) * Math.Exp(x / 2.0));
                    }
                }
            }
        }

        [TestMethod]
        public void AssociatedLaguerreAgreement () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-3, 1.0E3, 4)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.LaguerreL(n, x), OrthogonalPolynomials.LaguerreL(n, 0.0, x)));
                }
            }
        }

        [TestMethod]
        public void AssociatedLaguerreInvalidArguments () {
            Assert.ThrowsException<ArgumentOutOfRangeException>(() => OrthogonalPolynomials.LaguerreL(-1, 0.0, 0.0));
            Assert.ThrowsException<ArgumentOutOfRangeException>(() => OrthogonalPolynomials.LaguerreL(0, -1.1, 0.0));
            Assert.ThrowsException<ArgumentOutOfRangeException>(() => OrthogonalPolynomials.LaguerreL(0, 0.0, -0.1));
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
                        IntegrationSettings e = new IntegrationSettings();
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
        public void ChebyshevTLowOrders () {
            foreach (double x in TestUtilities.GenerateUniformRealValues(-1.0, +1.0, 4)) {
                Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(0, x) == 1.0);
                Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(1, x) == x);
            }
        }

        [TestMethod]
        public void ChebyshevTStandardization () {
            // A&S 22.2.4
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 8)) {
                Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(n, 1.0) == 1.0);
            }
        }

        [TestMethod]
        public void ChebyshevTSpecialCases () {
            // Values at 0 and -1
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 8)) {
                if (n % 2 == 0) {
                    if (n % 4 == 0) {
                        Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(n, 0.0) == 1.0);
                    } else {
                        Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(n, 0.0) == -1.0);
                    }
                    Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(n, -1.0) == 1.0);
                } else {
                    Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(n, 0.0) == 0.0);
                    Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(n, -1.0) == -1.0);
                }
            }
        }

        [TestMethod]
        public void ChebyshevTInequality () {
            // A&S 22.14.1
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 10)) {
                foreach (double x in TestUtilities.GenerateUniformRealValues(-1.0, 1.0, 10)) {
                    Assert.IsTrue(Math.Abs(OrthogonalPolynomials.ChebyshevT(n, x)) <= 1.0);
                }
            }
        }

        [TestMethod]
        public void ChebyshevTTuranInequality () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 10)) {
                foreach (double x in TestUtilities.GenerateUniformRealValues(-1.0, 1.0, 10)) {
                    Assert.IsTrue(
                        MoreMath.Sqr(OrthogonalPolynomials.ChebyshevT(n, x)) >= OrthogonalPolynomials.ChebyshevT(n + 1, x) * OrthogonalPolynomials.ChebyshevT(n - 1, x)
                    );
                }
            }
        }

        [TestMethod]
        public void ChebyshevTRecurrence () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 10)) {
                foreach (double x in bArguments) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ChebyshevT(n + 1, x), 2.0 * x * OrthogonalPolynomials.ChebyshevT(n, x) - OrthogonalPolynomials.ChebyshevT(n - 1, x)));
                }
            }
        }

        [TestMethod]
        public void ChebyshevTReflection () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 10)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 4)) {
                    if (n % 2 == 0) {
                        Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(n, -x) == OrthogonalPolynomials.ChebyshevT(n, x));
                    } else {
                        Assert.IsTrue(OrthogonalPolynomials.ChebyshevT(n, -x) == -OrthogonalPolynomials.ChebyshevT(n, x));
                    }
                }
            }
        }

        [TestMethod]
        public void ChebyshevTOrderDoubling () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1 , 100, 6)) {
                foreach (double x in bArguments) {
                    if (n % 2 == 0) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ChebyshevT(2 * n, x), OrthogonalPolynomials.ChebyshevT(n, 1.0 - 2 * x * x)));
                    } else {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ChebyshevT(2 * n, x), -OrthogonalPolynomials.ChebyshevT(n, 1.0 - 2 * x * x)));
                    }
                }
            }
        }

        [TestMethod]
        public void ChebyshevTProduct () {
            // 2 T_m(x) T_n(x) = T_{m + n}(x) + T_{|m - n|}(x)
            Random rng = new Random(1);
            foreach (int m in TestUtilities.GenerateIntegerValues(1, 100, 6, rng)) {
                foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5, rng)) {
                    foreach (double x in TestUtilities.GenerateUniformRealValues(-1.0, 1.0, 4)) {
                        Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                            OrthogonalPolynomials.ChebyshevT(m + n, x), OrthogonalPolynomials.ChebyshevT(Math.Abs(m - n), x),
                            2.0 * OrthogonalPolynomials.ChebyshevT(m, x) * OrthogonalPolynomials.ChebyshevT(n, x)
                        ));
                    }
                }
            }
        }

        [TestMethod]
        public void ChebyshevTPermutability () {
            // An unusual feature of Chebyshev polynomials is that T_n(T_m(x)) = T_m(T_n(x)).
            // This is a consequence of T_n(T_m(x)) = T_{n m}(x).
            Random rng = new Random(1);
            foreach (double x in TestUtilities.GenerateUniformRealValues(-1.0, 1.0, 4)) {
                foreach (int n in TestUtilities.GenerateIntegerValues(1, 16, 3, rng)) {
                    foreach (int m in TestUtilities.GenerateIntegerValues(1, 16, 3, rng)) {
                        double y1 = OrthogonalPolynomials.ChebyshevT(n, OrthogonalPolynomials.ChebyshevT(m, x));
                        double y2 = OrthogonalPolynomials.ChebyshevT(m, OrthogonalPolynomials.ChebyshevT(n, x));
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(y1, y2, 1.0E-8));
                    }
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
        public void ChebyshevTCosine () {
            // A&S 22.3.15: T_n(cos t) = cos(n t)
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double t in TestUtilities.GenerateUniformRealValues(-Math.PI, Math.PI, 5)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ChebyshevT(n, Math.Cos(t)), Math.Cos(n * t)));
                }
            }
        }


        [TestMethod]
        public void ChebyshevTInvalidArguments () {
            Assert.ThrowsException<ArgumentOutOfRangeException>(() => OrthogonalPolynomials.ChebyshevT(-1, 0.5));
            Assert.ThrowsException<ArgumentOutOfRangeException>(() => OrthogonalPolynomials.ChebyshevT(2, -1.1));
        }

        // Orthonormality test for chebyshev fails because endpoint singularity of weight causes
        // the numerical integral not to converge

        [TestMethod]
        public void ChebyshevULowOrders () {
            foreach (double x in TestUtilities.GenerateUniformRealValues(-1.0, 1.0, 8)) {
                Assert.IsTrue(OrthogonalPolynomials.ChebyshevU(0, x) == 1.0);
                Assert.IsTrue(OrthogonalPolynomials.ChebyshevU(1, x) == 2.0 * x);
            }
        }

        [TestMethod]
        public void ChebyshevUStandardization() {
            // A&S 2.2.5
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 8)) {
                Assert.IsTrue(OrthogonalPolynomials.ChebyshevU(n, 1.0) == n + 1);
            }
        }

        [TestMethod]
        public void ChebyshevUOrthonormality () {
            // A&S 22.2.5
            int[] orders = TestUtilities.GenerateIntegerValues(1, 100, 4);
            foreach (int m in orders) {
                foreach (int n in orders) {
                    double I = FunctionMath.Integrate((double x) => Math.Sqrt((1.0 + x) * (1.0 - x)) * OrthogonalPolynomials.ChebyshevU(m, x) * OrthogonalPolynomials.ChebyshevU(n, x), -1.0, +1.0);
                    if (m == n) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(I, Math.PI / 2.0));
                    } else {
                        Assert.IsTrue(Math.Abs(I) < TestUtilities.TargetPrecision);
                    }
                }
            }
        }

        [TestMethod]
        public void ChebyshevUReflection() {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 10)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 4)) {
                    if (n % 2 == 0) {
                        Assert.IsTrue(OrthogonalPolynomials.ChebyshevU(n, -x) == OrthogonalPolynomials.ChebyshevU(n, x));
                    } else {
                        Assert.IsTrue(OrthogonalPolynomials.ChebyshevU(n, -x) == -OrthogonalPolynomials.ChebyshevU(n, x));
                    }
                }
            }
        }

        [TestMethod]
        public void ChebyshevUTuranEquality() {
            // For Chebyshev U, the difference of the two sides of the Turan inequality is known to be 1.
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (double x in TestUtilities.GenerateUniformRealValues(-1.0, +1.0, 4)) {
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        OrthogonalPolynomials.ChebyshevU(n + 1, x) * OrthogonalPolynomials.ChebyshevU(n - 1, x), 1.0,
                        MoreMath.Sqr(OrthogonalPolynomials.ChebyshevU(n, x))
                    ));
                }
            }
        }

        [TestMethod]
        public void ChebyshevUChebyshevTRelationship () {
            // A&S 22.5.6
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (double x in TestUtilities.GenerateUniformRealValues(-1.0, +1.0, 4)) {
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        OrthogonalPolynomials.ChebyshevU(n, x), - x * OrthogonalPolynomials.ChebyshevU(n - 1, x),
                        OrthogonalPolynomials.ChebyshevT(n, x)
                    ));
                }
            }
        }

        [TestMethod]
        public void ChebyshevUInvalidArguments() {
            Assert.ThrowsException<ArgumentOutOfRangeException>(() => OrthogonalPolynomials.ChebyshevU(-1, 0.0));
            Assert.ThrowsException<ArgumentOutOfRangeException>(() => OrthogonalPolynomials.ChebyshevU(0, -1.1));
        }

        // Legendre

        [TestMethod]
        public void LegendreLowOrders () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-3, 1.0, 4)) {
                Assert.IsTrue(OrthogonalPolynomials.LegendreP(0, x) == 1.0);
                Assert.IsTrue(OrthogonalPolynomials.LegendreP(1, x) == x);
            }
        }

        [TestMethod]
        public void LegendreStandardization() {
            // A&S 22.2.10
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 8)) {
                Assert.IsTrue(OrthogonalPolynomials.LegendreP(n, 1.0) == 1.0);
            }
        }

        [TestMethod]
        public void LegendreSpecialCases () {
            // Values at 0, A&S 22.4.6
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 8)) {
                if (n % 2 == 0) {
                    int m = n / 2;
                    int s = (m % 2 == 0) ? 1 : -1;
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.LegendreP(n, 0.0), s * AdvancedIntegerMath.BinomialCoefficient(n, m) / MoreMath.Pow(2.0, n)));
                } else {
                    Assert.IsTrue(OrthogonalPolynomials.LegendreP(n, 0.0) == 0.0);
                }
            }
        }

        [TestMethod]
        public void LegendreInequality () {
            // A&S 22.14.7
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 8)) {
                foreach (double x in TestUtilities.GenerateUniformRealValues(-1.0, 1.0, 8)) {
                    Assert.IsTrue(Math.Abs(OrthogonalPolynomials.LegendreP(n, x)) <= 1.0);
                }
            }
        }

        [TestMethod]
        public void LegendreTuranInequality () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 8)) {
                foreach (double x in TestUtilities.GenerateUniformRealValues(-1.0, 1.0, 8)) {
                    Assert.IsTrue(
                        MoreMath.Sqr(OrthogonalPolynomials.LegendreP(n, x)) >= OrthogonalPolynomials.LegendreP(n + 1, x) * OrthogonalPolynomials.LegendreP(n - 1, x)
                    );
                }
            }
        }

        [TestMethod]
        public void LegendreRecurrence () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 8)) {
                foreach (double x in TestUtilities.GenerateUniformRealValues(-1.0, 1.0, 8)) {
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        (2 * n + 1) * x * OrthogonalPolynomials.LegendreP(n, x), -n * OrthogonalPolynomials.LegendreP(n - 1, x),
                        (n + 1) * OrthogonalPolynomials.LegendreP(n + 1, x)
                    ));
                }
            }
        }

        [TestMethod]
        public void LegendreReflection () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-3, 1.0, 4)) {
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
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 8)) {
                foreach (double x in TestUtilities.GenerateUniformRealValues(-1.0, 1.0, 8)) {
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
            int[] orders = TestUtilities.GenerateIntegerValues(1, 100, 4);
            foreach (int n in orders) {
                foreach (int m in orders) {
                    Func<double, double> f = (double x) => OrthogonalPolynomials.LegendreP(n, x) * OrthogonalPolynomials.LegendreP(m, x);
                    Interval r = Interval.FromEndpoints(-1.0, 1.0);
                    double I = FunctionMath.Integrate(f, r);
                    if (n == m) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(I, 2.0 / (2 * n + 1)));
                    } else {
                        Assert.IsTrue(Math.Abs(I) < TestUtilities.TargetPrecision);
                    }
                }
            }
        }

        [TestMethod]
        public void GegenbauerLowOrders() {
            foreach (double alpha in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 4))
                foreach (double x in TestUtilities.GenerateUniformRealValues(-1.0, 1.0, 4)) {
                    Assert.IsTrue(OrthogonalPolynomials.GegenbauerC(0, alpha, x) == 1.0);
                    Assert.IsTrue(OrthogonalPolynomials.GegenbauerC(1, alpha, x) == 2.0 * alpha * x);
                }
        }

        [TestMethod]
        public void GegenbauerOrthonormality () {
            int[] orders = TestUtilities.GenerateIntegerValues(1, 100, 4);
            foreach (int n in orders) {
                foreach (int m in orders) {
                    foreach (double alpha in TestUtilities.GenerateUniformRealValues(0.0, 2.0, 4)) {

                        IntegrationResult r = FunctionMath.Integrate(
                            x => OrthogonalPolynomials.GegenbauerC(n, alpha, x) * OrthogonalPolynomials.GegenbauerC(m, alpha, x) * Math.Pow((1.0 - x) * (1.0 + x), alpha - 0.5),
                            -1.0, 1.0
                        );

                        if (n == m) {
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(r.Value,
                                Math.PI * Math.Pow(2.0, 1.0 - 2.0 * alpha) / (n + alpha) * AdvancedMath.Gamma(n + 2.0 * alpha) / MoreMath.Sqr(AdvancedMath.Gamma(alpha)) / AdvancedIntegerMath.Factorial(n),
                                Math.Sqrt(TestUtilities.TargetPrecision)
                            ));
                        } else {
                            Assert.IsTrue(Math.Abs(r.Value) < Math.Sqrt(TestUtilities.TargetPrecision));
                        }

                    }
                }
            }
        } 

        [TestMethod]
        public void GegenbauerAgreement () {
            double alpha = 1.0 / Math.Sqrt(Double.MaxValue);
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (double x in TestUtilities.GenerateUniformRealValues(-1.0, 1.0, 8)) {

                    // A&S 22.5.36
                    // P_n(x) = C^{(1/2)})_n(x)
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.LegendreP(n, x), OrthogonalPolynomials.GegenbauerC(n, 1.0 / 2.0, x)));

                    // A&S 22.5.33
                    // T_n(x) = \lim_{\alpha \rightarrow 0} \frac{n}{2\alpha} C^{(\alpha)}_n(x)
                    // Since C^{(0)}_n(x) vanishes, this relationship must be defined via a limit.
                    // We take a value of alpha such that its square will always vanish.
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ChebyshevT(n, x), n / 2.0 * OrthogonalPolynomials.GegenbauerC(n, alpha, x) / alpha));

                    // Bizarrely, A&S 22.5.4 defines C^{(0)}_n(x) = \lim_{\alpha \rightarrow 0} C^{(\alpha)}_n(x) / \alpha, which means
                    // C^{(0)}_n(x) (read one way) does not equal C^{(0)}_n(x) (read the other way). It is only because of this re-definition
                    // that A&S 22.5.33 holds, which confused me for a while.

                    // A&S 22.5.34
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ChebyshevU(n, x), OrthogonalPolynomials.GegenbauerC(n, 1.0, x)));
                }
            }

        }

        [TestMethod]
        public void GegenbauerInvalidArguments () {
            Assert.ThrowsException<ArgumentOutOfRangeException>(() => OrthogonalPolynomials.GegenbauerC(-1, 0.0, 0.0));
            Assert.ThrowsException<ArgumentOutOfRangeException>(() => OrthogonalPolynomials.GegenbauerC(0, -1.0, 0.0));
            Assert.ThrowsException<ArgumentOutOfRangeException>(() => OrthogonalPolynomials.GegenbauerC(0, 0, -1.1));
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
        public void ZernikeLowOrders () {
            foreach (double x in TestUtilities.GenerateUniformRealValues(0.0, 1.0, 4)) {
                Assert.IsTrue(OrthogonalPolynomials.ZernikeR(0, 0, x) == 1.0);
                Assert.IsTrue(OrthogonalPolynomials.ZernikeR(1, 0, x) == 0.0);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ZernikeR(1, 1, x), x));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(OrthogonalPolynomials.ZernikeR(2, 0, x), 2.0 * x * x - 1.0));
                Assert.IsTrue(OrthogonalPolynomials.ZernikeR(2, 1, x) == 0.0);
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

            // \int_0^{1} \! R^{m}_{n}(\rho) J_{m}(x \rho) \rho \, d\rho = (-1)^{(n-m)/2} J_{n+1}(x) / x

            foreach (int n in TestUtilities.GenerateIntegerValues(2, 32, 4)) {
                foreach (int m in TestUtilities.GenerateUniformIntegerValues(0, n, 2)) {
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

                        // Near-zero integrals result from delicate cancelations, and thus
                        // cannot be expected to achieve the relative target precision;
                        // they can achieve the absolute target precisions, though.
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(I, J, new EvaluationSettings() { RelativePrecision = 0.0, AbsolutePrecision = TestUtilities.TargetPrecision }));

                    }
                }
            }

        }

        [TestMethod]
        public void ZernickeInvalidArguments () {
            Assert.ThrowsException<ArgumentOutOfRangeException>(() => OrthogonalPolynomials.ZernikeR(-1, 0, 0.0));
            Assert.ThrowsException<ArgumentOutOfRangeException>(() => OrthogonalPolynomials.ZernikeR(0, 0, 1.1));
        }

        [TestMethod]
        public void AssociatedLegendreLowOrders () {
            foreach (double x in TestUtilities.GenerateRealValues(0.01, 1.0, 10)) {
                Assert.IsTrue(OrthogonalPolynomials.LegendreP(0, 0, x) == 1.0);
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
