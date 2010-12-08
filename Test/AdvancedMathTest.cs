using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Test
{
    
    
    /// <summary>
    ///This is a test class for AdvancedMathTest and is intended
    ///to contain all AdvancedMathTest Unit Tests
    ///</summary>
    [TestClass()]
    public class AdvancedMathTest {


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


        // Bessel function tests

        private int[] orders = new int[] { 3, 30, 300 };
        private double[] arguments = new double[] { 0.11, 1.1, 11.1, 111.1 };

        private bool BesselYInRange (double n, double x) {

            if (x > n) {
                return (true);
            } else {
                if (AdvancedMath.LogGamma(n) - n * Math.Log(x / 2.0) > Math.Log(1.0e-16 * Double.MaxValue)) {
                    return (false);
                } else {
                    return (true);
                }
            }

        }

        [TestMethod]
        public void BesselNegativeOrderTest () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 15)) {
                    int s = 1 - 2*(n % 2);
                    Assert.IsTrue(AdvancedMath.BesselJ(-n, x) == s * AdvancedMath.BesselJ(n, x));
                    Assert.IsTrue(AdvancedMath.BesselY(-n, x) == s * AdvancedMath.BesselY(n, x));
                }
            }
        }

        [TestMethod]
        public void BesselRecurrenceTest () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1,100,5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E4,15)) {
                    Console.WriteLine("n = {0}, x={1}", n, x);
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(AdvancedMath.BesselJ(n - 1, x), AdvancedMath.BesselJ(n + 1, x), 2 * n / x * AdvancedMath.BesselJ(n, x)));
                    if (!BesselYInRange(n, x)) continue;
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(AdvancedMath.BesselY(n - 1, x), AdvancedMath.BesselY(n + 1, x), 2 * n / x * AdvancedMath.BesselY(n, x)), String.Format("n={0}, x={1} Y(n-1,x)={2}, Y(n,x)={3}, Y(n+1,x)={4}", n, x, AdvancedMath.BesselY(n - 1, x), AdvancedMath.BesselY(n, x), AdvancedMath.BesselY(n + 1, x)));
                }
            }
        }

        [TestMethod]
        public void BesselWronskianTest () {
            foreach (int n in orders) {
                foreach (double x in arguments) {
                    if ((x > n) || (Math.Exp(AdvancedMath.LogGamma(n) - n * Math.Log(x)) < Math.Sqrt(Double.MaxValue))) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(2 / Math.PI / x, AdvancedMath.BesselJ(n + 1, x) * AdvancedMath.BesselY(n, x) - AdvancedMath.BesselJ(n, x) * AdvancedMath.BesselY(n + 1, x)), String.Format("n={0}, x={1}, J(n,x)={2}, J(n+1,x) = {3}, Y(n,x)={4}, Y(n+1,x)={5}", n, x, AdvancedMath.BesselJ(n, x), AdvancedMath.BesselJ(n + 1, x), AdvancedMath.BesselY(n, x), AdvancedMath.BesselY(n + 1, x)));
                    }
                }
            }
        }

        [TestMethod]
        public void BesselJ0IntegralTest () {
            Interval r = Interval.FromEndpoints(0.0, Math.PI);
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 5)) {
                Function<double, double> f = delegate(double t) {
                    return (Math.Cos(x * Math.Sin(t)));
                };
                double J0 = FunctionMath.Integrate(f, r) / Math.PI;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.BesselJ(0, x), J0));
            }
        }

        [TestMethod]
        public void BesselY0IntegralTest () {
            Interval r = Interval.FromEndpoints(0.0, Math.PI / 2.0);
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 5)) {
                Function<double, double> f = delegate(double t) {
                    double s = Math.Sin(t);
                    return (Math.Cos(x * Math.Cos(t)) * (AdvancedMath.EulerGamma + Math.Log(2.0 * x * s * s)));
                };
                double Y0 = 4.0 / (Math.PI * Math.PI) * FunctionMath.Integrate(f, r);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.BesselY(0, x), Y0));
            }
        }

        [TestMethod]
        public void BesselJIntegralTest () {
            Interval r = Interval.FromEndpoints(0.0, Math.PI);
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E1, 5)) {
                foreach (int n in TestUtilities.GenerateIntegerValues(1, 10, 3)) {
                    Function<double, double> f = delegate(double t) {
                        return (Math.Cos(x * Math.Sin(t) - n * t));
                    };
                    double J = FunctionMath.Integrate(f, r) / Math.PI;
                    Console.WriteLine("{0} {1} {2} {3}", n, x, AdvancedMath.BesselJ(n, x), J);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.BesselJ(n, x), J));
                }
            }
        }

        [TestMethod]
        public void BesselKapteynIntegralTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 5)) {
                Function<double, double> f = delegate(double t) {
                    return (Math.Cos(x - t) * AdvancedMath.BesselJ(0, t));
                };
                Interval r = Interval.FromEndpoints(0.0, x);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(x * AdvancedMath.BesselJ(0, x), FunctionMath.Integrate(f, r)));
            }
        }

        [TestMethod]
        public void BesselLipshitzIntegralTest () {
            Function<double, double> f = delegate(double t) {
                return (Math.Exp(-t) * AdvancedMath.BesselJ(0, t));
            };
            Interval r = Interval.FromEndpoints(0.0, Double.PositiveInfinity);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(FunctionMath.Integrate(f, r), 1.0 / Math.Sqrt(2.0)));
        }

        [TestMethod]
        public void BesselWeberIntegralTest () {
            Function<double, double> f = delegate(double t) {
                return (Math.Exp(-t*t) * AdvancedMath.BesselJ(0, t) * t);
            };
            Interval r = Interval.FromEndpoints(0.0, Double.PositiveInfinity);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(FunctionMath.Integrate(f, r), Math.Exp(-1.0/4.0) / 2.0));
        }

        [TestMethod]
        public void BesselAgreementTest () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1,100,5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E4,15)) {
                    Console.WriteLine("n={0},x={1}", n, x);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.BesselJ(n, x), AdvancedMath.BesselJ((double) n, x)));
                    if (!BesselYInRange(n, x)) continue; // don't try to evaluate Y if it's too big
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.BesselY(n, x), AdvancedMath.BesselY((double) n, x)), String.Format("n={0}, x={1} YI={2} YR={3}", n, x, AdvancedMath.BesselY(n, x), AdvancedMath.BesselY((double) n, x)));
                }
            }
        }

        [TestMethod]
        public void RealBesselWronskianTest () {
            foreach (double nu in orders) {
                foreach (double x in arguments) {
                    if (BesselYInRange(nu, x)) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(2 / Math.PI / x, AdvancedMath.BesselJ(nu + 1.0, x) * AdvancedMath.BesselY(nu, x) - AdvancedMath.BesselJ(nu, x) * AdvancedMath.BesselY(nu + 1.0, x)), String.Format("n={0}, x={1}, J(n,x)={2}, J(n+1,x) = {3}, Y(n,x)={4}, Y(n+1,x)={5}", nu, x, AdvancedMath.BesselJ(nu, x), AdvancedMath.BesselJ(nu + 1, x), AdvancedMath.BesselY(nu, x), AdvancedMath.BesselY(nu + 1, x)));
                    }
                }
            }
        }

        [TestMethod]
        public void RealBesselInequalityTest () {
            foreach (double nu in TestUtilities.GenerateRealValues(0.5, 50.0, 3)) {
                foreach (double x in TestUtilities.GenerateRealValues(nu, 100.0, 3)) {
                    if (!BesselYInRange(nu, x)) continue;
                    double J = AdvancedMath.BesselJ(nu, x);
                    double Y = AdvancedMath.BesselY(nu, x);
                    double S = J * J + Y * Y;
                    Assert.IsTrue(2.0 / Math.PI / x <= S);
                    Assert.IsTrue(S <= 2.0 / Math.PI / Math.Sqrt(x * x - nu * nu));
                }
            }
        }

        [TestMethod]
        public void RealBesselFresnelTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2,1.0E2,5)) {
                Console.WriteLine("x={0}", x);
                double csum = 0.0;
                double ssum = 0.0;
                int k = 0;
                while (true) {
                    double csum_old = csum;
                    double ssum_old = ssum;
                    csum += AdvancedMath.BesselJ((4 * k + 1) / 2.0, x);
                    ssum += AdvancedMath.BesselJ((4 * k + 3) / 2.0, x);
                    if ((csum == csum_old) && (ssum == ssum_old)) break;
                    k += 1;
                    if (k > 100) throw new NonconvergenceException();
                }
                double C = AdvancedMath.FresnelC(Math.Sqrt(2.0 * x / Math.PI));
                double S = AdvancedMath.FresnelS(Math.Sqrt(2.0 * x / Math.PI));
                Console.WriteLine("k={0}, csum={1}, C={2}", k, csum, C);
                Console.WriteLine("k={0}, ssum={1}, S={2}", k, ssum, S);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(csum, C));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ssum, S));
            }
        }
        

        [TestMethod]
        public void BesselTowerTest () {

            // x shouldn't be too big or we will need too many terms before the series converges
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2,1.0E2,10) ) {

                // S1 = J0 + 2 J2 + 2 J4 + 2 J8 + ...
                // S2 = J0 - 2 J2 + 2 J4 - 2 J8 + ...
                // S3 = 2 J1 - 2 J3 + 2 J5 - 2 J7 + ...

                double J0 = AdvancedMath.BesselJ(0,x);
                double s0 = J0 * J0;
                double s1 = J0;
                double s2 = J0;
                double s3 = 0.0;

                int n = 1;
                while (true) {

                    double s0_old = s0;
                    double s1_old = s1;
                    double s2_old = s2;
                    double s3_old = s3;

                    double J1 = AdvancedMath.BesselJ(2*n - 1, x);
                    double J2 = AdvancedMath.BesselJ(2*n, x);
                    int s = 2 * (n % 2) - 1;

                    s0 += 2.0 * (J1 * J1 + J2 * J2);
                    s1 += 2.0 * J2;
                    s2 += -s * 2.0 * J2;
                    s3 += s * 2.0 * J1;

                    if ((s0 == s0_old) && (s1 == s1_old) && (s2 == s2_old) && (s3 == s3_old)) break;

                    n++;
                    if (n > 100) throw new NonconvergenceException();

                }
                Console.WriteLine("x={0}, n={1}", x, n);

                Console.WriteLine("{0} {1}", s0, 1.0);
                Console.WriteLine("{0} {1}", s1, 1.0);
                Console.WriteLine("{0} {1}", s2, Math.Cos(x));
                Console.WriteLine("{0} {1}", s3, Math.Sin(x));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(s0, 1.0));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s1, 1.0));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s2, Math.Cos(x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s3, Math.Sin(x)));

            }

        }

        [TestMethod]
        public void RealBesselJIntegralTest () {
            foreach (double nu in TestUtilities.GenerateRealValues(0.1, 10.0, 3)) {
                foreach (double x in TestUtilities.GenerateRealValues(0.01, 100.0, 3)) {
                    Function<double, double> f = delegate(double t) {
                        return (Math.Cos(x * Math.Cos(t)) * Math.Pow(Math.Sin(t), 2.0 * nu));
                    };
                    Interval r = Interval.FromEndpoints(0.0, Math.PI);
                    double I = FunctionMath.Integrate(f, r);
                    I = I * Math.Pow(x / 2.0, nu) / AdvancedMath.Gamma(nu + 0.5) / AdvancedMath.Gamma(0.5);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.BesselJ(nu, x), I));
                }
            }
        }

        [TestMethod]
        public void SphericalBesselSpecialCaseTest () {
            Assert.IsTrue(AdvancedMath.SphericalBesselJ(0, 0.0) == 1.0);
            Assert.IsTrue(AdvancedMath.SphericalBesselJ(1, 0.0) == 0.0);
            //Assert.IsTrue(AdvancedMath.SphericalBesselY(0, 0.0) == Double.NegativeInfinity);
            Assert.IsTrue(AdvancedMath.SphericalBesselY(1, 0.0) == Double.NegativeInfinity);
        }

        [TestMethod]
        public void SphericalBesselRecurrenceTest () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1,100,5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E4,15)) {
                    Console.WriteLine("n={0} x={1}", n, x);
                    double jm1 = AdvancedMath.SphericalBesselJ(n - 1, x);
                    double jp1 = AdvancedMath.SphericalBesselJ(n + 1, x);
                    double j = AdvancedMath.SphericalBesselJ(n, x);
                    Console.WriteLine("  {0:g16} + {1:g16} = {2:g16} ?= {3:g16}", jm1, jp1, jm1 + jp1, (2 * n + 1) * j / x);
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(jm1, jp1, (2*n+1)*j/x));
                    if (BesselYInRange(n, x)) {
                        double ym1 = AdvancedMath.SphericalBesselY(n - 1, x);
                        double yp1 = AdvancedMath.SphericalBesselY(n + 1, x);
                        double y = AdvancedMath.SphericalBesselY(n, x);
                        Console.WriteLine("  {0} + {1} = {2} ?= {3}", ym1, yp1, ym1 + yp1, (2 * n + 1) * y / x);
                        Assert.IsTrue(TestUtilities.IsSumNearlyEqual(AdvancedMath.SphericalBesselY(n - 1, x), AdvancedMath.SphericalBesselY(n + 1, x), (2 * n + 1) / x * AdvancedMath.SphericalBesselY(n, x)), "Y");
                    }

                }
            }
        }

        [TestMethod]
        public void SphericalBesselNegativeOrderTest () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 15)) {
                    int s = 2 * (n % 2) - 1;
                    Assert.IsTrue(AdvancedMath.SphericalBesselY(n, x) == s * AdvancedMath.SphericalBesselJ(-n - 1, x));
                    Assert.IsTrue(AdvancedMath.SphericalBesselY(-n, x) == s * AdvancedMath.SphericalBesselJ(n - 1, x));
                }
            }
        }

        [TestMethod]
        public void SphericalBesselWronskianTest () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 15)) {
                    if (BesselYInRange(n, x)) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.SphericalBesselJ(n, x) * AdvancedMath.SphericalBesselY(n - 1, x) - AdvancedMath.SphericalBesselJ(n - 1, x) * AdvancedMath.SphericalBesselY(n, x), 1.0 / (x * x)));
                    }
                }
            }
        }

        [TestMethod]
        public void SphericalBesselAgreementTest () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 15)) {
                    Console.WriteLine("n={0} x={1}", n, x);
                    Console.WriteLine("Sp = {0}", AdvancedMath.SphericalBesselJ(n, x));
                    Console.WriteLine("Re = {0}", Math.Sqrt(Math.PI / 2.0 / x) * AdvancedMath.BesselJ(n + 0.5, x));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.SphericalBesselJ(n, x), Math.Sqrt(Math.PI / 2.0 / x) * AdvancedMath.BesselJ(n + 0.5, x)));
                    if (BesselYInRange(n, x)) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.SphericalBesselY(n, x), Math.Sqrt(Math.PI / 2.0 / x) * AdvancedMath.BesselY(n + 0.5, x)));
                    }
                }
            }
        }

        [TestMethod]
        public void SphericalBesselTowerTest () {

            // x shouldn't be too big or we will need too many terms before the series converges
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 5)) {

                double s1 = 0.0;
                double s2 = 0.0;
                double s3 = 0.0;

                int n = 0;
                while (true) {

                    double j = AdvancedMath.SphericalBesselJ(n, x);
                    double ds1 = j * j;
                    double ds2 = (2 * n + 1) * ds1;
                    double ds3 = (1 - 2 * (n%2)) * ds2;

                    s1 += ds1;
                    s2 += ds2;
                    s3 += ds3;

                    if ((s1 + ds1 == s1) && (s2 + ds2 == s2) && (s3 + ds3 == s3)) break;

                    n++;
                    if (n > 100) throw new NonconvergenceException();

                }
                Console.WriteLine("x={0} n={0}", x, n);

                double x2 = 2 * x;
                Console.WriteLine("{0} {1}", s1, AdvancedMath.IntegralSi(x2)/x2);
                Console.WriteLine("{0} {1}", s2, 1.0);
                Console.WriteLine("{0} {1}", s3, Math.Sin(x2)/x2);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(s1, AdvancedMath.IntegralSi(x2)/x2));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s2, 1.0));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s3, Math.Sin(x2)/x2));

            }

        }

        // Gamma functions

        [TestMethod]
        public void GammaRecurrsionTest () {
            // limit x to avoid overflow
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E2,20)) {
                Console.WriteLine("x={0}, G={1}", x, AdvancedMath.Gamma(x));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(x + 1.0), x * AdvancedMath.Gamma(x)));
            }
        }

        [TestMethod]
        public void GammaSpecialValuesTest () {
            // it would be nice to be able to make some of these comparisons exact
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(0.5), Math.Sqrt(Math.PI)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(1.0), 1.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(1.5), Math.Sqrt(Math.PI) / 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(2.0), 1.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(3.0), 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(4.0), 6.0));
        }

        [TestMethod]
        public void GammaReflectionTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E2,20)) {
                // don't try a large near-integer, because trig functions of large numbers aren't accurate enough 
                if ((x > 20) && (Math.Abs(x - Math.Round(x)) < 0.05)) continue;
                double GP = AdvancedMath.Gamma(x);
                double GN = AdvancedMath.Gamma(-x);
                Console.WriteLine("x={0}, G(x)={1}, G(-x)={2}", x, GP, GN);
                Console.WriteLine("{0} ?= {1}", -x * GN * GP, Math.PI / Math.Sin(Math.PI * x));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(-x * GN * GP, Math.PI / Math.Sin(Math.PI * x)));
            }
        }

        [TestMethod]
        public void GammaInequalityTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0,1.0E2,10)) {
                if (x >= 2.0) {
                    double lower = Math.Pow(x / Math.E, x - 1.0);
                    double upper = Math.Pow(x / 2.0, x - 1.0);
                    double value = AdvancedMath.Gamma(x);
                    Assert.IsTrue((lower <= value) && (value <= upper));
                }
            }
        }

        [TestMethod]
        public void GammaIntegralTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0, 1.0E2, 3)) {
                Function<double,double> f = delegate (double t) {
                    return( Math.Pow(t, x - 1.0) * Math.Exp(-t) );
                };
                Interval r = Interval.FromEndpoints(0.0, Double.PositiveInfinity);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(x), FunctionMath.Integrate(f, r)));
            }
        }

        [TestMethod]
        public void LogGammaDuplicationTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E4,20)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.LogGamma(2.0 * x), AdvancedMath.LogGamma(x) + AdvancedMath.LogGamma(x+0.5) + (2.0*x-0.5) * Math.Log(2.0) - 0.5 * Math.Log(2.0*Math.PI)));
            }
        }

        [TestMethod]
        public void LogGammaTriplicationTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 20)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.LogGamma(3.0 * x), AdvancedMath.LogGamma(x) + AdvancedMath.LogGamma(x + 1.0 / 3.0) + AdvancedMath.LogGamma(x + 2.0 / 3.0) + (3.0 * x - 0.5) * Math.Log(3.0) - Math.Log(2.0 * Math.PI)));
            }
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void LogGammaNegativeArgumentTest () {
            AdvancedMath.LogGamma(-0.5);
        }

        [TestMethod]
        public void GammaRatioInequalityTest () {

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 15)) {

                double LB = Math.Log(Math.Sqrt(x + 0.25));
                double UB = Math.Log((x + 0.5) / Math.Sqrt(x + 0.75));
                double R = AdvancedMath.LogGamma(x + 1.0) - AdvancedMath.LogGamma(x + 0.5);

                Assert.IsTrue(LB <= R);
                Assert.IsTrue(R <= UB);

            }

        }

        [TestMethod]
        public void PsiSpecialCaseTest () {
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(0.5), -AdvancedMath.EulerGamma - 2.0 * Math.Log(2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(1.0), -AdvancedMath.EulerGamma));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(2.0), -AdvancedMath.EulerGamma + 1.0));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(1.0 / 3.0), -AdvancedMath.EulerGamma - 3.0 * Math.Log(3.0) / 2.0 - Math.PI / 2.0 / Math.Sqrt(3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(1.0 / 4.0), -AdvancedMath.EulerGamma - 3.0 * Math.Log(2.0) - Math.PI / 2.0));

        }

        [TestMethod]
        public void PsiRecurrenceTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E4,20)) {
                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(AdvancedMath.Psi(x), 1.0 / x, AdvancedMath.Psi(x + 1.0)));
            }
        }

        [TestMethod]
        public void PsiReflectionTest () {
            // don't get x get too big or inaccuracy of System.Math trig functions for large x will cause failure
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E2,15)) {
                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(AdvancedMath.Psi(x), Math.PI / Math.Tan(Math.PI * x), AdvancedMath.Psi(1.0 - x)));
            }
        }

        [TestMethod]
        public void PsiDuplicationTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 20)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(2 * x), AdvancedMath.Psi(x) / 2.0 + AdvancedMath.Psi(x + 0.5) / 2.0 + Math.Log(2.0)));
            }
        }

        // Incomplete gamma

        [TestMethod]
        public void RegularizedIncompleteGammaRecurrence () {
            foreach (double a in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 5)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.LeftRegularizedGamma(a + 1.0, x) + Math.Pow(x, a) * Math.Exp(-x) / AdvancedMath.Gamma(a + 1.0),
                        AdvancedMath.LeftRegularizedGamma(a, x)
                    ));
                }
            }
        }

        // add a=0 test when we allow a=0

        [TestMethod]
        public void RegularizedIncompleteGammaExponential () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.RightRegularizedGamma(1.0, x), Math.Exp(-x)
                ));
            }
        }

        [TestMethod]
        public void RegularizedIncompleteGammaUnitarity () {
            foreach (double a in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 10)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 10)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.LeftRegularizedGamma(a, x) + AdvancedMath.RightRegularizedGamma(a, x),
                        1.0
                    ));
                }
            }
        }

        [TestMethod]
        public void IncompleteGammaInequality () {
            foreach (double a in TestUtilities.GenerateRealValues(1.0E-4, 1, 10)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 10.0, 6)) {
                    Assert.IsTrue(
                        AdvancedMath.Gamma(a, x) <= Math.Exp(-x) * Math.Pow(x, a - 1.0)
                    );
                }
            }
        }

        [TestMethod]
        public void IncompleteGammaIntegerInequality () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 10000, 10)) {
                Assert.IsTrue(AdvancedMath.RightRegularizedGamma(n, n) <= 0.5);
                Assert.IsTrue(0.5 <= AdvancedMath.RightRegularizedGamma(n, n - 1));
            }
        }

        [TestMethod]
        public void IncompleteGammaErfcTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 20)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(Math.Sqrt(Math.PI) * AdvancedMath.Erfc(x), AdvancedMath.Gamma(0.5, x * x)));
            }
        }

        [TestMethod]
        public void IncompleteGammaExpTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 20)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(1.0, x), Math.Exp(-x)));
            }
        }

        [TestMethod]
        public void BetaIntegralTest () {
            // a and b values less than one correspond to singular integrals
            foreach (double a in TestUtilities.GenerateRealValues(1.0, 10.0, 3)) {
                foreach (double b in TestUtilities.GenerateRealValues(1.0, 10.0, 3)) {
                    double B = AdvancedMath.Beta(a, b);
                    Function<double, double> f = delegate(double t) {
                        return (Math.Pow(t, a - 1.0) * Math.Pow(1.0 - t, b - 1.0));
                    };
                    Interval r = Interval.FromEndpoints(0.0, 1.0);
                    double I = FunctionMath.Integrate(f, r);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(B, I));
                }
            }
        }

        [TestMethod]
        public void BetaTransformationTest () {
            foreach (double a in TestUtilities.GenerateRealValues(0.1, 10.0, 3)) {
                foreach (double b in TestUtilities.GenerateUniformRealValues(0.0, 1.0, 3)) {
                    double B1 = AdvancedMath.Beta(a, b);
                    double B2 = AdvancedMath.Beta(a + b, 1.0 - b);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(B1 * B2, Math.PI / a / Math.Sin(Math.PI * b)));
                }
            }
        }

        // Error function

        [TestMethod]
        public void ErfNegativeArgumentsTest () {
            foreach (double x in arguments) {
                Assert.IsTrue(AdvancedMath.Erf(-x) == -AdvancedMath.Erf(x));
            }
        }

        [TestMethod]
        public void ErfSpecialCasesTest () {
            Assert.IsTrue(AdvancedMath.Erf(0.0) == 0.0);
            Assert.IsTrue(AdvancedMath.Erfc(0.0) == 1.0);
        }

        [TestMethod]
        public void InverseErfTest () {
            //double[] probabilities = new double[] { 0.00005, 0.05, 0.5, 0.95, 0.99995 };
            foreach (double P in TestUtilities.GenerateRealValues(1.0E-5,1,10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Erf(AdvancedMath.InverseErf(P)), P));
            }
        }

        [TestMethod]
        public void InverseErfIntegralTest () {
            Function<double, double> f = delegate(double t) {
                return (AdvancedMath.InverseErf(t));
            };
            Interval r = Interval.FromEndpoints(0.0, 1.0);
            double I = FunctionMath.Integrate(f, r);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(I, 1.0 / Math.Sqrt(Math.PI)));
        }

        [TestMethod]
        public void ErfIntegralTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2,10.0,3)) {
                Function<double, double> f = delegate(double t) {
                    return (Math.Exp(-t * t));
                };
                Interval r = Interval.FromEndpoints(0.0, x);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Erf(x), 2.0 / Math.Sqrt(Math.PI) * FunctionMath.Integrate(f, r)));
            }
        }

        [TestMethod]
        public void ErfcComplementarityTest () {
            foreach (double x in arguments) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Erf(x) + AdvancedMath.Erfc(x), 1.0));
            }
        }

        [TestMethod]
        public void ErfcInequalityTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-5,1.0E5,20)) {
                double erfc = AdvancedMath.Erfc(x);
                double factor = 2.0 / Math.Sqrt(Math.PI) * Math.Exp(-x * x);
                double lower = factor / (x + Math.Sqrt(x * x + 2.0));
                double upper = factor / (x + Math.Sqrt(x * x + 4.0 / Math.PI));
                Assert.IsTrue((lower <= erfc) && (erfc <= upper));
            }

        }

        [TestMethod]
        public void ErfcIntegralTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 10.0, 3)) {
                Function<double, double> f = delegate(double t) {
                    return (Math.Exp(-t * t));
                };
                Interval r = Interval.FromEndpoints(x, Double.PositiveInfinity);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Erfc(x), 2.0 / Math.Sqrt(Math.PI) * FunctionMath.Integrate(f, r)));
            }
        }

        [TestMethod]
        public void IntegralESpecialCaseTest () {
            foreach (int n in orders) {
                Assert.IsTrue(AdvancedMath.IntegralE(n, 0.0) == 1.0 / (n - 1.0));
            }
            foreach (double x in arguments) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.IntegralE(0, x), Math.Exp(-x) / x));
            }
        }

        [TestMethod]
        public void IntegralE1InequalityTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E3,15)) {
                if (x >= Math.Log(Double.MaxValue / 10.0)) continue; // keep Exp(x) from overflowing
                double lower = Math.Log(1.0 + 2.0 / x) / 2.0;
                double upper = Math.Log(1.0 + 1.0 / x);
                double value = Math.Exp(x) * AdvancedMath.IntegralE(1, x);
                Console.WriteLine("{0}: {1} {2} {3}", x, lower, value, upper);
                Assert.IsTrue((lower < value) && (value < upper));
            }
        }

        [TestMethod]
        public void IntegralEInequalityTest () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1,100,5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E3,15)) {
                    if (x >= Math.Log(Double.MaxValue / 10.0)) continue; // keep Exp(x) from overflowing
                    double lower = 1.0 / (x + n);
                    double upper = 1.0 / (x + n - 1);
                    double value = Math.Exp(x) * AdvancedMath.IntegralE(n, x);
                    Console.WriteLine("n={0} x={1} {2} {3} {4}", n, x, lower, value, upper);
                    Assert.IsTrue((lower < value) && (value <= upper));
                }
            }
        }


        [TestMethod]
        public void IntegralERecurrenceTest () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E3, 15)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(n * AdvancedMath.IntegralE(n + 1, x) + x * AdvancedMath.IntegralE(n, x), Math.Exp(-x)));
                }
            }
        }

        [TestMethod]
        public void IntegralEIntegralTest () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 10, 2)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 10.0, 3)) {
                    Function<double, double> f = delegate(double t) {
                        return (Math.Exp(-x * t) / Math.Pow(t, n));
                    };
                    Interval r = Interval.FromEndpoints(1.0, Double.PositiveInfinity);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.IntegralE(n, x), FunctionMath.Integrate(f, r)));
                }
            }
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void IntegralEInvalidOrderTest () {
            AdvancedMath.IntegralE(-1, 1.1);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void IntegralEInvalidArgumentTest () {
            AdvancedMath.IntegralE(1, -1.1);
        }



        // Riemann Zeta

        [TestMethod]
        public void RiemannZetaSpecialCaseTest () {
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.RiemannZeta(-3.0), 1.0 / 120.0));
            //Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.RiemannZeta(-2.0), 0.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.RiemannZeta(-1.0), -1.0 / 12.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.RiemannZeta(0.0), -0.5));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.RiemannZeta(2.0), Math.PI * Math.PI / 6.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.RiemannZeta(4.0), Math.Pow(Math.PI, 4) / 90.0));
        }

        [TestMethod]
        public void RiemannZetaReflectionTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2,1.0E2,15)) {
                Console.WriteLine("x={0}", x);
                double zx = AdvancedMath.RiemannZeta(x);
                double gx = AdvancedMath.Gamma(x);
                double zr = AdvancedMath.RiemannZeta(1.0 - x);
                Console.WriteLine("  {0} vs. {1}", zr, 2.0 * Math.Pow(2.0 * Math.PI, -x) * Math.Cos(Math.PI * x / 2.0) * gx * zx);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(zr, 2.0 * Math.Pow(2.0 * Math.PI, -x) * Math.Cos(Math.PI * x / 2.0) * gx * zx));
            }
        }


        [TestMethod]
        public void ReimannZetaPrimesTest () {
            // pick high enough values so that p^-x == 1 within double precision before we reach the end of our list of primes
            foreach (double x in TestUtilities.GenerateRealValues(10.0, 1.0E3, 5)) {
                double f = 1.0;
                foreach (int p in primes) {
                    double t = 1.0 - Math.Pow(p, -x);
                    if (t == 1.0) break;
                    f = f * t;
                }
                Assert.IsTrue(TestUtilities.IsNearlyEqual(1.0 / AdvancedMath.RiemannZeta(x), f));
            }
        }

        private int[] primes = new int[] { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };

        [TestMethod]
        public void DirichletEtaSpecialCaseTest () {
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.DirichletEta(1.0), Math.Log(2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.DirichletEta(2.0), Math.PI * Math.PI / 12.0));
        }

        // Dawson

        // decide range of Dawson
        /*
        [TestMethod]
        public void DawsonReflectionTest () {
            foreach (double x in arguments) {
                Assert.IsTrue(AdvancedMath.Dawson(-x) == -AdvancedMath.Dawson(x));
            }
        }
        */

        [TestMethod]
        public void DawsonSpecialCaseTest () {
            Assert.IsTrue(AdvancedMath.Dawson(0.0) == 0.0);
        }

        [TestMethod]
        public void DawsonInequalityTest () {
            // this is a pretty lame inequality
            foreach (double x in arguments) {
                double F = AdvancedMath.Dawson(x);
                Assert.IsTrue((0 < F) && (F < 0.5410442246));
            }
        }

        [TestMethod]
        public void DawsonIntegralTest () {
            Function<double, double> f = delegate(double t) {
                return( Math.Exp(t * t) );
            };
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 4)) {
                Interval r = Interval.FromEndpoints(0.0, x);
                double DF = AdvancedMath.Dawson(x);
                double DI = Math.Exp(-x * x) * FunctionMath.Integrate(f, r);
                Console.WriteLine("{0} {1} {2}", x, DF, DI);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(DF, DI));
            }
        }
        // detects error in dawson function
        // x=2.51786267356676 DF=0.225166478728979 DI=0.221042742228873
        // integral is right!

        [TestMethod]
        public void FresnelReflectionTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 10)) {
                Console.WriteLine("x={0}", x);
                Assert.IsTrue(AdvancedMath.Fresnel(-x) == -AdvancedMath.Fresnel(x));
            }
        }

        [TestMethod]
        public void FresnelSIntegralTest () {
            Function<double, double> f = delegate(double t) {
                return ( Math.Sin(Math.PI / 2.0 * t * t) );
            };
            // if x gets too high, the integral has too many oscilations to converge
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E1, 3)) {
                Interval r = Interval.FromEndpoints(0.0, x);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.FresnelS(x), FunctionMath.Integrate(f, r)));
            }
        }

        [TestMethod]
        public void FresnelCIntegralTest () {
            Function<double, double> f = delegate(double t) {
                return (Math.Cos(Math.PI / 2.0 * t * t));
            };
            // if x gets too high, the integral has too many oscilations to converge
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E1, 3)) {
                Interval r = Interval.FromEndpoints(0.0, x);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.FresnelC(x), FunctionMath.Integrate(f, r)));
            }
        }

        [TestMethod]
        public void IntegralSiIntegralTest () {
            Function<double, double> f = delegate(double t) {
                return (Math.Sin(t) / t);
            };
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 5)) {
                Interval r = Interval.FromEndpoints(0.0, x);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.IntegralSi(x), FunctionMath.Integrate(f, r)));
            }
        }
        /*
        // nonconvergent
        [TestMethod]
        public void IntegralCiIntegralTest () {
            Function<double, double> f = delegate(double t) {
                return (Math.Cos(t) / t);
            };
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 5)) {
                Interval r = Interval.FromEndpoints(x, Double.PositiveInfinity);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.IntegralCi(x), -FunctionMath.Integrate(f, r)));
            }
        }
        */

        [TestMethod]
        public void LambertTest () {

            foreach (double x in TestUtilities.GenerateUniformRealValues(-1.0 / Math.E, 1.0, 6)) {
                double W = AdvancedMath.LambertW(x);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(W * Math.Exp(W), x));
            }

            foreach (double x in TestUtilities.GenerateRealValues(1.0, 1.0E3, 3)) {
                double W = AdvancedMath.LambertW(x);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(W * Math.Exp(W), x));
            }

        }

        [TestMethod]
        public void LambertSpecialCaseTest () {

            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.LambertW(-1.0 / Math.E), -1.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.LambertW(-Math.Log(2.0) / 2.0), -Math.Log(2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.LambertW(0.0), 0.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.LambertW(Math.E), 1.0));

        }


        [TestMethod]
        public void CoulombEtaZeroTest () {

            // can't go too far out because Math.Sin becomes unreliable (before our functions do)
            foreach (double rho in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 10)) {

                double F = AdvancedMath.CoulombF(0, 0.0, rho);
                Console.WriteLine("{0} {1} {2}", rho, F, Math.Sin(rho));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(F, Math.Sin(rho)));

                double G = AdvancedMath.CoulombG(0, 0.0, rho);
                Console.WriteLine("{0} {1} {2}", rho, G, Math.Cos(rho));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(G, Math.Cos(rho)));

            }

        }

        [TestMethod]
        public void CoulombWronskianTest () {

            foreach (int L in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double eta in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 10)) {
                    foreach (double rho in TestUtilities.GenerateRealValues(1.0E-3, 1.0E3, 20)) {

                        Console.WriteLine("L={0} eta={1} rho={2}", L, eta, rho);

                        CoulombWronskianHelper(L, eta, rho);
                        CoulombWronskianHelper(L, -eta, rho);

                    }
                }
            }
        }

        private static void CoulombWronskianHelper (int L, double eta, double rho) {

            double FM = AdvancedMath.CoulombF(L - 1, eta, rho);
            double GM = AdvancedMath.CoulombG(L - 1, eta, rho);

            double F = AdvancedMath.CoulombF(L, eta, rho);
            double G = AdvancedMath.CoulombG(L, eta, rho);

            //Console.WriteLine("FM={0} GM={1} F={2} G={3}", FM, GM, F, G);
            //Console.WriteLine("FM * G - F * GM = {0}", FM * G - F * GM);
            //Console.WriteLine("W = {0}", L / Math.Sqrt(L * L + eta * eta));

            Assert.IsTrue(TestUtilities.IsSumNearlyEqual(FM * G, -F * GM, L / Math.Sqrt(L * L + eta * eta), 8.0 * TestUtilities.TargetPrecision));

        }

        [TestMethod]
        public void CoulombRecursionTest () {
            
            foreach (int L in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double eta in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 10)) {
                    foreach (double rho in TestUtilities.GenerateRealValues(1.0E-3, 1.0E3, 20)) {

                        Console.WriteLine("L={0} eta={1} rho={2}", L, eta, rho);

                        CoulombRecursionTestHelper(L, eta, rho,
                            AdvancedMath.CoulombF(L - 1, eta, rho), AdvancedMath.CoulombF(L, eta, rho), AdvancedMath.CoulombF(L + 1, eta, rho)
                        );

                        CoulombRecursionTestHelper(L, -eta, rho,
                            AdvancedMath.CoulombF(L - 1, -eta, rho), AdvancedMath.CoulombF(L, -eta, rho), AdvancedMath.CoulombF(L + 1, -eta, rho)
                        );

                        CoulombRecursionTestHelper(L, eta, rho,
                            AdvancedMath.CoulombG(L - 1, eta, rho), AdvancedMath.CoulombG(L, eta, rho), AdvancedMath.CoulombG(L + 1, eta, rho)
                        );

                        CoulombRecursionTestHelper(L, -eta, rho,
                            AdvancedMath.CoulombG(L - 1, -eta, rho), AdvancedMath.CoulombG(L, -eta, rho), AdvancedMath.CoulombG(L + 1, -eta, rho)
                        );

                    }
                }
            }
            

        }

        private static void CoulombRecursionTestHelper (double L, double eta, double rho, double UM, double U, double UP) {

            double am = (L + 1) * Math.Sqrt(L * L + eta * eta);
            double a = (2 * L + 1) * (eta + L * (L + 1) / rho);
            double ap = L * Math.Sqrt((L + 1) * (L + 1) + eta * eta);

            Assert.IsTrue(TestUtilities.IsSumNearlyEqual(am * UM, ap * UP, a * U, 4.0 * TestUtilities.TargetPrecision));

        }

        [TestMethod]
        public void ModifiedBesselArgumentZeroTest () {

            Assert.IsTrue(AdvancedMath.ModifiedBesselI(0.0, 0.0) == 1.0);
            Assert.IsTrue(AdvancedMath.ModifiedBesselI(1.0, 0.0) == 0.0);
            Assert.IsTrue(AdvancedMath.ModifiedBesselK(0.0, 0.0) == Double.PositiveInfinity);
            //Assert.IsTrue(AdvancedMath.ModifiedBesselK(1.0, 0.0) == Double.PositiveInfinity);

        }

        [TestMethod]
        public void ModifiedBesselWronskianTest () {

            foreach (double nu in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 8)) {

                    Console.WriteLine("nu={0} x={1}", nu, x);

                    double I = AdvancedMath.ModifiedBesselI(nu, x);
                    double Ip1 = AdvancedMath.ModifiedBesselI(nu + 1.0, x);
                    double K = AdvancedMath.ModifiedBesselK(nu, x);
                    double Kp1 = AdvancedMath.ModifiedBesselK(nu + 1.0, x);

                    Console.WriteLine("{0} ?= {1}", I * Kp1 + Ip1 * K, 1.0 / x);

                    // no need to use IsSumNearlyEqual; both terms are positive so no cancelation is possible
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(I * Kp1 + Ip1 * K, 1.0 / x));

                }
            }

        }

        [TestMethod]
        public void ModifiedBesselHalfIntegerOrderTest () {

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 8)) {
                //Console.WriteLine("x={0}", x);

                double F = Math.Sqrt(Math.PI / 2.0 / x);

                //Console.WriteLine("{0} ?= {1}", F * AdvancedMath.ModifiedBesselI(0.5, x), Math.Sinh(x) / x);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(F * AdvancedMath.ModifiedBesselI(0.5, x), Math.Sinh(x) / x));

                //Console.WriteLine("{0} ?= {1}", F * AdvancedMath.ModifiedBesselI(-0.5, x), Math.Cosh(x) / x);
                //Assert.IsTrue(TestUtilities.IsNearlyEqual(F * AdvancedMath.ModifiedBesselI(-0.5, x), Math.Cosh(x) / x));

                //Console.WriteLine("{0} ?= {1}", AdvancedMath.ModifiedBesselK(0.5, x), F * Math.Exp(-x));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.ModifiedBesselK(0.5, x), F * Math.Exp(-x)));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.ModifiedBesselK(1.5, x), F * Math.Exp(-x) * (1.0 + 1.0 / x)));

            }

        }

        [TestMethod]
        public void ModifiedBesselIntegralTest () {

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 5)) {
                //Console.WriteLine("x={0}", x);

                Function<double, double> f = delegate (double t) {
                    return (Math.Exp(x * Math.Cos(t)));
                };
                Interval r = Interval.FromEndpoints(0.0, Math.PI);
                double I = FunctionMath.Integrate(f, r);

                //Console.WriteLine("{0} ?= {1}", I / Math.PI, AdvancedMath.ModifiedBesselI(0.0, x));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(I / Math.PI, AdvancedMath.ModifiedBesselI(0.0, x)));

            }

        }

        [TestMethod]
        public void AiryIntegralTest () {

            Function<double, double> f = delegate(double t) {
                return (AdvancedMath.AiryAi(t));
            };
            Interval r = Interval.FromEndpoints(0.0, Double.PositiveInfinity);
            double I = FunctionMath.Integrate(f, r);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(I, 1.0 / 3.0));

        }

        [TestMethod]
        public void AiryZeroArgumentTest () {

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.AiryAi(0.0), 1.0 / Math.Pow(3.0, 2.0 / 3.0) / AdvancedMath.Gamma(2.0 / 3.0)
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.AiryBi(0.0), 1.0 / Math.Pow(3.0, 1.0 / 6.0) / AdvancedMath.Gamma(2.0 / 3.0)
            ));
        }

        [TestMethod]
        public void GoldenRatioTest () {

            Assert.IsTrue(TestUtilities.IsNearlyEqual((1.0 + AdvancedMath.GoldenRatio) / AdvancedMath.GoldenRatio, AdvancedMath.GoldenRatio));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(1.0 / AdvancedMath.GoldenRatio, AdvancedMath.GoldenRatio - 1.0));
        }

        [TestMethod]
        public void CatalanIntegralTest () {

            Function<double, double> f1 = delegate(double t) {
                return (t / Math.Sin(t) / Math.Cos(t));
            };
            Interval r1 = Interval.FromEndpoints(0.0, Math.PI / 4.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(FunctionMath.Integrate(f1, r1), AdvancedMath.Catalan));

            Function<double, double> f2 = delegate(double t) {
                return (-Math.Log(t) / (1 + t * t));
            };
            Interval r2 = Interval.FromEndpoints(0.0, 1.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(FunctionMath.Integrate(f2, r2), AdvancedMath.Catalan));

        }

        [TestMethod]
        public void DiLogSpecialCaseTest () {

            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.DiLog(-1.0), -Math.PI * Math.PI / 12.0));
            Assert.IsTrue(AdvancedMath.DiLog(0.0) == 0.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.DiLog(0.5), Math.PI * Math.PI / 12.0 - Math.Log(2.0) * Math.Log(2.0) / 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.DiLog(1.0 / AdvancedMath.GoldenRatio), Math.PI * Math.PI / 10.0 - Math.Log(AdvancedMath.GoldenRatio) * Math.Log(AdvancedMath.GoldenRatio)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.DiLog(1.0), Math.PI * Math.PI / 6.0));

        }

        [TestMethod]
        public void DiLogDuplicationTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0, 10)) {
                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                    AdvancedMath.DiLog(x), AdvancedMath.DiLog(-x), AdvancedMath.DiLog(x * x) / 2.0
                ));
            }
        }

        [TestMethod]
        public void DiLogBailyIdentityTest () {

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                36.0 * AdvancedMath.DiLog(1.0 / 2.0) - 36.0 * AdvancedMath.DiLog(1.0 / 4.0) -
                12.0 * AdvancedMath.DiLog(1.0 / 8.0) + 6.0 * AdvancedMath.DiLog(1.0 / 64.0),
                Math.PI * Math.PI
            ));

        }

        [TestMethod]
        public void DiLogExpIntegralTest () {

            foreach (double x in TestUtilities.GenerateRealValues(0.1,10.0,5)) {

                Function<double, double> f = delegate (double t) {
                    return (t / (Math.Exp(t) + x));
                };
                Interval r = Interval.FromEndpoints(0.0, Double.PositiveInfinity);
                double I = FunctionMath.Integrate(f, r);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.DiLog(-x), -x * I));

            }

        }

        [TestMethod]
        public void DiLogLogIntegralTest () {

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 5)) {

                Function<double, double> f = delegate(double t) {
                    return (Math.Log(t) / t / (t - x));
                };
                Interval r = Interval.FromEndpoints(1.0, Double.PositiveInfinity);
                double I = FunctionMath.Integrate(f, r);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.DiLog(x), x * I));

            }

        }

        [TestMethod]
        public void EllipticKSpecialCases () {

            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticK(0.0), Math.PI / 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticK(1.0 / Math.Sqrt(2.0)), Math.Pow(AdvancedMath.Gamma(1.0 / 4.0), 2.0) / Math.Sqrt(Math.PI) / 4.0));

        }

        [TestMethod]
        public void EllipticKInequality () {
            foreach (double k in TestUtilities.GenerateUniformRealValues(0.0, 1.0, 4)) {
                double S = AdvancedMath.EllipticK(k) + Math.Log(Math.Sqrt(1.0 - k * k));
                Assert.IsTrue((Math.Log(4.0) <= S) && (S <= Math.PI / 2.0));
            }
        }

        [TestMethod]
        public void EllipticKIntegral () {

            Interval i = Interval.FromEndpoints(0.0, Math.PI / 2.0);

            foreach (double k in TestUtilities.GenerateRealValues(0.01, 1.0, 4)) {

                Function<double, double> f = delegate(double t) {
                    double z = k * Math.Sin(t);
                    return (1.0 / Math.Sqrt(1.0 - z * z));
                };

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    FunctionMath.Integrate(f, i), AdvancedMath.EllipticK(k)
                ));

            }

        }

        [TestMethod]
        public void EllipticKCatalanIntegral () {

            System.Diagnostics.Stopwatch timer = System.Diagnostics.Stopwatch.StartNew();


            Interval i = Interval.FromEndpoints(0.0, 1.0);

            Function<double, double> f1 = delegate(double k) {
                return (AdvancedMath.EllipticK(k));
            };
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate(f1, i), 2.0 * AdvancedMath.Catalan
            ));

            Function<double, double> f2 = delegate(double k) {
                return (AdvancedMath.EllipticK(k) * k);
            };
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate(f2, i), 1.0
            ));

            timer.Stop();
            Console.WriteLine(timer.ElapsedTicks);

        }

        [TestMethod]
        public void EllipticFSpecialCases () {
            foreach (double phi in TestUtilities.GenerateUniformRealValues(0.0, Math.PI / 2.0, 4)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticF(phi, 0.0), phi));
            }
        }

        [TestMethod]
        public void EllipticFIntegration () {

            foreach (double k in TestUtilities.GenerateRealValues(0.01, 1.0, 4)) {

                Function<double, double> f = delegate(double t) {
                    double z = k * Math.Sin(t);
                    return (1.0 / Math.Sqrt(1.0 - z * z));
                };

                foreach (double phi in TestUtilities.GenerateUniformRealValues(0.0, Math.PI/2.0, 4)) {

                    Interval i = Interval.FromEndpoints(0.0, phi);

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        FunctionMath.Integrate(f, i), AdvancedMath.EllipticF(phi, k)
                    ));

                }

            }

        }

        [TestMethod]
        public void CarlsonFDuplication () {

            double[] args = TestUtilities.GenerateRealValues(0.01, 100.0, 5);
            for (int i = 0; i < args.Length; i++) {
                double x = args[i];
                for (int j = 0; j <= i; j++) {
                    double y = args[j];
                    for (int k = 0; k <= j; k++) {
                        double z = args[k];

                        double a = Math.Sqrt(x * y) + Math.Sqrt(x * z) + Math.Sqrt(y * z);
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(
                            AdvancedMath.CarlsonF(x, y, z), 2.0 * AdvancedMath.CarlsonF(x+a, y+a, z+a)
                        ));

                    }
                }
            }
        }

        [TestMethod]
        public void CarlsonSpecialCases () {

            foreach (double x in TestUtilities.GenerateRealValues(0.01, 1.0E4, 8)) {

                // triply degenerate case
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.CarlsonF(x, x, x), 1.0 / Math.Sqrt(x)));


            }

        }


        [TestMethod]
        public void ETest () {
            double k = 1.0E-4;
            Console.WriteLine(AdvancedMath.EllipticK(k));
            //Console.WriteLine(AdvancedMath.EllipticK_Series(0.5));
            Console.WriteLine(AdvancedMath.CarlsonF(0.0, 1.0 - k * k, 1.0));

        }



#if FUTURE

        [TestMethod]
        public void ModifiedBesselSpecialCaseTest () {
            Assert.IsTrue(AdvancedMath.BesselI(0, 0.0) == 1.0);
            Assert.IsTrue(AdvancedMath.BesselI(1, 0.0) == 0.0);
            Assert.IsTrue(AdvancedMath.BesselI(10, 0.0) == 0.0);
        }

        [TestMethod]
        public void ModifiedBesselRecurrenceTest () {
            foreach (int n in TestUtilities.GenerateIntegerValues(0, 2, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(0.03, 300.0, 4)) {
                    // don't let x get too big or I(x) will explode

                    double IM = AdvancedMath.BesselI(n - 1, x);
                    double I0 = AdvancedMath.BesselI(n, x);
                    double IP = AdvancedMath.BesselI(n + 1, x);

                    Console.WriteLine("{0} {1} IM={2} I0={3} IP={4}", n, x, IM, I0, IP);
                    Console.WriteLine("  {0} ?= {1}", IM, (2.0 * n / x) * I0 + IP);

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(IM, (2.0 * n / x) * I0 + IP));
                }
            }
        }

        [TestMethod]
        public void ModifiedBesselTowerTest () {

            // x shouldn't be too big or we will need too many terms before the series converges
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 10)) {

                // S1 = I0 + 2 I1 + 2 I2 + 2 I3 + ... = e^x
                // S2 = I0 + 2 I2 + 2 I4 + 2 I8 + ... = cosh x
                // S3 = 2 I1 + 2 I3 + 2 I5 + 2 I7 + ... = sinh x 

                double I0 = AdvancedMath.BesselI(0, x);

                double s1 = I0;
                double s2 = I0;
                double s3 = 0.0;

                int n = 1;
                while (true) {

                    double s1_old = s1;
                    double s2_old = s2;
                    double s3_old = s3;

                    double I1 = AdvancedMath.BesselI(2 * n - 1, x);
                    double I2 = AdvancedMath.BesselI(2 * n, x);

                    s1 += 2.0 * (I1 + I2);
                    s2 += 2.0 * I2;
                    s3 += 2.0 * I1;

                    if ((s1 == s1_old) && (s2 == s2_old) && (s3 == s3_old)) break;

                    n++;
                    if (n > 100) throw new NonconvergenceException();

                }
                Console.WriteLine("x={0}, n={1}", x, n);

                Console.WriteLine("{0} {1}", s1, Math.Exp(x));
                Console.WriteLine("{0} {1}", s2, Math.Cosh(x));
                Console.WriteLine("{0} {1}", s3, Math.Sinh(x));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(s1, Math.Exp(x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s2, Math.Cosh(x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s3, Math.Sinh(x)));

            }

        }

        [TestMethod]
        public void ModifiedBesselI0IntegralTest () {
            foreach (double x in TestUtilities.GenerateRealValues(0.003, 300.0, 6)) {
                // don't let x get too big or I(x) will explode
                Function<double,double> f = delegate (double t) {
                    return( Math.Cosh(x * Math.Cos(t) ) );
                };
                Interval r = Interval.FromEndpoints(0.0, Math.PI);
                double I = FunctionMath.Integrate(f, r) / Math.PI;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.BesselI(0, x), I));
            }
        }

        [TestMethod]
        public void ModifiedBesselInIntegralTest () {
            foreach (int n in TestUtilities.GenerateIntegerValues(0, 2, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(0.03, 300.0, 4)) {
                    // don't let x get too big or I(x) will explode
                    Function<double, double> f = delegate(double t) {
                        return (Math.Exp(x * Math.Cos(t)) * Math.Cos(n * t));
                    };
                    Interval r = Interval.FromEndpoints(0.0, Math.PI);
                    double I = FunctionMath.Integrate(f, r) / Math.PI;
                    Console.WriteLine("{0} {1} {2} {3}", n, x, I, AdvancedMath.BesselI(n, x));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.BesselI(n, x), I));
                }
            }
        }

#endif

    }
}
