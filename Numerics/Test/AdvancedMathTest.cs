using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Test {
    
    [TestClass]
    public class AdvancedMathTest {

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
        public void IntegerBesselSpecialCase () {
            Assert.IsTrue(AdvancedMath.BesselJ(0, 0.0) == 1.0);
            Assert.IsTrue(AdvancedMath.BesselJ(1, 0.0) == 0.0);
            Assert.IsTrue(AdvancedMath.BesselY(0, 0.0) == Double.NegativeInfinity);
            Assert.IsTrue(AdvancedMath.BesselY(1, 0.0) == Double.NegativeInfinity);
        }

        [TestMethod]
        public void IntegerBesselNegativeOrder () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 16)) {
                    int s = 1 - 2*(n % 2);
                    Assert.IsTrue(AdvancedMath.BesselJ(-n, x) == s * AdvancedMath.BesselJ(n, x));
                    Assert.IsTrue(AdvancedMath.BesselY(-n, x) == s * AdvancedMath.BesselY(n, x));
                }
            }
        }

        [TestMethod]
        public void IntegerBesselNegativeArgument () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 16)) {
                    int s = (n % 2 == 0) ? 1 : -1;
                    Assert.IsTrue(AdvancedMath.BesselJ(n, -x) == s * AdvancedMath.BesselJ(n, x));
                }
            }
        }

        [TestMethod]
        public void IntegerBesselRecurrence () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1,100,4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E4,16)) {
                    Console.WriteLine("n = {0}, x={1}", n, x);
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        AdvancedMath.BesselJ(n - 1, x), AdvancedMath.BesselJ(n + 1, x), 2 * n / x * AdvancedMath.BesselJ(n, x)
                    ));
                    if (!BesselYInRange(n, x)) continue;
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        AdvancedMath.BesselY(n - 1, x), AdvancedMath.BesselY(n + 1, x), 2 * n / x * AdvancedMath.BesselY(n, x)
                    ));
                }
            }
        }

        [TestMethod]
        public void IntegerBesselWronskian () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 16)) {
                    if (BesselYInRange(n, x)) {
                        if ((x > n) || (Math.Exp(AdvancedMath.LogGamma(n) - n * Math.Log(x)) < Math.Sqrt(Double.MaxValue))) {
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(2 / Math.PI / x, AdvancedMath.BesselJ(n + 1, x) * AdvancedMath.BesselY(n, x) - AdvancedMath.BesselJ(n, x) * AdvancedMath.BesselY(n + 1, x)), String.Format("n={0}, x={1}, J(n,x)={2}, J(n+1,x) = {3}, Y(n,x)={4}, Y(n+1,x)={5}", n, x, AdvancedMath.BesselJ(n, x), AdvancedMath.BesselJ(n + 1, x), AdvancedMath.BesselY(n, x), AdvancedMath.BesselY(n + 1, x)));
                        }
                    }
                }
            }
        }

        [TestMethod]
        public void BesselJ0Integral () {
            // Abromowitz & Stegun 9.1.18
            // J_0(x) = \int_{0}^{\pi} \cos( x \sin(t) ) \, dt
            // don't let x get too big, or the integral becomes too oscilatory to do accurately
            Interval r = Interval.FromEndpoints(0.0, Math.PI);
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 8)) {
                Func<double, double> f = delegate(double t) {
                    return (Math.Cos(x * Math.Sin(t)));
                };
                double J0 = FunctionMath.Integrate(f, r) / Math.PI;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.BesselJ(0, x), J0));
            }
        }

        [TestMethod]
        public void BesselY0Integral () {
            // Abromowitz & Stegen 9.1.19
            Interval r = Interval.FromEndpoints(0.0, Math.PI / 2.0);
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 8)) {
                Func<double, double> f = delegate(double t) {
                    double s = Math.Sin(t);
                    return (Math.Cos(x * Math.Cos(t)) * (AdvancedMath.EulerGamma + Math.Log(2.0 * x * s * s)));
                };
                double Y0 = 4.0 / (Math.PI * Math.PI) * FunctionMath.Integrate(f, r);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.BesselY(0, x), Y0));
            }
        }

        [TestMethod]
        public void IntegerBesselJIntegral () {
            // Abromowitz & Stegun 9.1.21
            Interval r = Interval.FromEndpoints(0.0, Math.PI);
            foreach (double x in TestUtilities.GenerateRealValues(0.1, 10.0, 8)) {
                foreach (int n in TestUtilities.GenerateIntegerValues(1, 10, 4)) {
                    Func<double, double> f = delegate(double t) {
                        return (Math.Cos(x * Math.Sin(t) - n * t));
                    };
                    double J = FunctionMath.Integrate(f, r) / Math.PI;
                    Console.WriteLine("n={0} x={1} J={2} I={3}", n, x, AdvancedMath.BesselJ(n, x), J);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.BesselJ(n, x), J));
                }
            }
        }

        [TestMethod]
        public void BesselKapteynIntegral () {
            // don't pick x too big, so as not to have a problem with the inaccuracy of trig functions for large arguments
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 8)) {
                Func<double, double> f = delegate(double t) {
                    return (Math.Cos(x - t) * AdvancedMath.BesselJ(0, t));
                };
                Interval r = Interval.FromEndpoints(0.0, x);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(x * AdvancedMath.BesselJ(0, x), FunctionMath.Integrate(f, r)));
            }
        }

        [TestMethod]
        public void BesselLipshitzIntegral () {
            // \int_{0}^{\infty} J_0(x) \, dx = \frac{1}{\sqrt{2}}
            Func<double, double> f = delegate(double t) {
                return (Math.Exp(-t) * AdvancedMath.BesselJ(0, t));
            };
            Interval r = Interval.FromEndpoints(0.0, Double.PositiveInfinity);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(FunctionMath.Integrate(f, r), 1.0 / Math.Sqrt(2.0)));
        }

        [TestMethod]
        public void BesselWeberIntegral () {
            Func<double, double> f = delegate(double t) {
                return (Math.Exp(-t*t) * AdvancedMath.BesselJ(0, t) * t);
            };
            Interval r = Interval.FromEndpoints(0.0, Double.PositiveInfinity);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(FunctionMath.Integrate(f, r), Math.Exp(-1.0/4.0) / 2.0));
        }

        [TestMethod]
        public void IntegerBesselRealBesselAgreement () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1,100,8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E4,16)) {
                    //Console.WriteLine("n={0},x={1}", n, x);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.BesselJ(n, x),
                        AdvancedMath.BesselJ((double) n, x
                    )));
                    if (!BesselYInRange(n, x)) continue; // don't try to evaluate Y if it's too big
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.BesselY(n, x),
                        AdvancedMath.BesselY((double) n, x)
                    ));
                }
            }
        }

        [TestMethod]
        public void FullBesselRealBesselAgreement () {
            foreach (double nu in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-3, 1.0E4, 16)) {
                    if (!BesselYInRange(nu, x)) continue;
                    SolutionPair s = AdvancedMath.Bessel(nu, x);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        s.FirstSolutionValue, AdvancedMath.BesselJ(nu, x)
                    ));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        s.SecondSolutionValue, AdvancedMath.BesselY(nu, x)
                    ));
                }
            }
        }

        [TestMethod]
        public void FullBesselWronskian () {
            foreach (double nu in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-3, 1.0E4, 16)) {
                    if (!BesselYInRange(nu, x)) continue;
                    SolutionPair s = AdvancedMath.Bessel(nu, x);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        s.FirstSolutionValue * s.SecondSolutionDerivative - s.SecondSolutionValue * s.FirstSolutionDerivative,
                        2.0 / Math.PI / x
                    ));
                }
            }
        }

        [TestMethod]
        public void FullBesselDerivative () {
            // A&S 9.1.27
            foreach (double nu in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-3, 1.0E4, 16)) {
                    if (!BesselYInRange(nu, x)) continue;
                    SolutionPair s0 = AdvancedMath.Bessel(nu, x);
                    SolutionPair s1 = AdvancedMath.Bessel(nu + 1.0, x);
                    Console.WriteLine("nu={0} x={1} J={2} Jp1={3} J'={4}", nu, x, s0.FirstSolutionValue, s1.FirstSolutionValue, s0.FirstSolutionDerivative);
                    Console.WriteLine("                           J'={0}", nu / x * s0.FirstSolutionValue - s1.FirstSolutionValue);
                    
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        nu / x * s0.FirstSolutionValue, -s1.FirstSolutionValue, s0.FirstSolutionDerivative
                    ));
                    
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        nu / x * s0.SecondSolutionValue, -s1.SecondSolutionValue, s0.SecondSolutionDerivative
                    ));
                    
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        s0.FirstSolutionValue, -(nu + 1.0) / x * s1.FirstSolutionValue, s1.FirstSolutionDerivative
                    ));
                    
                    /*
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        s0.SecondSolutionValue, -(nu + 1.0) / x * s1.SecondSolutionValue, s1.SecondSolutionDerivative
                    ));
                    */
                }
            }
        }

        [TestMethod]
        public void FullModifiedBesselDerivitive () {
            // A&S 9.6.26
            foreach (double nu in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 8)) {

                    SolutionPair s0 = AdvancedMath.ModifiedBessel(nu, x);
                    SolutionPair s1 = AdvancedMath.ModifiedBessel(nu + 1.0, x);
                    
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        s1.FirstSolutionValue + nu / x * s0.FirstSolutionValue, s0.FirstSolutionDerivative
                    ));
                                        
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        -s1.SecondSolutionValue, nu / x * s0.SecondSolutionValue, s0.SecondSolutionDerivative
                    ));
                    
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        s0.FirstSolutionValue, -(nu + 1.0) / x * s1.FirstSolutionValue, s1.FirstSolutionDerivative
                    ));
                    
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        s0.SecondSolutionValue + (nu + 1.0) / x * s1.SecondSolutionValue, -s1.SecondSolutionDerivative
                    ));
                    
                    // when signs agree, we can just test directly; when signs disagree, we test sum, which allows for cancelation

                }
            }
        }

        [TestMethod]
        public void RealBesselWronskian () {
            foreach (double nu in orders) {
                foreach (double x in arguments) {
                    if (BesselYInRange(nu, x)) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(
                            2 / Math.PI / x,
                            AdvancedMath.BesselJ(nu + 1.0, x) * AdvancedMath.BesselY(nu, x) - AdvancedMath.BesselJ(nu, x) * AdvancedMath.BesselY(nu + 1.0, x)
                        ));
                    }
                }
            }
        }

        [TestMethod]
        public void RealBesselRecurrence () {
            foreach (double nu in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 16)) {
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        AdvancedMath.BesselJ(nu - 1.0, x), AdvancedMath.BesselJ(nu + 1.0, x),
                        2.0 * nu / x * AdvancedMath.BesselJ(nu, x)
                    ));
                    if (BesselYInRange(nu, x)) {
                        Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                            AdvancedMath.BesselY(nu - 1.0, x), AdvancedMath.BesselY(nu + 1.0, x),
                            2.0 * nu / x * AdvancedMath.BesselY(nu, x)
                        ));
                    }
                }
            }
        }

        [TestMethod]
        public void RealBesselInequality () {
            foreach (double nu in TestUtilities.GenerateRealValues(0.5, 50.0, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(nu, 1000.0, 8)) {
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
        public void RealBesselFresnel () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2,1.0E2,8)) {
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
        public void BesselTower () {

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
        public void RealBesselJIntegral () {
            foreach (double nu in TestUtilities.GenerateRealValues(0.1, 10.0, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(0.01, 100.0, 4)) {
                    Func<double, double> f = delegate(double t) {
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
        public void SphericalBesselSpecialCase () {
            Assert.IsTrue(AdvancedMath.SphericalBesselJ(0, 0.0) == 1.0);
            Assert.IsTrue(AdvancedMath.SphericalBesselJ(1, 0.0) == 0.0);
            Assert.IsTrue(AdvancedMath.SphericalBesselY(0, 0.0) == Double.NegativeInfinity);
            Assert.IsTrue(AdvancedMath.SphericalBesselY(1, 0.0) == Double.NegativeInfinity);
        }

        [TestMethod]
        public void SphericalBesselRecurrence () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1,100,5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E4,15)) {
                    Console.WriteLine("n={0} x={1}", n, x);
                    double jm1 = AdvancedMath.SphericalBesselJ(n - 1, x);
                    double jp1 = AdvancedMath.SphericalBesselJ(n + 1, x);
                    double j = AdvancedMath.SphericalBesselJ(n, x);
                    //Console.WriteLine("  {0:g16} + {1:g16} = {2:g16} ?= {3:g16}", jm1, jp1, jm1 + jp1, (2 * n + 1) * j / x);
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(jm1, jp1, (2*n+1)*j/x));
                    if (BesselYInRange(n, x)) {
                        //double ym1 = AdvancedMath.SphericalBesselY(n - 1, x);
                        //double yp1 = AdvancedMath.SphericalBesselY(n + 1, x);
                        //double y = AdvancedMath.SphericalBesselY(n, x);
                        //Console.WriteLine("  {0} + {1} = {2} ?= {3}", ym1, yp1, ym1 + yp1, (2 * n + 1) * y / x);
                        Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                            AdvancedMath.SphericalBesselY(n - 1, x),
                            AdvancedMath.SphericalBesselY(n + 1, x), (2 * n + 1) / x * AdvancedMath.SphericalBesselY(n, x)
                        ));
                    }

                }
            }
        }

        [TestMethod]
        public void SphericalBesselNegativeOrder () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 15)) {
                    int s = 2 * (n % 2) - 1;
                    Assert.IsTrue(AdvancedMath.SphericalBesselY(n, x) == s * AdvancedMath.SphericalBesselJ(-n - 1, x));
                    Assert.IsTrue(AdvancedMath.SphericalBesselY(-n, x) == s * AdvancedMath.SphericalBesselJ(n - 1, x));
                }
            }
        }

        [TestMethod]
        public void SphericalBesselWronskian () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 15)) {
                    if (BesselYInRange(n, x)) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.SphericalBesselJ(n, x) * AdvancedMath.SphericalBesselY(n - 1, x) - AdvancedMath.SphericalBesselJ(n - 1, x) * AdvancedMath.SphericalBesselY(n, x), 1.0 / (x * x)));
                    }
                }
            }
        }

        [TestMethod]
        public void SphericalBesselRealBesselAgreement () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 16)) {
                    Console.WriteLine("n={0} x={1}", n, x);
                     Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.SphericalBesselJ(n, x),
                        Math.Sqrt(Math.PI / 2.0 / x) * AdvancedMath.BesselJ(n + 0.5, x)
                    ));
                    if (BesselYInRange(n, x)) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(
                            AdvancedMath.SphericalBesselY(n, x),
                            Math.Sqrt(Math.PI / 2.0 / x) * AdvancedMath.BesselY(n + 0.5, x)
                         ));
                    }
                }
            }
        }

        [TestMethod]
        public void SphericalBesselTower () {

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
        public void GammaRecurrsion () {
            // limit x to avoid overflow
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E2,16)) {
                Console.WriteLine("x={0}, G={1}", x, AdvancedMath.Gamma(x));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(x + 1.0), x * AdvancedMath.Gamma(x)));
            }
        }

        [TestMethod]
        public void GammaSpecialCases () {
            // it would be nice to be able to make some of these comparisons exact
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(0.5), Math.Sqrt(Math.PI)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(1.0), 1.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(1.5), Math.Sqrt(Math.PI) / 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(2.0), 1.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(3.0), 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(4.0), 6.0));
        }

        [TestMethod]
        public void GammaReflection () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E2,16)) {
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
        public void GammaInequality () {
            foreach (double x in TestUtilities.GenerateRealValues(2.0,1.0E2,16)) {
                // for x >= 2
                double lower = Math.Pow(x / Math.E, x - 1.0);
                double upper = Math.Pow(x / 2.0, x - 1.0);
                double value = AdvancedMath.Gamma(x);
                Assert.IsTrue((lower <= value) && (value <= upper));
            }
        }

        [TestMethod]
        public void GammaIntegral () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0, 1.0E2, 4)) {
                Func<double,double> f = delegate (double t) {
                    return( Math.Pow(t, x - 1.0) * Math.Exp(-t) );
                };
                Interval r = Interval.FromEndpoints(0.0, Double.PositiveInfinity);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Gamma(x), FunctionMath.Integrate(f, r)));
            }
        }

        [TestMethod]
        public void GammaTrottIdentity () {

            double G1 = AdvancedMath.Gamma(1.0 / 24.0);
            double G5 = AdvancedMath.Gamma(5.0 / 24.0);
            double G7 = AdvancedMath.Gamma(7.0 / 24.0);
            double G11 = AdvancedMath.Gamma(11.0 / 24.0);

            Assert.IsTrue(TestUtilities.IsNearlyEqual((G1 * G11) / (G5 * G7), Math.Sqrt(3.0) * Math.Sqrt(2.0 + Math.Sqrt(3.0)))); 
        }

        [TestMethod]
        public void LogGammaDuplication () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E4,24)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.LogGamma(2.0 * x),
                    AdvancedMath.LogGamma(x) + AdvancedMath.LogGamma(x + 0.5) + (2.0 * x - 0.5) * Math.Log(2.0) - 0.5 * Math.Log(2.0 * Math.PI)
                ));
            }
        }

        [TestMethod]
        public void LogGammaTriplication () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 24)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.LogGamma(3.0 * x),
                    AdvancedMath.LogGamma(x) + AdvancedMath.LogGamma(x + 1.0 / 3.0) + AdvancedMath.LogGamma(x + 2.0 / 3.0) + (3.0 * x - 0.5) * Math.Log(3.0) - Math.Log(2.0 * Math.PI)
                ));
            }
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void LogGammaNegativeArgument () {
            AdvancedMath.LogGamma(-0.5);
        }

        [TestMethod]
        public void GammaRatioInequality () {

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 16)) {

                double LB = Math.Log(Math.Sqrt(x + 0.25));
                double UB = Math.Log((x + 0.5) / Math.Sqrt(x + 0.75));
                double R = AdvancedMath.LogGamma(x + 1.0) - AdvancedMath.LogGamma(x + 0.5);

                Assert.IsTrue(LB <= R);
                Assert.IsTrue(R <= UB);

            }

        }

        [TestMethod]
        public void PsiSpecialCases () {
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(1.0), -AdvancedMath.EulerGamma));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(2.0), -AdvancedMath.EulerGamma + 1.0));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(1.0 / 2.0), -AdvancedMath.EulerGamma - 2.0 * Math.Log(2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(1.0 / 3.0), -AdvancedMath.EulerGamma - 3.0 * Math.Log(3.0) / 2.0 - Math.PI / 2.0 / Math.Sqrt(3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(1.0 / 4.0), -AdvancedMath.EulerGamma - 3.0 * Math.Log(2.0) - Math.PI / 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(1.0 / 6.0), -AdvancedMath.EulerGamma - 3.0 * Math.Log(3.0) / 2.0 - 2.0 * Math.Log(2.0) - Math.PI / 2.0 * Math.Sqrt(3.0)));

        }

        [TestMethod]
        public void PsiRecurrence () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E4,24)) {
                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(AdvancedMath.Psi(x), 1.0 / x, AdvancedMath.Psi(x + 1.0)));
            }
        }

        [TestMethod]
        public void PsiReflection () {
            // don't let x get too big or inaccuracy of System.Math trig functions for large x will cause failure
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E2,16)) {
                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(AdvancedMath.Psi(x), Math.PI / Math.Tan(Math.PI * x), AdvancedMath.Psi(1.0 - x)));
            }
        }

        [TestMethod]
        public void PsiDuplication () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 24)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.Psi(2 * x), AdvancedMath.Psi(x) / 2.0 + AdvancedMath.Psi(x + 0.5) / 2.0 + Math.Log(2.0)
                ));
            }
        }


        [TestMethod]
        public void TriGammaSpecialCases () {
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(1, 1.0 / 4.0), Math.PI * Math.PI + 8.0 * AdvancedMath.Catalan));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(1, 1.0 / 2.0), Math.PI * Math.PI / 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Psi(1, 1.0), Math.PI * Math.PI / 6.0));
        }

        [TestMethod]
        public void TriGammaReflection () {
            // don't let x get too big, or the problem of trig functions with large arguments will occur
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 16)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.Psi(1, 1.0 - x) + AdvancedMath.Psi(1, x),
                    MoreMath.Pow( Math.PI / Math.Sin(Math.PI * x), 2)
                ));
            }
        }

        [TestMethod]
        public void TetraGammaReflection () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 16)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.Psi(2, 1.0 - x) - AdvancedMath.Psi(2, x),
                    2.0 * MoreMath.Pow(Math.PI / Math.Sin(Math.PI * x), 3) * Math.Cos(Math.PI * x)
                ));
            }
        }

        [TestMethod]
        public void PolyGammaRiemann () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 16)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.Psi(n, 1.0),
                    -MoreMath.Pow(-1, n) * AdvancedIntegerMath.Factorial(n) * AdvancedMath.RiemannZeta(n + 1)
                ));
            }
        }

        [TestMethod]
        public void PolyGammaDuplication () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 16)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        MoreMath.Pow(2, n + 1) * AdvancedMath.Psi(n, 2.0 * x),
                        AdvancedMath.Psi(n, x) + AdvancedMath.Psi(n, x + 0.5)
                    ));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        MoreMath.Pow(3, n + 1) * AdvancedMath.Psi(n, 3.0 * x),
                        AdvancedMath.Psi(n, x) + AdvancedMath.Psi(n, x + 1.0 / 3.0) + AdvancedMath.Psi(n, x + 2.0 / 3.0)
                    ));
                }
            }
        }

        
        [TestMethod]
        public void PolyGammaIntegral () {
            // don't let these n or x values get too extreme, or the very small result will make the integral return without full precision
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 10, 3)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 10)) {
                    double F = FunctionMath.Integrate(t => {
                        return (Math.Exp(-t * x) * MoreMath.Pow(t, n) / MoreMath.ExpMinusOne(-t));
                    }, Interval.FromEndpoints(0.0, Double.PositiveInfinity));
                    if (n % 2 != 0) F = -F;
                    //Console.WriteLine("{0} {1} {2} {3}", n, x, F, AdvancedMath.Psi(n, x));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(F, AdvancedMath.Psi(n, x)));
                }
            }
        }
        

        [TestMethod]
        public void PolyGammaRecurrence () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 5)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 15)) {
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                        AdvancedMath.Psi(n, x), MoreMath.Pow(-1, n) * AdvancedIntegerMath.Factorial(n) / MoreMath.Pow(x, n+1),
                        AdvancedMath.Psi(n, x + 1.0)
                    ));
                }
            }
        }

        [TestMethod]
        public void DigammaAgreement () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 24)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.Psi(0, x), AdvancedMath.Psi(x)
                ));
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
                    Func<double, double> f = delegate(double t) {
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
            Func<double, double> f = delegate(double t) {
                return (AdvancedMath.InverseErf(t));
            };
            Interval r = Interval.FromEndpoints(0.0, 1.0);
            double I = FunctionMath.Integrate(f, r);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(I, 1.0 / Math.Sqrt(Math.PI)));
        }

        [TestMethod]
        public void ErfIntegralTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2,10.0,3)) {
                Func<double, double> f = delegate(double t) {
                    return (Math.Exp(-t * t));
                };
                Interval r = Interval.FromEndpoints(0.0, x);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Erf(x), 2.0 / Math.Sqrt(Math.PI) * FunctionMath.Integrate(f, r)));
            }
        }

        [TestMethod]
        public void ErfcComplementarity () {
            foreach (double x in arguments) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Erf(x) + AdvancedMath.Erfc(x), 1.0));
            }
        }

        [TestMethod]
        public void ErfcInequality () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-5,1.0E5,32)) {
                double erfc = AdvancedMath.Erfc(x);
                double factor = 2.0 / Math.Sqrt(Math.PI) * Math.Exp(-x * x);
                double lower = factor / (x + Math.Sqrt(x * x + 2.0));
                double upper = factor / (x + Math.Sqrt(x * x + 4.0 / Math.PI));
                Assert.IsTrue((lower <= erfc) && (erfc <= upper));
            }

        }

        [TestMethod]
        public void ErfcIntegral () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 10.0, 4)) {
                Func<double, double> f = delegate(double t) {
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
        public void IntegralE1Inequality () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E3,16)) {
                if (x >= Math.Log(Double.MaxValue / 10.0)) continue; // keep Exp(x) from overflowing
                double lower = Math.Log(1.0 + 2.0 / x) / 2.0;
                double upper = Math.Log(1.0 + 1.0 / x);
                double value = Math.Exp(x) * AdvancedMath.IntegralE(1, x);
                Console.WriteLine("{0}: {1} {2} {3}", x, lower, value, upper);
                Assert.IsTrue((lower < value) && (value < upper));
            }
        }

        [TestMethod]
        public void IntegralEInequality () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1,100,8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E3,16)) {
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
        public void IntegralERecurrence () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 24)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        n * AdvancedMath.IntegralE(n + 1, x) + x * AdvancedMath.IntegralE(n, x),
                        Math.Exp(-x)
                    ));
                }
            }
        }

        [TestMethod]
        public void IntegralEIntegral () {
            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 10.0, 8)) {
                    Func<double, double> f = delegate(double t) {
                        return (Math.Exp(-x * t) / MoreMath.Pow(t, n));
                    };
                    Interval r = Interval.FromEndpoints(1.0, Double.PositiveInfinity);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.IntegralE(n, x), FunctionMath.Integrate(f, r)));
                }
            }
        }

        [TestMethod]
        public void IntegralEIncompleteGamma () {

            foreach (int n in TestUtilities.GenerateIntegerValues(1, 100, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 16)) {

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.IntegralE(-n, x),
                        AdvancedMath.Gamma(n + 1, x) / MoreMath.Pow(x, n + 1)
                    ));

                }
            }

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
            Func<double, double> f = delegate(double t) {
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
            Func<double, double> f = delegate(double t) {
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
            Func<double, double> f = delegate(double t) {
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
            Func<double, double> f = delegate(double t) {
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
        public void ModifiedBesselTower () {

            // e^z = I_0 + 2 I_1 + 2 I_2 + 2 I_3 + ...
            // cosh z = I_0 2 + 2 I_2 + 2 I_4 + 2 I_6 + ...
            // sinh z = 2 I_1 + 2 I_3 + 2 I_5

            // it takes about n ~ x terms to converge, so don't let x get too big

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 16)) {

                double I = AdvancedMath.ModifiedBesselI(0, x);

                double exp = I;
                double cosh = I;
                double sinh = 0.0;

                for (int n = 1; n < 100; n++) {

                    double exp_old = exp;
                    I = 2.0 * AdvancedMath.ModifiedBesselI(n, x);
                    exp += I;
                    if (n % 2 == 0) {
                        cosh += I;
                    } else {
                        sinh += I;
                    }
                    if (exp == exp_old) {
                        //Console.WriteLine("{0} {1}", x, n);
                        break;
                    }
                }

                Assert.IsTrue(TestUtilities.IsNearlyEqual(exp, Math.Exp(x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(cosh, Math.Cosh(x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(sinh, Math.Sinh(x)));

            }

        }

        [TestMethod]
        public void ModifiedBesselIntegralTest () {

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 5)) {
                //Console.WriteLine("x={0}", x);

                Func<double, double> f = delegate (double t) {
                    return (Math.Exp(x * Math.Cos(t)));
                };
                Interval r = Interval.FromEndpoints(0.0, Math.PI);
                double I = FunctionMath.Integrate(f, r);

                //Console.WriteLine("{0} ?= {1}", I / Math.PI, AdvancedMath.ModifiedBesselI(0.0, x));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(I / Math.PI, AdvancedMath.ModifiedBesselI(0.0, x)));

            }

        }

        [TestMethod]
        public void ModifiedBesselAgreement () {

            foreach (double nu in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 16)) {

                    Console.WriteLine("{0} {1}", nu, x);

                    SolutionPair s = AdvancedMath.ModifiedBessel(nu, x);
                    double I = AdvancedMath.ModifiedBesselI(nu, x);
                    double K = AdvancedMath.ModifiedBesselK(nu, x);

                    Console.WriteLine("  {0} {1}", s.FirstSolutionValue, I);
                    Console.WriteLine("  {0} {1}", s.SecondSolutionValue, K);

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(s.FirstSolutionValue, AdvancedMath.ModifiedBesselI(nu, x)));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(s.SecondSolutionValue, AdvancedMath.ModifiedBesselK(nu, x)));

                }
            }

        }

        [TestMethod]
        public void ModifiedBesselWronskian () {

            foreach (double nu in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 8)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 16)) {

                    SolutionPair s = AdvancedMath.ModifiedBessel(nu, x);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        s.SecondSolutionValue * s.FirstSolutionDerivative - s.FirstSolutionValue * s.SecondSolutionDerivative,
                        1.0 / x
                    ));
                    // there is no cancelation, because I, I', K > 0 and K' < 0

                }
            }

        }

        [TestMethod]
        public void BesselModifiedBesselRelationship () {

            // This is a very interesting test because it relates (normal Bessel) J and (modified Bessel) I
            //   I_{\nu}(x) = \sum_{k=0}^{\infty} J_{\nu + k}(x) \frac{x^k}{k!}
            // We want x not too big, so that \frac{x^k}{k!} converges, and x <~ \nu, so that J_{\nu+k} decreases rapidly

            foreach (double nu in TestUtilities.GenerateRealValues(0.1, 10.0, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(0.1, nu, 4)) {

                    double I = 0.0;
                    for (int k = 0; k < 100; k++) {

                        double I_old = I;
                        I += AdvancedMath.BesselJ(nu + k, x) * MoreMath.Pow(x, k) / AdvancedIntegerMath.Factorial(k);
                        if (I == I_old) {
                            Console.WriteLine("{0} {1} {2}", nu, x, k);
                            break;
                        }

                    }

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(I, AdvancedMath.ModifiedBesselI(nu, x)));

                }
            }

        }

        [TestMethod]
        public void AiryIntegral () {

            Func<double, double> f = delegate(double t) {
                return (AdvancedMath.AiryAi(t));
            };
            Interval r = Interval.FromEndpoints(0.0, Double.PositiveInfinity);
            double I = FunctionMath.Integrate(f, r);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(I, 1.0 / 3.0));

        }

        [TestMethod]
        public void AiryZeroArgument () {

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.AiryAi(0.0), 1.0 / Math.Pow(3.0, 2.0 / 3.0) / AdvancedMath.Gamma(2.0 / 3.0)
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.AiryBi(0.0), 1.0 / Math.Pow(3.0, 1.0 / 6.0) / AdvancedMath.Gamma(2.0 / 3.0)
            ));
        }

        [TestMethod]
        public void AiryIntegrals () {
            
            //double I3 = FunctionMath.Integrate(t => MoreMath.Pow(AdvancedMath.AiryAi(t), 3), Interval.FromEndpoints(Double.NegativeInfinity, Double.PositiveInfinity));
            //Assert.IsTrue(TestUtilities.IsNearlyEqual(I3, MoreMath.Pow(AdvancedMath.Gamma(1.0 / 3.0) / Math.PI, 2) / 4.0));

            double I4 = FunctionMath.Integrate(t => MoreMath.Pow(AdvancedMath.AiryAi(t), 4), Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(I4, Math.Log(3.0) / (24.0 * Math.PI * Math.PI)));
            
        }

        [TestMethod]
        public void AiryBairyIntegral () {

            // for small x: the integral will not converse for x <~ 5 because of the 1/sqrt(t) divergence near 0
            // for large x: Bi explodes exponentially
            // but in between, it is a decent test of Ai and Bi together
            foreach (double x in TestUtilities.GenerateRealValues(5.0, 50.0, 8)) {

                double I = FunctionMath.Integrate(
                    t => { return (Math.Exp(x * t - t * t * t / 12.0) / Math.Sqrt(t)); },
                    Interval.FromEndpoints(0.0, Double.PositiveInfinity)
                );

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    MoreMath.Pow(AdvancedMath.AiryAi(x), 2) + MoreMath.Pow(AdvancedMath.AiryBi(x), 2), I / Math.Pow(Math.PI, 3.0 / 2.0)
                ));

            }

        }

        [TestMethod]
        public void GoldenRatioTest () {

            Assert.IsTrue(TestUtilities.IsNearlyEqual((1.0 + AdvancedMath.GoldenRatio) / AdvancedMath.GoldenRatio, AdvancedMath.GoldenRatio));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(1.0 / AdvancedMath.GoldenRatio, AdvancedMath.GoldenRatio - 1.0));
        }

        [TestMethod]
        public void CatalanIntegralTest () {

            Func<double, double> f1 = delegate(double t) {
                return (t / Math.Sin(t) / Math.Cos(t));
            };
            Interval r1 = Interval.FromEndpoints(0.0, Math.PI / 4.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(FunctionMath.Integrate(f1, r1), AdvancedMath.Catalan));

            Func<double, double> f2 = delegate(double t) {
                return (-Math.Log(t) / (1 + t * t));
            };
            Interval r2 = Interval.FromEndpoints(0.0, 1.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(FunctionMath.Integrate(f2, r2), AdvancedMath.Catalan));

        }

        [TestMethod]
        public void DiLogSpecialCases () {

            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.DiLog(-1.0), -Math.PI * Math.PI / 12.0));
            Assert.IsTrue(AdvancedMath.DiLog(0.0) == 0.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.DiLog(0.5), Math.PI * Math.PI / 12.0 - Math.Log(2.0) * Math.Log(2.0) / 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.DiLog(1.0 / AdvancedMath.GoldenRatio), Math.PI * Math.PI / 10.0 - Math.Log(AdvancedMath.GoldenRatio) * Math.Log(AdvancedMath.GoldenRatio)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.DiLog(1.0), Math.PI * Math.PI / 6.0));

        }

        [TestMethod]
        public void DiLogDuplication () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0, 10)) {
                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                    AdvancedMath.DiLog(x), AdvancedMath.DiLog(-x), AdvancedMath.DiLog(x * x) / 2.0
                ));
            }
        }

        [TestMethod]
        public void DiLogBailyIdentity () {

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                36.0 * AdvancedMath.DiLog(1.0 / 2.0) - 36.0 * AdvancedMath.DiLog(1.0 / 4.0) -
                12.0 * AdvancedMath.DiLog(1.0 / 8.0) + 6.0 * AdvancedMath.DiLog(1.0 / 64.0),
                Math.PI * Math.PI
            ));

        }

        [TestMethod]
        public void DiLogExpIntegral () {

            foreach (double x in TestUtilities.GenerateRealValues(0.1,10.0,5)) {

                Func<double, double> f = delegate (double t) {
                    return (t / (Math.Exp(t) + x));
                };
                Interval r = Interval.FromEndpoints(0.0, Double.PositiveInfinity);
                double I = FunctionMath.Integrate(f, r);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.DiLog(-x), -x * I));

            }

        }

        [TestMethod]
        public void DiLogLogIntegral () {

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 5)) {

                Func<double, double> f = delegate(double t) {
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
            // infinity at 1

        }

        [TestMethod]
        public void EllipticKInequality () {
            // DLMF 19.9.1
            foreach (double k in TestUtilities.GenerateUniformRealValues(0.0, 1.0, 4)) {
                double S = AdvancedMath.EllipticK(k) + Math.Log(Math.Sqrt(1.0 - k * k));
                Assert.IsTrue((Math.Log(4.0) <= S) && (S <= Math.PI / 2.0));
            }
        }

        [TestMethod]
        public void EllipticKIntegration () {

            // The defining integral for K

            Interval i = Interval.FromEndpoints(0.0, Math.PI / 2.0);

            foreach (double k in TestUtilities.GenerateRealValues(0.01, 1.0, 8)) {

                Func<double, double> f = delegate(double t) {
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

            Func<double, double> f1 = delegate(double k) {
                return (AdvancedMath.EllipticK(k));
            };
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                FunctionMath.Integrate(f1, i), 2.0 * AdvancedMath.Catalan
            ));

            Func<double, double> f2 = delegate(double k) {
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

            foreach (double k in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 8)) {

                Func<double, double> f = delegate(double t) {
                    double z = k * Math.Sin(t);
                    return (1.0 / Math.Sqrt(1.0 - z * z));
                };

                foreach (double phi in TestUtilities.GenerateUniformRealValues(0.0, Math.PI/2.0, 8)) {

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
        public void CarlsonFSpecialCases () {

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 8)) {

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.CarlsonF(x, x, x),
                    1.0 / Math.Sqrt(x)
                ));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.CarlsonF(0, x, x),
                    Math.PI / 2.0 / Math.Sqrt(x)
                ));

            }

        }

        [TestMethod]
        public void CarlsonFInequality () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-3, 1.0E3, 8)) {
                foreach (double y in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 8)) {
                    foreach (double z in TestUtilities.GenerateRealValues(1.0E-1, 1.0E5, 8)) {

                        double min = 3.0 / (Math.Sqrt(x) + Math.Sqrt(y) + Math.Sqrt(z));
                        double R_F = AdvancedMath.CarlsonF(x, y, z);
                        double max = 1.0 / Math.Pow(x * y * z, 1.0 / 6.0);

                        Assert.IsTrue(min <= R_F && R_F <= max);

                    }
                }
            }

        }

        [TestMethod]
        public void CarlsonDSymmetrizedSum () {
            // Carlson 1994 equation (54)
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-3, 1.0E3, 8)) {
                foreach (double y in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 8)) {
                    foreach (double z in TestUtilities.GenerateRealValues(1.0E-1, 1.0E5, 8)) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(
                            AdvancedMath.CarlsonD(x, y, z) + AdvancedMath.CarlsonD(z, x, y) + AdvancedMath.CarlsonD(y, z, x),
                            3.0 / Math.Sqrt(x * y * z)
                        ));
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(
                            z * AdvancedMath.CarlsonD(x, y, z) + y * AdvancedMath.CarlsonD(z, x, y) + x * AdvancedMath.CarlsonD(y, z, x),
                            3.0 * AdvancedMath.CarlsonF(x, y, z)
                        ));
                    }
                }
            }

        }

        [TestMethod]
        public void CarslonDSpecialCases () {

            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 8)) {

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.CarlsonD(x, x, x), Math.Pow(x, -3.0 / 2.0)
                ));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.CarlsonD(0, x, x), 3.0 * Math.PI / 4.0 * Math.Pow(x, -3.0 / 2.0)
                ));

            }

        }

        [TestMethod]
        public void CarlsonLegendreRelation () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 16)) {

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.CarlsonF(0.0, x + 1.0, 1.0) * AdvancedMath.CarlsonD(0.0, x + 1.0, x) +
                    AdvancedMath.CarlsonD(0.0, x + 1.0, 1.0) * AdvancedMath.CarlsonF(0.0, x + 1.0, x),
                    3.0 / 2.0 * Math.PI / x
                ));

            }

        }

        [TestMethod]
        public void CarlsonLemniscaticValues () {
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.CarlsonF(0.0, 1.0, 2.0), MoreMath.Pow(AdvancedMath.Gamma(1.0 / 4.0), 2) / Math.Sqrt(2.0 * Math.PI) / 4.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.CarlsonD(0.0, 2.0, 1.0), 3.0 * MoreMath.Pow(AdvancedMath.Gamma(3.0 / 4.0), 2) / Math.Sqrt(2.0 * Math.PI)));
        }

        [TestMethod]
        public void EllipticESPecialCases () {
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticE(0.0), Math.PI / 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticE(1.0), 1.0));
            
            double g2 = Math.Pow(AdvancedMath.Gamma(1.0 / 4.0), 2.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.EllipticE(1.0 / Math.Sqrt(2.0)),
                Math.Pow(Math.PI, 3.0 / 2.0) / g2 + g2 / 8.0 / Math.Sqrt(Math.PI)
            ));
             
        }

        [TestMethod]
        public void EllipticLegendreRelation () {
            foreach (double k in TestUtilities.GenerateRealValues(0.01, 1.0, 8)) {
                double kp = Math.Sqrt(1.0 - k * k);
                double E = AdvancedMath.EllipticE(k);
                double EP = AdvancedMath.EllipticE(kp);
                double K = AdvancedMath.EllipticK(k);
                double KP = AdvancedMath.EllipticK(kp);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(K * EP + KP * E, Math.PI / 2.0 + K * KP));
            }
        }

        [TestMethod]
        public void EllipticEIntegration () {

            foreach (double k in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 8)) {

                // define the integrand
                Func<double, double> f = delegate(double t) {
                    double z = k * Math.Sin(t);
                    return (Math.Sqrt(1.0 - z * z));
                };

                // test complete integral
                double C = FunctionMath.Integrate(f, Interval.FromEndpoints(0.0, Math.PI / 2.0));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticE(k), C));

                // test incomplete integral
                foreach (double phi in TestUtilities.GenerateUniformRealValues(0.0, Math.PI / 2.0, 8)) {
                    double I = FunctionMath.Integrate(f, Interval.FromEndpoints(0.0, phi));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.EllipticE(phi, k), I));
                }

            }

        }

        [TestMethod]
        public void EllipticLandenTransform () {

            foreach (double k in TestUtilities.GenerateUniformRealValues(0.0, 1.0, 8)) {

                double kp = Math.Sqrt(1.0 - k * k);
                double k1 = (1.0 - kp) / (1.0 + kp);

                // For complete 1st kind
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.EllipticK(k),
                    (1.0 + k1) * AdvancedMath.EllipticK(k1)
                ));

                // For complete 2nd kind
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.EllipticE(k) + kp * AdvancedMath.EllipticK(k),
                    (1.0 + kp) * AdvancedMath.EllipticE(k1)
                ));

                /*
                foreach (double t in TestUtilities.GenerateUniformRealValues(0.0, Math.PI / 2.0, 4)) {

                    double s = Math.Sin(t);
                    double t1 = Math.Asin((1.0 + kp) * s * Math.Sqrt(1.0 - s * s) / Math.Sqrt(1.0 - k * k * s * s));

                    // Since t1 can be > pi/2, we need to handle a larger range of angles before we can test the incomplete cases
 

                }

                */
            }

        }

        [TestMethod]
        public void EllipticKCarlsonFRelationship () {

            foreach (double k in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 8)) {

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.EllipticK(k),
                    AdvancedMath.CarlsonF(0.0, 1.0 - k * k, 1.0)
                ));

            }

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
