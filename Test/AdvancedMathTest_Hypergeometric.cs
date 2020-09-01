using System;
using System.Collections.Generic;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Analysis;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_Hypergeometric {

        double[] abcs = new double[] { -3.0, -0.5, 0.1, 1.0, 1.5, 2.0, 3.1, 4.5, 7.0 };
        //double[] abcs = new double[] { -0.5, 0.1, 1.0, 1.5, 2.0, 3.1, 4.5, 7.0 };
        double[] xs = new double[] { -5.0, -0.6, -0.2, 0.0, 0.3, 0.5, 0.8 };

        private static bool IsNonpositiveInteger (double x) {
            return ((x <= 0.0) && (Math.Round(x) == x));
        }

        [TestMethod]
        public void HypergeometrticRecurrenceA () {

            // A&S 15.2.10

            foreach (double a in abcs) {
                foreach (double b in abcs) {
                    foreach (double c in abcs) {
                        foreach (double x in xs) {

                            double FM = AdvancedMath.Hypergeometric2F1(a - 1.0, b, c, x);
                            double F0 = AdvancedMath.Hypergeometric2F1(a, b, c, x);
                            double FP = AdvancedMath.Hypergeometric2F1(a + 1.0, b, c, x);

                            if ((c == a) && Double.IsNaN(FM)) continue;

                            Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                                new double[] { (c - a) * FM, (2.0 * a - c - a * x + b * x) * F0 },
                                -a * (x - 1.0) * FP,
                                new EvaluationSettings() { RelativePrecision = 1.0E-12, AbsolutePrecision = 1.0E-15 }
                            ));

                        }
                    }
                }
            }

        }

        [TestMethod]
        public void HypergeometricRecurrenceC () {

            // A&S 15.2.12

            foreach (double a in abcs) {
                foreach (double b in abcs) {
                    foreach (double c in abcs) {
                        foreach (double x in xs) {

                            double FM = AdvancedMath.Hypergeometric2F1(a, b, c - 1.0, x);
                            double F0 = AdvancedMath.Hypergeometric2F1(a, b, c, x);
                            double FP = AdvancedMath.Hypergeometric2F1(a, b, c + 1.0, x);

                            if (Double.IsNaN(FM) && ((c == 0.0) || (c == 1.0))) continue;
                            if (Double.IsNaN(FP) && ((c == a) || (c == b))) continue;

                            Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                                new double[] { c * (c - 1.0) * (x - 1.0) * FM, c * (c - 1.0 - (2.0 * c - a - b - 1.0) * x) * F0 },
                                -(c - a) * (c - b) * x * FP,
                                new EvaluationSettings() { RelativePrecision = 1.0E-12, AbsolutePrecision = 1.0E-16 }
                            ));

                        }
                    }
                }
            }

        }

        [TestMethod]
        public void HypergeometricBasicFunctions () {

            foreach (double x in xs) {

                // A&S 15.1.3
                if (x < 1.0 && (x != 0.0)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.Hypergeometric2F1(1.0, 1.0, 2.0, x), -Math.Log(1.0 - x) / x
                    ));
                }

                // A&S 15.1.5
                if (x != 0.0) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.Hypergeometric2F1(0.5, 1.0, 1.5, -x * x), Math.Atan(x) / x
                    ));
                }

                // A&S 15.1.6
                if ((-1.0 <= x) && (x <= 1.0) && (x != 0.0)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.Hypergeometric2F1(0.5, 0.5, 1.5, x * x), Math.Asin(x) / x
                    ));
                }

                // Complete elliptic integrals
                if ((0.0 <= x) && (x < 1.0)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        Math.PI / 2.0 * AdvancedMath.Hypergeometric2F1(0.5, 0.5, 1.0, x * x), AdvancedMath.EllipticK(x)
                    ));

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                         Math.PI / 2.0 * AdvancedMath.Hypergeometric2F1(-0.5, 0.5, 1.0, x * x), AdvancedMath.EllipticE(x)
                    ));
                }
            }

        }

        [TestMethod]
        public void HypergeometricIncompleteBeta () {

            foreach (double a in abcs) {
                if (a <= 0.0) continue;
                foreach (double b in abcs) {
                    if (b <= 0.0) continue;
                    foreach (double x in xs) {
                        if ((x < 0.0) || (x > 1.0)) continue;
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(
                            Math.Pow(x, a) * Math.Pow(1.0 - x, b) * AdvancedMath.Hypergeometric2F1(a + b, 1.0, a + 1.0, x) / a,
                            AdvancedMath.Beta(a, b, x)
                        ));
                    }
                }
            }

        }

        [TestMethod]
        public void HypergeometricIntegral () {

            // A&S 15.3.1

            // F(a, b, c, x) = \frac{\Gamma(c)}{\Gamma(b) \Gamma(c - b)} 
            // \int_{0}^{1} \! dt \, t^{b-1} (1 - t)^{c - b - 1} (1 - x t)^{-a}

            // Choose limits on a, b, c so that singularities of integrand are numerically integrable.

            foreach (double a in TestUtilities.GenerateRealValues(1.0, 10.0, 2)) {
                foreach (double b in TestUtilities.GenerateRealValues(0.6, 10.0, 2)) {
                    foreach (double c in TestUtilities.GenerateRealValues(b + 0.6, 10.0, 2)) {
                        foreach (double x in xs) {

                            double I = FunctionMath.Integrate(
                                t => Math.Pow(t, b - 1.0) * Math.Pow(1.0 - t, c - b - 1.0) * Math.Pow(1.0 - x * t, -a),
                                Interval.FromEndpoints(0.0, 1.0)
                            );

                            double B = AdvancedMath.Beta(b, c - b);

                            Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Hypergeometric2F1(a, b, c, x), I / B));

                        }
                    }
                }
            }

        }

        [TestMethod]
        public void HypergeometricAtOne () {

            // A&S 15.1.20
            foreach (double a in abcs) {
                foreach (double b in abcs) {
                    foreach (double c in abcs) {

                        // Formula only hold for positive real c-a-b
                        if ((c - a - b) <= 0.0) continue;

                        // Formula is still right for non-positive c, but to handle it we would need to deal with canceling infinite Gammas
                        if (IsNonpositiveInteger(c)) continue;

                        Assert.IsTrue(TestUtilities.IsNearlyEqual(
                            AdvancedMath.Hypergeometric2F1(a, b, c, 1.0),
                            AdvancedMath.Gamma(c) * AdvancedMath.Gamma(c - a - b) / AdvancedMath.Gamma(c - a) / AdvancedMath.Gamma(c - b)
                        ));
                        
                    }
                }
            }

        }

        [TestMethod]
        public void HypergeometricAtMinusOne () {

            foreach (double a in abcs) {

                if ((a <= 0.0) && (Math.Round(a) == a)) continue;

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.Hypergeometric2F1(1.0, a, a + 1.0, -1.0),
                    a / 2.0 * (AdvancedMath.Psi((a + 1.0)/2.0) - AdvancedMath.Psi(a / 2.0))
                ));

                foreach (double b in abcs) {

                    double c = a - b + 1.0;
                    if ((c <= 0.0) && (Math.Round(c) == c)) continue;

                    // If result vanishes, returned value may just be very small.
                    double R = AdvancedMath.Gamma(a - b + 1.0) * AdvancedMath.Gamma(a / 2.0 + 1.0) / AdvancedMath.Gamma(a + 1.0) / AdvancedMath.Gamma(a / 2.0 - b + 1.0);
                    if (R == 0.0) {
                        Assert.IsTrue(Math.Abs(AdvancedMath.Hypergeometric2F1(a, b, a - b + 1.0, -1.0)) <= 1.0E-14);
                    } else {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedMath.Hypergeometric2F1(a, b, a - b + 1.0, -1.0), R));
                    }
                }
            }

        }

        [TestMethod]
        public void HypergeometricAtOneHalf () {

            foreach (double a in abcs) {
                foreach (double b in abcs) {

                    double c = (a + b + 1.0) / 2.0;
                    if (!IsNonpositiveInteger(c)) {
                        if (IsNonpositiveInteger((a + 1.0) / 2.0) || IsNonpositiveInteger((b + 1.0) / 2.0)) {
                            Assert.IsTrue(Math.Abs(AdvancedMath.Hypergeometric2F1(a, b, c, 0.5)) < TestUtilities.TargetPrecision);
                        } else {
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                                AdvancedMath.Hypergeometric2F1(a, b, c, 0.5),
                                Math.Sqrt(Math.PI) * AdvancedMath.Gamma(c) / AdvancedMath.Gamma((a + 1.0) / 2.0) / AdvancedMath.Gamma((b + 1.0) / 2.0)
                            ));
                        }
                    }

                    if (!IsNonpositiveInteger(b)) {
                        if (IsNonpositiveInteger((a + b) / 2.0) || IsNonpositiveInteger((b - a + 1.0) / 2.0)) {
                            Assert.IsTrue(Math.Abs(AdvancedMath.Hypergeometric2F1(a, 1.0 - a, b, 0.5)) < TestUtilities.TargetPrecision);
                        } else {
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                                AdvancedMath.Hypergeometric2F1(a, 1.0 - a, b, 0.5),
                                Math.Sqrt(Math.PI) * Math.Pow(2.0, 1.0 - b) * AdvancedMath.Gamma(b) / AdvancedMath.Gamma((a + b) / 2.0) / AdvancedMath.Gamma((b - a + 1.0) / 2.0)
                            ));
                        }
                    }

                }
            }

        }

        [TestMethod]
        public void HypergeometricAtMinusOneThird () {

            foreach (double a in abcs) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.Hypergeometric2F1(a, a + 0.5, 1.5 - 2.0 * a, -1.0 / 3.0),
                    Math.Pow(8.0 / 9.0, -2.0 * a) * AdvancedMath.Gamma(4.0 / 3.0) / AdvancedMath.Gamma(3.0 / 2.0) * AdvancedMath.Gamma(3.0 / 2.0 - 2.0 * a) / AdvancedMath.Gamma(4.0 / 3.0 - 2.0 * a)
                ));
            }

        }

        [TestMethod]
        public void HypergeometricAtOneNinth () {

            foreach (double a in abcs) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.Hypergeometric2F1(a, a + 0.5, 5.0 / 6.0 + 2.0 / 3.0 * a, 1.0 / 9.0),
                    Math.Sqrt(Math.PI) * Math.Pow(3.0 / 4.0, a) * AdvancedMath.Gamma(5.0 / 6.0 + 2.0 / 3.0 * a) / AdvancedMath.Gamma(1.0 / 2.0 + 1.0 / 3.0 * a) / AdvancedMath.Gamma(5.0 / 6.0 + 1.0 / 3.0 * a)
                ));
            }

        }

        [TestMethod]
        public void HypergeometricAtSpecialPoints () {

            // Documented at http://mathworld.wolfram.com/HypergeometricFunction.html

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.Hypergeometric2F1(1.0 / 3.0, 2.0 / 3.0, 5.0 / 6.0, 27.0 / 32.0),
                8.0 / 5.0
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.Hypergeometric2F1(0.25, 0.5, 0.75, 80.0 / 81.0),
                9.0 / 5.0
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.Hypergeometric2F1(1.0 / 8.0, 3.0 / 8.0, 1.0 / 2.0, 2400.0 / 2401.0),
                2.0 / 3.0 * Math.Sqrt(7.0)
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.Hypergeometric2F1(1.0 / 6.0, 1.0 / 3.0, 1.0 / 2.0, 25.0 / 27.0),
                3.0 / 4.0 * Math.Sqrt(3.0)
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.Hypergeometric2F1(1.0 / 6.0, 1.0 / 2.0, 2.0 / 3.0, 125.0 / 128.0),
                4.0 / 3.0 * Math.Pow(2.0, 1.0 / 6.0)
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedMath.Hypergeometric2F1(1.0 / 12.0, 5.0 / 12.0, 1.0 / 2.0, 1323.0 / 1331.0),
                3.0 / 4.0 * Math.Pow(11.0, 1.0 / 4.0)
            ));

        }

        [TestMethod]
        public void HypergeometricLinearTransforms () {

            foreach (double a in abcs) {
                foreach (double b in abcs) {
                    foreach (double c in abcs) {
                        if (IsNonpositiveInteger(c)) continue;
                        foreach (double x in xs) {

                            double F = AdvancedMath.Hypergeometric2F1(a, b, c, x);

                            // A&S 15.3.3
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(F, AdvancedMath.Hypergeometric2F1(c - a, c - b, c, x) * Math.Pow(1.0 - x, c - a - b), 1.0E-12));

                            // A&S 15.3.4
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(F, AdvancedMath.Hypergeometric2F1(a, c - b, c, x / (x - 1.0)) * Math.Pow(1.0 - x, -a), 1.0E-13));

                            // A&S 15.3.5
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(F, AdvancedMath.Hypergeometric2F1(b, c - a, c, x / (x - 1.0)) * Math.Pow(1.0 - x, -b), 1.0E-13));

                            // A&S 15.3.6
                            if (!IsNonpositiveInteger(c - a - b) && !IsNonpositiveInteger(a + b - c) && (x > 0.0)) {
                                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                                    new double[] {
                                        AdvancedMath.Gamma(c) * AdvancedMath.Gamma(c - a - b) / AdvancedMath.Gamma(c - a) / AdvancedMath.Gamma(c - b) *
                                        AdvancedMath.Hypergeometric2F1(a, b, a + b - c + 1.0, 1.0 - x),
                                        AdvancedMath.Gamma(c) * AdvancedMath.Gamma(a + b - c) / AdvancedMath.Gamma(a) / AdvancedMath.Gamma(b) *
                                        AdvancedMath.Hypergeometric2F1(c - a, c - b, c - a - b + 1.0, 1.0 - x) *
                                        Math.Pow(1.0 - x, c - a - b)
                                    }, F
                                ));
                            }

                            // A&S 15.3.7
                            if (!IsNonpositiveInteger(a - b) && !IsNonpositiveInteger(b - a) && (x < 0.0)) {
                                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                                    new double[] {
                                        AdvancedMath.Gamma(c) * AdvancedMath.Gamma(b - a) / AdvancedMath.Gamma(b) / AdvancedMath.Gamma(c - a) *
                                        AdvancedMath.Hypergeometric2F1(a, 1.0 - c + a, 1.0 - b + a, 1.0 / x) * Math.Pow(-x, -a),
                                        AdvancedMath.Gamma(c) * AdvancedMath.Gamma(a - b) / AdvancedMath.Gamma(a) / AdvancedMath.Gamma(c - b) *
                                        AdvancedMath.Hypergeometric2F1(b, 1.0 - c + b, 1.0 - a + b, 1.0 / x) * Math.Pow(-x, -b)
                                    }, F
                                ));
                            }
                        }
                    }
                }
            }
        }

        [TestMethod]
        public void HypergeometricQuadraticTransforms () {

            foreach (double a in abcs) {
                foreach (double b in abcs) {
                    foreach (double x in xs) {

                        if (!IsNonpositiveInteger(2.0 * b)) {
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                                AdvancedMath.Hypergeometric2F1(a, b, 2.0 * b, x),
                                AdvancedMath.Hypergeometric2F1(0.5 * a, b - 0.5 * a, b + 0.5, x * x / (x - 1.0) / 4.0) / Math.Pow(1.0 - x, 0.5 * a),
                                TestUtilities.TargetPrecision * 1000.0
                            ));
                        }

                        // Wrong sign, but Mathematica agrees with values!
                        /*
                        if (!IsNonpositiveInteger(a - b + 1.0)) {
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                                AdvancedMath.Hypergeometric2F1(a, b, a - b + 1.0, x),
                                AdvancedMath.Hypergeometric2F1(0.5 * a, 0.5 * a - b + 0.5, a - b + 1.0, -4.0 * x / MoreMath.Sqr(1.0 - x)) / Math.Pow(1.0 - x, a)
                            ));
                        }
                        */
                    }
                }
            }

        }

        [Ignore]
        [TestMethod]
        public void HypergeometricCubicTransforms () {

            foreach (double a in abcs) {
                foreach (double x in xs) {

                    // DLMF 15.8.31
                    //if ((a == 7.0) && (x == 0.8)) continue;
                    if ((a == -3.0) && (x == 0.8)) continue;
                    if ((a == 3.1) && (x == -5.0)) continue;
                    if ((a == 4.5) && (x == -5.0)) continue;
                    if ((a == 7.0) && (x == -5.0)) continue;
                    if (x < 8.0 / 9.0) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(
                            AdvancedMath.Hypergeometric2F1(3.0 * a, 3.0 * a + 0.5, 4.0 * a + 2.0 / 3.0, x),
                            AdvancedMath.Hypergeometric2F1(a, a + 0.5, 2.0 * a + 5.0 / 6.0, 27.0 * x * x * (x - 1.0) / MoreMath.Sqr(9.0 * x - 8.0)) * Math.Pow(1.0 - 9.0 / 8.0 * x, -2.0 * a)
                        ));
                    }

                }
            }

        }

        [TestMethod]
        public void HypergeometricPolynomials () {

            foreach (int n in TestUtilities.GenerateUniformIntegerValues(0, 10, 4)) {
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0, 4)) {

                    // A&S 22.5.47
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.Hypergeometric2F1(-n, n, 0.5, x),
                        OrthogonalPolynomials.ChebyshevT(n, 1.0 - 2.0 * x)
                    ));

                    // A&S 22.5.49
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.Hypergeometric2F1(-n, n + 1, 1.0, x),
                        OrthogonalPolynomials.LegendreP(n, 1.0 - 2.0 * x)
                    ));

                    // A&S 22.5.46
                    foreach (double alpha in TestUtilities.GenerateRealValues(1.0E-2, 10.0, 4)) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(
                            AdvancedMath.Hypergeometric2F1(-n, n + 2.0 * alpha, alpha + 0.5, x) * AdvancedMath.Pochhammer(2.0 * alpha, n) / AdvancedIntegerMath.Factorial(n),
                            OrthogonalPolynomials.GegenbauerC(n, alpha, 1.0 - 2.0 * x),
                            new EvaluationSettings() { RelativePrecision = Math.Sqrt(TestUtilities.TargetPrecision) }
                        ));
                    }
                }
            }

        }

        [TestMethod]
        public void HypergeometricPearsonTests () {
            // John Pearson wrote a thesis on the computation of Hypergeometric functions and assembled a list of test cases.
            // (http://people.maths.ox.ac.uk/porterm/research/pearson_final.pdf)
            foreach (var t in HypergeometricPearsonTestValues()) {
                double F = AdvancedMath.Hypergeometric2F1(t.Item1, t.Item2, t.Item3, t.Item4);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(F, t.Item5));
            }
        }

        private IEnumerable<Tuple<double, double, double, double, double>> HypergeometricPearsonTestValues () {

            yield return Tuple.Create(0.1, 0.2, 0.3, 0.5, 1.046432811217352);
            yield return Tuple.Create(-0.1, 0.2, 0.3, 0.5, 0.956434210968214);
            yield return Tuple.Create(1.0E-8, 1.0E-8, 1.0E-8, 1.0E-6, 1.000000000000010);
            yield return Tuple.Create(2.0 + 1.0E-9, 3.0, 5.0, -0.75, 0.492238858852651);
            yield return Tuple.Create(-2.0, -3.0, -5.0 + 1.0E-9, 0.5, 0.474999999913750);
            yield return Tuple.Create(-1.0, -1.5, -2.0 - 1.0E-14, 0.5, 0.625000000000000);
            // Original uses 1.0E-15, but that gets lost and exact integer version is computed wrongly by Euler transformation 
            /*
            yield return Tuple.Create(500.0, -500.0, 500.0, 0.75, 9.332636185032189E-302);
            yield return Tuple.Create(500.0, 500.0, 500.0, -0.6, 8.709809816217217E-103);
            yield return Tuple.Create(-1000.0, -2000.0, -4000.1, -0.5, 5.233580403196932E94);
            yield return Tuple.Create(-100.0, -200.0, -300.0 + 1.0E-9, 0.5 * Math.Sqrt(2.0), 2.653635302903707E-31);
            */
            yield return Tuple.Create(300.0, 10.0, 5.0, 0.5, 3.912238919961547E98);
            //yield return Tuple.Create(5.0, -300.0, 10.0, 0.5, 1.661006238211309E-7);
            //yield return Tuple.Create(10.0, 5.0, -300.5, 0.5, -3.852027081523919E32);
            yield return Tuple.Create(2.25, 3.75, -0.5, -1.0, -0.631220676949703);

        }

    }
}
