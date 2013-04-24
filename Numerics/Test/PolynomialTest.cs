using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public class PolynomialTest {

        [TestMethod]
        public void PolynomialProperties () {

            Polynomial p = Polynomial.FromCoefficients(1.0, 2.0, 3.0);
            Assert.IsTrue(p.Degree == 2);
            Assert.IsTrue(p.Coefficient(0) == 1.0);
            Assert.IsTrue(p.Coefficient(1) == 2.0);
            Assert.IsTrue(p.Coefficient(2) == 3.0);

        }

        [TestMethod]
        public void PolynomialArithmetic () {

            Polynomial p1 = Polynomial.FromCoefficients(1.0, 2.0, 0.0, 3.0);
            Polynomial p2 = Polynomial.FromCoefficients(3.0, -2.0, 1.0);

            Polynomial p12s = p1 + p2;
            Assert.IsTrue(p12s.Degree == Math.Max(p1.Degree, p2.Degree));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(p12s.Evaluate(1.0), p1.Evaluate(1.0) + p2.Evaluate(1.0)));

            Polynomial p12d = p1 - p2;
            Assert.IsTrue(p12d.Degree == Math.Max(p1.Degree, p2.Degree));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(p12d.Evaluate(1.0), p1.Evaluate(1.0) - p2.Evaluate(1.0)));

            Polynomial p12p = p1 * p2;
            Assert.IsTrue(p12p.Degree == p1.Degree + p2.Degree);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(p12p.Evaluate(1.0), p1.Evaluate(1.0) * p2.Evaluate(1.0)));

        }

        [TestMethod]
        public void PolynomialDivisionIdentity () {

            // (x - 1) is always a factor of (x^n - 1), with quotient equal to the polynomial of degree (n-1) with all unit coefficients
            //   (x^n - 1) / (x - 1) = x^{n-1} + x^{n-2} + \cdots + x^2 + x + 1

            Polynomial v = Polynomial.FromCoefficients(-1, 1);

            foreach (int n in TestUtilities.GenerateIntegerValues(2, 16, 4)) {

                double[] uc = new double[n + 1]; uc[0] = -1; uc[n] = 1;
                Polynomial u = Polynomial.FromCoefficients(uc);

                Polynomial q, r; q = Polynomial.Divide(u, v, out r);

                Assert.IsTrue(q.Degree == n - 1);
                for (int i = 0; i <= q.Degree; i++) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(q.Coefficient(i), 1.0));
                }

                Assert.IsTrue(r.Degree == 0);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(r.Coefficient(0), 0.0));

            }

        }

        [TestMethod]
        public void PolynomialRemainderTheorem () {

            // Polynomial remainder theorem, aka little Bezout theorem, says p(x) / (x - a) has remainder p(a) 


            foreach (int n in TestUtilities.GenerateIntegerValues(2, 16, 4)) {

                Polynomial p = Polynomial.FromCoefficients(TestUtilities.GenerateUniformRealValues(-16.0, 16.0, n));

                foreach (double x in TestUtilities.GenerateUniformRealValues(-16.0, 16.0, 4)) {

                    Polynomial r; Polynomial.Divide(p, Polynomial.FromCoefficients(-x, 1.0), out r);

                    Assert.IsTrue(r.Degree == 0);
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(r.Coefficient(0), p.Evaluate(x)));

                }

            }

        }

        [TestMethod]
        public void PolynomialCalculus () {

            Polynomial p = Polynomial.FromCoefficients(4.0, -3.0, 2.0, -1.0);

            Polynomial i = p.Integrate(3.14159);

            double i1 = i.Evaluate(1.0) - i.Evaluate(0.0);
            double i2 = FunctionMath.Integrate(p.Evaluate, Interval.FromEndpoints(0.0, 1.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(i1, i2));

            Polynomial id = i.Differentiate();


        }

    }
}
