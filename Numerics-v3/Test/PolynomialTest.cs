using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;

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
        public void PolynomialNegation () {

            Polynomial p = Polynomial.FromCoefficients(-1.0, 0.0, 1.0);

            Polynomial np = -p;

            Assert.IsTrue(p.Degree == np.Degree);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(p.Evaluate(2.0), -np.Evaluate(2.0)));

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

        [TestMethod]
        public void PolynomialFitFrom2dArray () {

            Polynomial p0 = Polynomial.FromCoefficients(5.0, -6.0, 7.0, -8.0);
            Console.WriteLine(p0);

            double[,] xy = new double[p0.Degree + 1, 2];
            double[] x = TestUtilities.GenerateUniformRealValues(-1.0, 2.0, p0.Degree + 1);
            for (int i = 0; i < x.Length; i++) {
                xy[i, 0] = x[i];
                xy[i, 1] = p0.Evaluate(x[i]);
            }

            Polynomial p1 = Polynomial.FromPoints(xy);
            Console.WriteLine(p1);

            Assert.IsTrue(p1.Degree == p0.Degree);

            for (int i = 0; i < p0.Degree; i++) {
                Console.WriteLine("{0} {1}", p0.Coefficient(i), p1.Coefficient(i));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(p0.Coefficient(i), p1.Coefficient(i)));
            }

            Assert.IsTrue(TestUtilities.IsNearlyEqual(p0.Evaluate(-1.0), p1.Evaluate(-1.0)));

        }

        [TestMethod]
        public void PolynomialFitFromXY () {

            Polynomial p0 = Polynomial.FromCoefficients(1.0, 0.0, -3.0, 2.0);
            Console.WriteLine(p0);

            XY[] points1 = new XY[p0.Degree + 1];
            for (int x = 0; x < points1.Length; x++) points1[x] = new XY(x, p0.Evaluate(x));

            Polynomial p1 = Polynomial.FromPoints(points1);
            Console.WriteLine(p1);

            for (int i = 0; i < p0.Degree; i++) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(p0.Coefficient(i), p1.Coefficient(i)));
            }
        }

    }
}
