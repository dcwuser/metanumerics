using System;
using Meta.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Test {
    

    [TestClass]
    public class ComplexMathTest {

        private Complex a = new Complex(3.0, 4.0);
        private Complex b = new Complex(1.5, -2.3);

        [TestMethod]
        public void ComplexIParts () {
            Assert.IsTrue(Complex.I.Re == 0.0);
            Assert.IsTrue(Complex.I.Im == 1.0);
        }

        [TestMethod]
        public void ComplexIDefinition () {
            Assert.IsTrue(ComplexMath.Sqrt(-1.0) == Complex.I);
        }

        // Square root

        [TestMethod]
        public void ComplexSqrt () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 16)) {
                Complex sz = ComplexMath.Sqrt(z);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Sqr(sz), z));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Arg(z) / 2.0, ComplexMath.Arg(sz)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(Math.Sqrt(ComplexMath.Abs(z)), ComplexMath.Abs(sz)));
            }
        }

        [TestMethod]
        public void ComplexSqrtSpecialCases () {
            Assert.IsTrue(ComplexMath.Sqrt(0.0) == 0.0);
            Assert.IsTrue(ComplexMath.Sqrt(1.0) == 1.0);
            Assert.IsTrue(ComplexMath.Sqrt(-1.0) == Complex.I);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Sqrt(Complex.I), (ComplexMath.I + 1.0) / Math.Sqrt(2.0)));
        }

        // Trig functions

        [TestMethod]
        public void ComplexPythagorean () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 16)) {
                if (Math.Abs(z.Im) > Math.Log(Double.MaxValue / 10.0)) continue; // sin and cos blow up in imaginary part of argument gets too big 
                Complex sin = ComplexMath.Sin(z);
                Complex cos = ComplexMath.Cos(z);
                Console.WriteLine("z={0} s={1} c={2} s^2+c^2={3}", z, sin, cos, sin*sin + cos*cos);
                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(sin * sin, cos * cos, 1.0));
            }
        }

        [TestMethod]
        public void ComplexDoubleAngle () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 16)) {
                if (Math.Abs(z.Im) > Math.Log(Double.MaxValue / 10.0)) continue;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Sin(2.0 * z), 2.0 * ComplexMath.Sin(z) * ComplexMath.Cos(z)));
                //Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Cos(2.0 * z), ComplexMath.Sqr(ComplexMath.Cos(z)) - ComplexMath.Sqr(ComplexMath.Sin(z))));
            }
        }

        [TestMethod]
        public void ComplexSecTanTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 10)) {
                if (Math.Abs(z.Im) > Math.Log(Double.MaxValue / 10.0)) continue; // sin and cos blow up in imaginary part of argument gets too big 
                Complex sec = 1.0 / ComplexMath.Cos(z);
                Complex tan = ComplexMath.Tan(z);
                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(sec * sec, - tan * tan, 1));
            }
        }

        [TestMethod]
        public void ComplexTrigPeriodicity () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 8)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Sin(z), ComplexMath.Sin(z + 2.0 * Math.PI)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Cos(z), ComplexMath.Cos(z + 2.0 * Math.PI)));
            }
        }

        [TestMethod]
        public void ComplexNegativeAngles () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 16)) {
                if (Math.Abs(z.Im) > Math.Log(Double.MaxValue / 10.0)) continue; // sin and cos blow up in imaginary part of argument gets too big 
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Sin(-z), -ComplexMath.Sin(z)), String.Format("remainder {0}", ComplexMath.Sin(-a) + ComplexMath.Sin(a)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Cos(-z), ComplexMath.Cos(z)), String.Format("{0} vs. {1}", ComplexMath.Cos(-a), ComplexMath.Cos(a)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Tan(-z), -ComplexMath.Tan(z)));
            }
        }

        [TestMethod]
        public void ComplexAngleAdditionTest () {
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Sin(a + b), ComplexMath.Sin(a) * ComplexMath.Cos(b) + ComplexMath.Cos(a) * ComplexMath.Sin(b)), String.Format("{0} != {1}", ComplexMath.Sin(a + b), ComplexMath.Sin(a) * ComplexMath.Cos(b) + ComplexMath.Cos(a) * ComplexMath.Sin(b)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Cos(a + b), ComplexMath.Cos(a) * ComplexMath.Cos(b) - ComplexMath.Sin(a) * ComplexMath.Sin(b)));
        }

        // Hyperbolic tests

        [TestMethod]
        public void ComplexSinSinh () {
            foreach (Complex x in TestUtilities.GenerateComplexValues(1.0E-2,1.0E2, 10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Sinh(x), -ComplexMath.I * ComplexMath.Sin(ComplexMath.I * x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.I * ComplexMath.Sinh(-ComplexMath.I * x), ComplexMath.Sin(x)));
            }
        }

        [TestMethod]
        public void ComplexCosCosh () {
            foreach (Complex x in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Cosh(x), ComplexMath.Cos(ComplexMath.I * x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Cosh(-ComplexMath.I * x), ComplexMath.Cos(x)));
            }
        }

        [TestMethod]
        public void ComplexTanTanh () {
            foreach (Complex x in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Tanh(x), -ComplexMath.I * ComplexMath.Tan(ComplexMath.I * x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Tanh(-ComplexMath.I * x), -ComplexMath.I * ComplexMath.Tan(x)));
            }
        }


        [TestMethod]
        public void ComplexSinhCoshRelation () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 10)) {
                Complex sinh = ComplexMath.Sinh(z);
                Complex cosh = ComplexMath.Cosh(z);
                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(cosh * cosh, - sinh * sinh, 1));
            }
        }

        // Powers

        [TestMethod()]
        public void ComplexPowMultiplication () {
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Pow(a, 2.0) * ComplexMath.Pow(a, 3.0), ComplexMath.Pow(a, 5.0)));
        }


        [TestMethod]
        public void ComplexPowOneHalf () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 8)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Pow(z, 0.5), ComplexMath.Sqrt(z)));
            }

        }

        // log and exp

        [TestMethod]
        public void ComplexLogSpecialCase () {
            Assert.IsTrue(ComplexMath.Log(Math.E) == 1);
            Assert.IsTrue(ComplexMath.Log(1.0) == 0.0);
            Assert.IsTrue(ComplexMath.Log(-1.0) == ComplexMath.I * Math.PI);
            Assert.IsTrue(ComplexMath.Log(ComplexMath.I) == ComplexMath.I * Math.PI / 2.0);
            Assert.IsTrue(ComplexMath.Log(-ComplexMath.I) == -ComplexMath.I * Math.PI / 2.0);
        }

        [TestMethod]
        public void ComplexLogExp () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Exp(ComplexMath.Log(z)), z));
            }
        }

        [TestMethod]
        public void ComplexLogSum () {
            Complex[] z = TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 10);
            for (int i = 0; i < z.Length; i++) {
                for (int j = 0; j < i; j++) {
                    if (Math.Abs(ComplexMath.Arg(z[i]) + ComplexMath.Arg(z[j])) >= Math.PI) continue;
                    Console.WriteLine("{0} {1}", z[i], z[j]);
                    Console.WriteLine("ln(z1) + ln(z2) = {0}", ComplexMath.Log(z[i]) + ComplexMath.Log(z[j]));
                    Console.WriteLine("ln(z1*z2) = {0}", ComplexMath.Log(z[i] * z[j]));
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(ComplexMath.Log(z[i]), ComplexMath.Log(z[j]), ComplexMath.Log(z[i] * z[j])));
                }
            }
        }

        [TestMethod]
        public void ComplexExpSpecialCase () {
            Assert.IsTrue(ComplexMath.Exp(0.0) == 1.0);
            Assert.IsTrue(ComplexMath.Exp(1.0) == Math.E);
        }

        [TestMethod]
        public void ComplexExpReflection () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Exp(-z), 1.0 / ComplexMath.Exp(z)));
            }
        }

        [TestMethod]
        public void ComplexExpSumTest () {
            foreach (Complex x in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 6)) {
                foreach (Complex y in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 6)) {
                    // don't overflow exp
                    if ((Math.Abs(x.Re) + Math.Abs(y.Re)) > Math.Log(Double.MaxValue / 10.0)) continue;
                    Console.WriteLine("x,y = {0},{1}", x, y);
                    Complex ex = ComplexMath.Exp(x);
                    Complex ey = ComplexMath.Exp(y);
                    Complex xy = x + y;
                    Complex exy = ComplexMath.Exp(xy);

                    // if components of x and y differ by orders of magnitude,
                    // trailing digits will be lost in (x+y), thus changing e^(x+y)
                    // but they will not be lost in e^x and e^y, so agreement will fail due to rounding error
                    // thus we should reduce expected agreement by this ratio
                    double rr = Math.Abs(Math.Log(Math.Abs(x.Re / y.Re)));
                    double ri = Math.Abs(Math.Log(Math.Abs(x.Im / y.Im)));
                    double r = rr;
                    if (ri > r) r = ri;
                    //if (r < 0.0) r = 0.0;
                    double e = TestUtilities.TargetPrecision * Math.Exp(r);
                    Console.WriteLine("e={0}", e);

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(exy, ex*ey, e), String.Format("x={0} y={1} E(x)={2} E(y)={3} E(x)*E(y)={4} E(x+y)={5}", x, y, ex, ey, ex * ey, exy));
                }
            }

        }

        [TestMethod]
        public void ComplexMagnitude () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 10)) {
                double abs = ComplexMath.Abs(z);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(abs * abs, z * z.Conjugate));
            }
        }

        [TestMethod]
        public void ComplexPowExponent () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 6)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Pow(Math.E, z), ComplexMath.Exp(z)));
            }
        }

        [TestMethod]
        public void ComplexPowSpecialCases () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Pow(z, -1), 1.0 / z));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Pow(z, 0), 1.0));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Pow(z, 1), z));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Pow(z, 2), z * z));
            }
        }

        [TestMethod]
        public void ComplexPow () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 8)) {
                Complex p = 1.0;
                for (int k = 1; k < 16; k++) {
                    p *= z;
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Pow(z, k), p));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Pow(z, -k), 1.0 / p));
                }
            }
        }

        [TestMethod]
        public void ComplexRealAgreement () {

            foreach (double x in TestUtilities.GenerateRealValues(0.01, 1000.0, 8)) {

                Assert.IsTrue(ComplexMath.Exp(x) == Math.Exp(x));
                Assert.IsTrue(ComplexMath.Log(x) == Math.Log(x));
                Assert.IsTrue(ComplexMath.Sqrt(x) == Math.Sqrt(x));

                Assert.IsTrue(ComplexMath.Abs(x) == Math.Abs(x));
                Assert.IsTrue(ComplexMath.Arg(x) == 0.0);

            }

        }

        [TestMethod]
        public void ComplexSqrtTest () {
            Console.WriteLine(ComplexSqrt(new Complex(0.0, 0.0)));
            Console.WriteLine(ComplexSqrt(new Complex(2.0, 0.0)));
            Console.WriteLine(ComplexSqrt(new Complex(0.0, 2.0)));
            Console.WriteLine(ComplexSqrt(new Complex(-2.0, 0.0)));
            Console.WriteLine(ComplexSqrt(new Complex(0.0, -2.0)));
            double r = 1.0 / 16.0;
            Console.WriteLine(ComplexSqrt(new Complex(3.0, r)));
            Console.WriteLine(ComplexSqrt(new Complex(3.0, -r)));
            Console.WriteLine(ComplexSqrt(new Complex(r, 3.0)));
            Console.WriteLine(ComplexSqrt(new Complex(-r, 3.0)));
            Console.WriteLine(ComplexSqrt(new Complex(-3.0, r)));
            Console.WriteLine(ComplexSqrt(new Complex(-3.0, -r)));
            Console.WriteLine(ComplexSqrt(new Complex(r, -3.0)));
            Console.WriteLine(ComplexSqrt(new Complex(-r, -3.0)));
        }

        private static Complex ComplexSqrt (Complex z) {

            // We could just call Pow(z, 0.5), but that is slow (involves multiple trig functions)
            // and looses a bit of symmetry (e.g. square roots of pure negative numbers may have tiny real parts because Pi is not exact)

            // Instead we expand (x + i y)^2 = a + i b to obtain a quadratic equations in x^2 and y^2. Solving (and thinking a bit about signs) yields
            //   x = \sqrt{\frac{m + a}{2}}
            //   y = \pm \sqrt{\frac{m - a}{2}}
            // where m = |a + i b| = \sqrt{a^2 + b^2}
            
            // Using this is very straightforward unless |b| << |a|, in which case m ~ |a| and there is a near-cancelation in one of the expressions.
            // In this case we pull out a factor a and use a series expansion of \sqrt{1 + w} - 1 where w = (b/a)^2.

            // In the end the whole algorithm is a little bit more complicated than I would like, but it is much faster than Pow(z, 0.5).

            double x, y;
            if (Math.Abs(z.Im) < Math.Abs(z.Re) / 8.0) {
                double s = Series1(MoreMath.Sqr(z.Im / z.Re));
                if (z.Re > 0.0) {
                    x = Math.Sqrt(z.Re * (s + 2.0) / 2.0);
                    y = Math.Sqrt(z.Re * s / 2.0);
                } else {
                    x = Math.Sqrt(-z.Re * s / 2.0);
                    y = Math.Sqrt(-z.Re * (s + 2.0) / 2.0);
                }
            } else {
                double m = MoreMath.Hypot(z.Re, z.Im);
                x = Math.Sqrt((m + z.Re) / 2.0);
                y = Math.Sqrt((m - z.Re) / 2.0);
            }
            if (z.Im < 0.0) y = -y;
            return (new Complex(x, y));
        }

        // Compute \sqrt{1 + x} - 1 via series

        private static double Series1 (double x) {

            double t = x / 2.0;
            double s = t;

            for (int k = 2; k < 16; k++) {
                double s_old = s;
                t *= - (x / 2.0) * (2 * k - 3) / k;
                s += t;
                if (s == s_old) return(s);
            }

            throw new NonconvergenceException();
        }

    }
}
