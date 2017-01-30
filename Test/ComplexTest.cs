using System;
using Meta.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Test
{
    
    [TestClass]
    public class ComplexTest {

        private Complex a = new Complex(1.5, -2.2);
        private Complex b = new Complex(-7.4, 3.8);
        private Complex c = new Complex(5.9, 6.0);

        [TestMethod]
        public void ComplexComponents () {
            Complex z = new Complex(-1.0, 2.0);
            Assert.IsTrue(z.Re == -1.0);
            Assert.IsTrue(z.Im == 2.0);
        }

        [TestMethod]
        public void ComplexEquality () {

            Complex ac = a;

            // Equality operator
            Assert.IsTrue(a == ac);
            Assert.IsTrue(ac == a);
            Assert.IsFalse(a == b);
            Assert.IsFalse(b == a);

            // Inequality operator
            Assert.IsFalse(a != ac);
            Assert.IsFalse(ac != a);
            Assert.IsTrue(a != b);
            Assert.IsTrue(b != a);

            // Equals method
            Assert.IsTrue(a.Equals(ac));
            Assert.IsFalse(a.Equals(b));
            Assert.IsFalse(a.Equals(new object()));
            Assert.IsFalse(a.Equals(null));

            // Hash
            Assert.IsTrue(a.GetHashCode() == ac.GetHashCode());

        }

        [TestMethod]
        public void ComplexConjugation () {
            Assert.IsTrue(a.Conjugate.Re == a.Re);
            Assert.IsTrue(a.Conjugate.Im == -a.Im);
            Assert.IsTrue((a * a.Conjugate).Im == 0.0);
        }

        [TestMethod]
        public void ComplexNegation () {
            Assert.IsTrue((a + (-a)) == 0);
        }

        [TestMethod]
        public void ComplexSubtraction () {
            Assert.IsTrue((a - a) == Complex.Zero);
            Assert.IsTrue((0 - a) == (-a));
            Assert.IsTrue((a - b) == -(b - a));
        }


        [TestMethod]
        public void ComplexMultiplication () {
            // Commutative
            Assert.IsTrue(a * b == b * a);
            // Associative
            Assert.IsTrue(TestUtilities.IsNearlyEqual(a * (b * c), (a * b) * c));
            // Distributive
            Assert.IsTrue(a * (b + c) == (a * b + a * c));
            // Multiplicative identity
            Assert.IsTrue(a * Complex.One == a);
            // Multiplication by zero
            Assert.IsTrue(a * Complex.Zero == Complex.Zero);
        }

        [TestMethod]
        public void ComplexDivision () {
            // Multiplicative inverse
            Assert.IsTrue(a / a == Complex.One);
            // Division identity
            Assert.IsTrue(a / Complex.One == a);
            // Associative
            Assert.IsTrue(TestUtilities.IsNearlyEqual(a / b / c, a / c / b));
            // Distributive
            Assert.IsTrue(TestUtilities.IsNearlyEqual((b + c) / a, b / a + c / a));
        }

        [TestMethod]
        public void ComplexExplicitCast () {
            Complex a = new Complex(1.0, 0.0);
            double b = (double) a;
            Assert.IsTrue(b == a.Re);
        }

        [TestMethod]
        [ExpectedException(typeof(InvalidCastException))]
        public void ComplexExplicitCastFailure () {
            Complex a = new Complex(1.0,1.0);
            double b = (double) a;
        }

        [TestMethod]
        public void ComplexAddition () {
            // Commutative
            Assert.IsTrue(a + b == b + a);
            // Associative
            Assert.IsTrue(TestUtilities.IsNearlyEqual((a + b) + c, a + (b + c)));
            // Relation of subtration to negation and additon
            Assert.IsTrue(a + (-b) == a - b);
        }


        [TestMethod]
        public void ComplexEquals () {
            Assert.IsTrue(Complex.Equals(a, a));
            Assert.IsFalse(Complex.Equals(a, b));
        }

        [TestMethod]
        public void ComplexDivisionDifficult () {

            // Difficult complex divisions from Baudin & Smith (https://arxiv.org/pdf/1210.4539v2.pdf)
            // The naive denominator would overflow for many of these.
 
            Complex r1 = (new Complex(1.0, 1.0)) / (new Complex(1.0, Double.MaxValue));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(r1.Re, 1.0 / Double.MaxValue));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(r1.Im, -1.0 / Double.MaxValue));

            Complex r2 = (new Complex(1.0, 1.0)) / (new Complex(Math.Pow(2.0, -1023), Math.Pow(2.0, -1023)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(r2.Re, Math.Pow(2.0, 1023)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(r2.Im, 0.0));

            Complex r3 = new Complex(Math.Pow(2.0, 1023), Math.Pow(2.0, -1023)) / (new Complex(Math.Pow(2.0, 677), Math.Pow(2.0, -677)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(r3.Re, Math.Pow(2.0, 346)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(r3.Im, -Math.Pow(2.0, -1008)));

            Complex r4  = (new Complex(Math.Pow(2.0, 1023), Math.Pow(2.0, 1023))) / (new Complex(1.0, 1.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(r4.Re, Math.Pow(2.0, 1023)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(r4.Im, 0.0));

            Complex r5 = (new Complex(Math.Pow(2.0, 1020), Math.Pow(2.0, -844))) / (new Complex(Math.Pow(2.0, 656), Math.Pow(2.0, -780)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(r5.Re, Math.Pow(2.0, 364)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(r5.Im, -Math.Pow(2.0, -1072)));

        }

        [TestMethod]
        public void ComplexPowSpecialValues () {

            Assert.IsTrue(ComplexMath.Pow(Complex.Zero, Complex.Zero) == Complex.One);
            Assert.IsTrue(ComplexMath.Pow(Complex.One, Complex.Zero) == Complex.One);
            Assert.IsTrue(ComplexMath.Pow(-Complex.One, Complex.Zero) == Complex.One);
            Assert.IsTrue(ComplexMath.Pow(Complex.I, Complex.Zero) == Complex.One);
            Assert.IsTrue(ComplexMath.Pow(a, Complex.Zero) == Complex.One);

            Assert.IsTrue(ComplexMath.Pow(Complex.Zero, Complex.One) == Complex.Zero);
            Assert.IsTrue(ComplexMath.Pow(Complex.One, Complex.One) == Complex.One);
            Assert.IsTrue(ComplexMath.Pow(-Complex.One, Complex.One) == -Complex.One);
            Assert.IsTrue(ComplexMath.Pow(Complex.I, Complex.One) == Complex.I);
            Assert.IsTrue(ComplexMath.Pow(a, Complex.One) == a);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Pow(Complex.I, Complex.I), Math.Exp(-Math.PI / 2.0)));

        }

    }
}
