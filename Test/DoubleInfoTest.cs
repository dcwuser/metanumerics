using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Extended;

namespace Test {

    [TestClass]
    public class DoubleInfoTests {

        [TestMethod]
        public void DoubleInfoZeroProperties () {

            DoubleInfo zero = new DoubleInfo(0.0);
            Assert.IsTrue(zero.IsFinite);
            Assert.IsFalse(zero.IsInfinite);
            Assert.IsFalse(zero.IsNaN);
            Assert.IsFalse(zero.IsNegative);
            Assert.IsFalse(zero.IsSubnormal); // This is debatable.
            Assert.IsTrue(zero.IsZero);

            Assert.IsTrue(zero.Value == 0.0);
            Assert.IsTrue(zero.Next.Value == Double.Epsilon);

            Assert.IsTrue(zero.Exponent == 0); // This isn't mathematically necessary, but IEEE754 requires it.
            Assert.IsTrue(zero.Mantissa == 0L);
        }

        [TestMethod]
        public void DoubleInfoOneProperties () {

            DoubleInfo one = new DoubleInfo(1.0);
            Assert.IsTrue(one.IsFinite);
            Assert.IsFalse(one.IsInfinite);
            Assert.IsFalse(one.IsNaN);
            Assert.IsFalse(one.IsNegative);
            Assert.IsFalse(one.IsSubnormal);
            Assert.IsFalse(one.IsZero);

            Assert.IsTrue(one.Value == 1.0);
            Assert.IsTrue(one.Next.Value != 1.0);
            Assert.IsTrue(one.Previous.Value != 1.0);

            Assert.IsTrue(one.Exponent == 0);
            Assert.IsTrue(one.Mantissa == 1L);

        }

        [TestMethod]
        public void DoubleInfoNegativeInfinityProperties () {

            DoubleInfo negInf = new DoubleInfo(Double.NegativeInfinity);
            Assert.IsFalse(negInf.IsFinite);
            Assert.IsTrue(negInf.IsInfinite);
            Assert.IsFalse(negInf.IsNaN);
            Assert.IsTrue(negInf.IsNegative);
            Assert.IsFalse(negInf.IsSubnormal);
            Assert.IsFalse(negInf.IsZero);

            Assert.IsTrue(negInf.Value == Double.NegativeInfinity);
            Assert.IsTrue(negInf.Next.Value == Double.MinValue);

        }

        [TestMethod]
        public void DoubleInfoNaNProperties () {

            DoubleInfo nan = new DoubleInfo(Double.NaN);
            Assert.IsFalse(nan.IsFinite);
            Assert.IsFalse(nan.IsInfinite);
            Assert.IsTrue(nan.IsNaN);
            Assert.IsFalse(nan.IsSubnormal);
            Assert.IsFalse(nan.IsZero);

            Assert.IsTrue(Double.IsNaN(nan.Value));

        }

        [TestMethod]
        public void DoubleInfoEpsilonProperties () {

            DoubleInfo epsilon = new DoubleInfo(Double.Epsilon);
            Assert.IsTrue(epsilon.IsFinite);
            Assert.IsFalse(epsilon.IsInfinite);
            Assert.IsFalse(epsilon.IsNaN);
            Assert.IsTrue(epsilon.IsSubnormal);
            Assert.IsFalse(epsilon.IsZero);

            Assert.IsTrue(epsilon.Value == Double.Epsilon);

        }

    }

}