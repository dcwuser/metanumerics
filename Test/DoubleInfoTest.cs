using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using FluentAssertions;

using Meta.Numerics.Extended;

namespace Test {

    [TestClass]
    public class DoubleInfoTests {

        [TestMethod]
        public void DoubleInfoZeroProperties () {

            DoubleInfo zero = new DoubleInfo(0.0);
            zero.IsFinite.Should().BeTrue();
            zero.IsInfinite.Should().BeFalse();
            zero.IsNaN.Should().BeFalse();
            zero.IsNegative.Should().BeFalse();
            zero.IsSubnormal.Should().BeFalse(); // This is debatable.
            zero.IsZero.Should().BeTrue();

            zero.Value.Should().Be(0.0);
            zero.Next.Value.Should().Be(Double.Epsilon);

            zero.Exponent.Should().Be(0);
            zero.Mantissa.Should().Be(0L);

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

        [TestMethod]
        public void DoubleInfoPowersOfTwo () {

            foreach (int n in TestUtilities.GenerateIntegerValues(1, 1000, 4)) {

                DoubleInfo dp = DoubleInfo.For(Math.Pow(2, n));
                dp.Mantissa.Should().Be(1);
                dp.Exponent.Should().Be(n);

                DoubleInfo dm = DoubleInfo.For(Math.Pow(2, -n));
                dm.Mantissa.Should().Be(1);
                dm.Exponent.Should().Be(-n);

            }

        }

    }

}