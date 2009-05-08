
using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Statistics;
using Meta.Numerics;

namespace Test
{
    
    
    /// <summary>
    ///This is a test class for UncertainValueTest and is intended
    ///to contain all UncertainValueTest Unit Tests
    ///</summary>
    [TestClass()]
    public class UncertainValueTest {


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

        UncertainValue a = new UncertainValue(1.0, 0.5);
        UncertainValue b = new UncertainValue(-2.0, 3.0);

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void UncertainValueConstructorTest () {
            UncertainValue u = new UncertainValue(1.0, -1.0);
        }

        [TestMethod]
        public void UncertainValueArithmeticValuesTest () {
            Assert.IsTrue((a + b).Value == a.Value + b.Value);
            Assert.IsTrue((a - b).Value == a.Value - b.Value);
            Assert.IsTrue((a * b).Value == a.Value * b.Value);
            Assert.IsTrue((a / b).Value == a.Value / b.Value);
            Assert.IsTrue((-a).Value == -a.Value);
        }

        [TestMethod]
        public void UncertainValueAdditionTriangleTest () {
            UncertainValue c = a + b;
            Assert.IsTrue(TestUtilities.IsNearlyEqual(Math.Pow(c.Uncertainty, 2), Math.Pow(a.Uncertainty, 2) + Math.Pow(b.Uncertainty, 2)));
        }

        [TestMethod]
        public void UncertainValueSubtractionTriangleTest () {
            UncertainValue c = a - b;
            Assert.IsTrue(TestUtilities.IsNearlyEqual(Math.Pow(c.Uncertainty, 2), Math.Pow(a.Uncertainty, 2) + Math.Pow(b.Uncertainty, 2)));
        }

        [TestMethod]
        public void UncertainValueMixedArithmeticTest () {
            double b = -2.0;

            UncertainValue c = a * b;
            Assert.IsTrue(c.Value == a.Value * b);
            Assert.IsTrue(c.Uncertainty == a.Uncertainty * Math.Abs(b));

            c = a / b;
            Assert.IsTrue(c.Value == a.Value / b);
            Assert.IsTrue(c.Uncertainty == a.Uncertainty / Math.Abs(b));

            c = a + b;
            Assert.IsTrue(c.Value == a.Value + b);
            Assert.IsTrue(c.Uncertainty == a.Uncertainty);

            c = a - b;
            Assert.IsTrue(c.Value == a.Value - b);
            Assert.IsTrue(c.Uncertainty == a.Uncertainty);

        }

        [TestMethod]
        public void UncertainValueZeroUncertaintiesTest () {

            UncertainValue p = new UncertainValue(3.141592653, 0.0);
            UncertainValue q = new UncertainValue(2.718281828, 0.0);

            UncertainValue sum = p + q;
            Assert.IsTrue(sum.Value == p.Value + q.Value);
            Assert.IsTrue(sum.Uncertainty == 0.0);

            UncertainValue difference = p - q;
            Assert.IsTrue(difference.Value == p.Value - q.Value);
            Assert.IsTrue(difference.Uncertainty == 0.0);

            UncertainValue product = p * q;
            Assert.IsTrue(product.Value == p.Value * q.Value);
            Assert.IsTrue(product.Uncertainty == 0.0);

            UncertainValue quotient = p / q;
            Assert.IsTrue(quotient.Value == p.Value / q.Value);
            Assert.IsTrue(quotient.Uncertainty == 0.0);

            UncertainValue exp = UncertainMath.Exp(p);
            Assert.IsTrue(exp.Value == Math.Exp(p.Value));
            Assert.IsTrue(exp.Uncertainty == 0.0);

            UncertainValue log = UncertainMath.Log(q);
            Assert.IsTrue(log.Value == Math.Log(q.Value));
            Assert.IsTrue(log.Uncertainty == 0.0);

            UncertainValue tan = UncertainMath.Tan(p);
            Assert.IsTrue(tan.Value == Math.Tan(p.Value));
            Assert.IsTrue(tan.Uncertainty == 0.0);

            UncertainValue atan = UncertainMath.Atan(q);
            Assert.IsTrue(atan.Value == Math.Atan(q.Value));
            Assert.IsTrue(atan.Uncertainty == 0.0);


        }

    }
}
