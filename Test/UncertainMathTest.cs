using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;


namespace Test {
    

    [TestClass]
    public class UncertainMathTest {

        UncertainValue a = new UncertainValue(2.0, 0.5);

        [TestMethod]
        public void UncertainMathValueAgreement () {
            // We need to move to IsNearlyEqual because exact equality fails in release mode, even though all-digit printed values are the same
            // I suspect this has to do with the optomized comparison including bits beyond those that strictly belong to the value and are random
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Cos(a).Value, Math.Cos(a.Value)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Sin(a).Value, Math.Sin(a.Value)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Tan(a).Value, Math.Tan(a.Value)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Cosh(a).Value, Math.Cosh(a.Value)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Sinh(a).Value, Math.Sinh(a.Value)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Tanh(a).Value, Math.Tanh(a.Value)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Log(a).Value, Math.Log(a.Value)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Exp(a).Value, Math.Exp(a.Value)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Sqrt(a).Value, Math.Sqrt(a.Value)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Atan(a).Value, Math.Atan(a.Value)));
        }

        [TestMethod]
        public void UncertainMathSqrtSquare () {
            UncertainValue sa = UncertainMath.Sqrt(a);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(sa.RelativeUncertainty, a.RelativeUncertainty / 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Pow(sa, 2.0).Value, a.Value));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(UncertainMath.Pow(sa, 2.0).Uncertainty, a.Uncertainty));
        }

        [TestMethod]
        public void UncertainMathLogExp () {
            UncertainValue la = UncertainMath.Log(a);
            Assert.IsTrue(la.Uncertainty == a.RelativeUncertainty);
            Assert.IsTrue(UncertainMath.Exp(la) == a);
        }

        [TestMethod]
        public void UncertainMathCosACos () {
            UncertainValue y = new UncertainValue(0.5, 0.1);
            UncertainValue x = UncertainMath.Acos(y);
            UncertainValue cx = UncertainMath.Cos(x);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(cx.Value, y.Value));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(cx.Uncertainty, y.Uncertainty));
        }

    }
}
