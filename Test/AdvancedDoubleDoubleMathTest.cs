using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Extended;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public class AdvancedDoubleDoubleMathTest {

        [TestMethod]
        public void DoubleDoubleErrorFunctionSpecialCases () {
            Assert.IsTrue(AdvancedDoubleDoubleMath.Erf(DoubleDouble.Zero) == DoubleDouble.Zero);
            Assert.IsTrue(AdvancedDoubleDoubleMath.Erfc(DoubleDouble.Zero) == DoubleDouble.One);
        }

        [TestMethod]
        public void DoubleDoubleErrorFunctionReflection () {
            foreach (DoubleDouble x in DoubleDoubleTest.GetRandomDoubleDoubles(1.0E-4, 1.0E4, 8)) {
                DoubleDouble xErf = AdvancedDoubleDoubleMath.Erf(x);
                Assert.IsTrue(-AdvancedDoubleDoubleMath.Erf(-x) == AdvancedDoubleDoubleMath.Erf(x));
            }
        }

        [TestMethod]
        public void DoubleDoubleErrorFunctionComplementiarity () {
            foreach (DoubleDouble x in DoubleDoubleTest.GetRandomDoubleDoubles(1.0E-2, 1.0E2, 8)) {
                DoubleDouble xErf = AdvancedDoubleDoubleMath.Erf(x);
                DoubleDouble xErfc = AdvancedDoubleDoubleMath.Erfc(x);
                Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(xErf + xErfc, DoubleDouble.One));
            }
        }

        [TestMethod]
        public void DoubleDoubleErrorFunctionAgreement () {
            foreach (DoubleDouble x in DoubleDoubleTest.GetRandomDoubleDoubles(1.0E-2, 1.0E2, 8)) {

                DoubleDouble xErf = AdvancedDoubleDoubleMath.Erf(x);
                double xErfAsDouble = (double) xErf;

                double xAsDouble = (double) x;
                double xAsDoubleErf = AdvancedMath.Erf(xAsDouble);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(xErfAsDouble, xAsDoubleErf));

                DoubleDouble xErfc = AdvancedDoubleDoubleMath.Erfc(x);
                Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(xErf + xErfc, DoubleDouble.One));

            }
        }

        [TestMethod]
        public void DoubleDoubleLogGammanSpecialCases () {
            Assert.IsTrue(AdvancedDoubleDoubleMath.LogGamma(0.0) == DoubleDouble.PositiveInfinity);
            Assert.IsTrue(AdvancedDoubleDoubleMath.LogGamma(1.0) == DoubleDouble.Zero);
            Assert.IsTrue(AdvancedDoubleDoubleMath.LogGamma(2.0) == DoubleDouble.Zero);
        }

        [TestMethod]
        public void DoubleDoubleLogGammaRecurrence () {
            foreach (DoubleDouble x in DoubleDoubleTest.GetRandomDoubleDoubles(1.0, 1.0E4, 8)) {
                Assert.IsTrue(DoubleDoubleTest.IsNearlyEqual(
                    DoubleDouble.Log(x) + AdvancedDoubleDoubleMath.LogGamma(x),
                    AdvancedDoubleDoubleMath.LogGamma(x + DoubleDouble.One)
                ));
            }
        }

        [TestMethod]
        public void DoubleDoubleLogGammaAgreement () {
            foreach (DoubleDouble x in DoubleDoubleTest.GetRandomDoubleDoubles(1.0E-2, 1.0E3, 12)) {

                DoubleDouble xLogGamma = AdvancedDoubleDoubleMath.LogGamma(x);
                double xLogGammaAsDouble = (double) xLogGamma;

                double xAsDouble = (double) x;
                double xAsDoubleLogGamma = AdvancedMath.LogGamma(xAsDouble);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(xLogGammaAsDouble, xAsDoubleLogGamma));

            }
        }

    }
}
