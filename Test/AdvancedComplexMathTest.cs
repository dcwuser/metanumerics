using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Test {


    [TestClass()]
    public class AdvancedComplexMathTest {


        private TestContext testContextInstance;

        public TestContext TestContext {
            get {
                return testContextInstance;
            }
            set {
                testContextInstance = value;
            }
        }

        #region Additional test attributes
        // 
        //You can use the following additional attributes as you write your tests:
        //
        //Use ClassInitialize to run code before running the first test in the class
        //[ClassInitialize()]
        //public static void MyClassInitialize(TestContext testContext)
        //{
        //}
        //
        //Use ClassCleanup to run code after all tests in a class have run
        //[ClassCleanup()]
        //public static void MyClassCleanup()
        //{
        //}
        //
        //Use TestInitialize to run code before running each test
        //[TestInitialize()]
        //public void MyTestInitialize()
        //{
        //}
        //
        //Use TestCleanup to run code after each test has run
        //[TestCleanup()]
        //public void MyTestCleanup()
        //{
        //}
        //
        #endregion

        // Complex error function

        [TestMethod]
        public void ComplexFaddeevaSpecialCaseTest () {
            Assert.IsTrue(AdvancedComplexMath.Faddeeva(0.0) == 1.0);
        }

        [TestMethod]
        public void ComplexFaddeevaConjugationTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(-4,4,50)) {
                if (z.Im * z.Im > Math.Log(Double.MaxValue / 10.0)) continue; // large imaginary values blow up
                Console.WriteLine("z={0}, w(z*)={1}, w*(-z)={2}", z, AdvancedComplexMath.Faddeeva(z.Conjugate), AdvancedComplexMath.Faddeeva(-z).Conjugate);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.Faddeeva(z.Conjugate), AdvancedComplexMath.Faddeeva(-z).Conjugate));
            }
        }

        [TestMethod]
        public void ComplexFaddeevaDawsonTest () {
            foreach (double x in TestUtilities.GenerateRealValues(-4,4,10)) {
                Complex w = AdvancedComplexMath.Faddeeva(x);
                double F = AdvancedMath.Dawson(x);
                Console.WriteLine("x={0} F={1} w={2}", x, F, w);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(w.Re, Math.Exp(-x * x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(w.Im, (2.0/Math.Sqrt(Math.PI)) * AdvancedMath.Dawson(x)));
            }
        }

        // Complex Gamma Function

        [TestMethod]
        public void ComplexGammaConjugationTest () {
            // limited to 10^-2 to 10^2 to avoid overflow
            foreach (Complex z in TestUtilities.GenerateComplexValues(-2, 2, 10)) {
                Console.WriteLine("z={0} G(z*)={1} G*(z)={2}", z, AdvancedComplexMath.Gamma(z.Conjugate), AdvancedComplexMath.Gamma(z).Conjugate);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.Gamma(z.Conjugate), AdvancedComplexMath.Gamma(z).Conjugate));
            }
        }

        [TestMethod]
        public void ComplexGammaRecurranceTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(-2, 2, 20)) {
                Complex G = AdvancedComplexMath.Gamma(z);
                Complex Gz = G * z;
                Complex GP = AdvancedComplexMath.Gamma(z + 1.0);
                Console.WriteLine("z={0:g16} G(z)={1:g16} G(z)*z={2:g16} G(z+1)={3:g16}", z, G, Gz, GP);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(Gz,GP));
            }
        }

        // we need to implement complex powers before we can do this test

        /*
        [TestMethod]
        public void ComplexGammaDuplicationTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(-2, 1.8, 10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.Gamma(2.0 * z), ComplexMath.Pow(2.0, 2.0 * z - 0.5) * AdvancedComplexMath.Gamma(z) * AdvancedComplexMath.Gamma(z + 0.5) / Math.Sqrt(2.0 * Math.PI)), String.Format("x={0}, Gamma(x)={1}, Gamma(2x) = {2}", z, AdvancedMath.Gamma(z), AdvancedMath.Gamma(2.0 * z)));
            }
        }
        */

        [TestMethod]
        public void ComplexGammaInequalityTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(-2, 2, 10)) {
                Console.WriteLine(z);
                Assert.IsTrue(ComplexMath.Abs(AdvancedComplexMath.Gamma(z)) <= Math.Abs(AdvancedMath.Gamma(z.Re)));
            }
        }

        [TestMethod]
        public void ComplexGammaKnownLinesTest () {
            // to avoid problems with trig functions of large arguments, don't let x get too big
            foreach (double x in TestUtilities.GenerateRealValues(-4, 2, 15)) {
                Console.WriteLine(x);
                double px = Math.PI * x;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.Gamma(new Complex(0.0, x)) * AdvancedComplexMath.Gamma(new Complex(0.0, -x)), Math.PI / x / Math.Sinh(px)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.Gamma(new Complex(0.5, x)) * AdvancedComplexMath.Gamma(new Complex(0.5, -x)), Math.PI / Math.Cosh(px)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.Gamma(new Complex(1.0, x)) * AdvancedComplexMath.Gamma(new Complex(1.0, -x)), px / Math.Sinh(px)));
            }
        }

        // Log Gamma

        [TestMethod]
        public void ComplexLogGammaConjugationTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(-4, 4, 20)) {
                if (z.Re <= 0.0) continue;
                Console.WriteLine("z={0} lnG(z*)={1} lnG*(z)={2}", z, AdvancedComplexMath.LogGamma(z.Conjugate), AdvancedComplexMath.LogGamma(z).Conjugate);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.LogGamma(z.Conjugate), AdvancedComplexMath.LogGamma(z).Conjugate));
            }
        }

        [TestMethod]
        public void ComplexLogGammaAgreementTest () {
            foreach (double x in TestUtilities.GenerateRealValues(-4, 4, 20)) {
                //Console.WriteLine("{0} {1}", AdvancedComplexMath.LogGamma(x), AdvancedMath.LogGamma(x));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.LogGamma(x), AdvancedMath.LogGamma(x)));
            }
        }

        // periodicity in imaginary part of ln z means that LogGamma recurrence and LogGamma duplication tests fail

    }
}
