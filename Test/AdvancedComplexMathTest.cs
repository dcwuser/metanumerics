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
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4,1.0E4,50)) {
                if (z.Im * z.Im > Math.Log(Double.MaxValue / 10.0)) continue; // large imaginary values blow up
                Console.WriteLine("z={0}, w(z*)={1}, w*(-z)={2}", z, AdvancedComplexMath.Faddeeva(z.Conjugate), AdvancedComplexMath.Faddeeva(-z).Conjugate);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.Faddeeva(z.Conjugate), AdvancedComplexMath.Faddeeva(-z).Conjugate));
            }
        }

        [TestMethod]
        public void ComplexFaddeevaDawsonTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4,1.0E4,10)) {
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
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 10)) {
                Console.WriteLine("z={0} G(z*)={1} G*(z)={2}", z, AdvancedComplexMath.Gamma(z.Conjugate), AdvancedComplexMath.Gamma(z).Conjugate);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.Gamma(z.Conjugate), AdvancedComplexMath.Gamma(z).Conjugate));
            }
        }

        [TestMethod]
        public void ComplexGammaRecurranceTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 20)) {
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
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 10)) {
                Console.WriteLine(z);
                Assert.IsTrue(ComplexMath.Abs(AdvancedComplexMath.Gamma(z)) <= Math.Abs(AdvancedMath.Gamma(z.Re)));
            }
        }

        [TestMethod]
        public void ComplexGammaKnownLinesTest () {
            // to avoid problems with trig functions of large arguments, don't let x get too big
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 100.0, 15)) {
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
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 20)) {
                if (z.Re <= 0.0) continue;
                Console.WriteLine("z={0} lnG(z*)={1} lnG*(z)={2}", z, AdvancedComplexMath.LogGamma(z.Conjugate), AdvancedComplexMath.LogGamma(z).Conjugate);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.LogGamma(z.Conjugate), AdvancedComplexMath.LogGamma(z).Conjugate));
            }
        }

        [TestMethod]
        public void ComplexLogGammaAgreementTest () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 20)) {
                //Console.WriteLine("{0} {1}", AdvancedComplexMath.LogGamma(x), AdvancedMath.LogGamma(x));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.LogGamma(x), AdvancedMath.LogGamma(x)));
            }
        }

        // periodicity in imaginary part of ln z means that LogGamma recurrence and LogGamma duplication tests fail

        [TestMethod]
        public void ComplexPsiRecurrenceTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 12)) {
                //Complex z = new Complex(0.398612, 888.865);
                //double t = Math.Exp(900.0);
                //Console.WriteLine("{0} {1}", t, 1.0 / t);
                Console.WriteLine("z={0}", z);
                //Console.WriteLine("Psi(1-z)={0}", AdvancedComplexMath.Psi(1.0 - z));
                //Console.WriteLine("Tan(pi z)={0}", ComplexMath.Tan(Math.PI * z));
                //Console.WriteLine("{0} {1}", z, AdvancedComplexMath.Psi(z));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedComplexMath.Psi(z + 1.0),
                    AdvancedComplexMath.Psi(z) + 1.0 / z
                ));
            }
        }

        [TestMethod]
        public void ComplexPsiDuplicationTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E3, 12)) {
                Console.WriteLine(z);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedComplexMath.Psi(2.0 * z),
                    AdvancedComplexMath.Psi(z) / 2.0 + AdvancedComplexMath.Psi(z + 0.5) / 2.0 + Math.Log(2.0)
                ));
            }
        }

        [TestMethod]
        public void ComplexPsiConjugationTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E3, 12)) {
                Console.WriteLine(z);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedComplexMath.Psi(z.Conjugate),
                    AdvancedComplexMath.Psi(z).Conjugate
                ));
            }
        }

        [TestMethod]
        public void ComplexPsiImaginaryPartsTest () {
            // don't get these values get too big, because we loose trig function accuracy
            foreach (double y in TestUtilities.GenerateRealValues(0.01, 100.0, 10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedComplexMath.Psi(new Complex(0.0, y)).Im,
                    (Math.PI / Math.Tanh(Math.PI * y) + 1.0 / y) / 2.0
                ));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedComplexMath.Psi(new Complex(0.5, y)).Im,
                    Math.PI * Math.Tanh(Math.PI * y) / 2.0
                ));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedComplexMath.Psi(new Complex(1.0, y)).Im,
                    (Math.PI / Math.Tanh(Math.PI * y) - 1.0 / y) / 2.0
                ));
            }
        }

        [TestMethod]
        public void ComplexDiLogSymmetryTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 12)) {

                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                    AdvancedComplexMath.DiLog(z), AdvancedComplexMath.DiLog(-z),
                    AdvancedComplexMath.DiLog(z * z) / 2.0
                ));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    - AdvancedComplexMath.DiLog(1.0 / z),
                    AdvancedComplexMath.DiLog(z) + ComplexMath.Log(-z) * ComplexMath.Log(-z) / 2.0 + Math.PI * Math.PI / 6.0
                ));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    - AdvancedComplexMath.DiLog(1.0 - z),
                    AdvancedComplexMath.DiLog(z) + ComplexMath.Log(z) * ComplexMath.Log(1.0 - z) - Math.PI * Math.PI / 6.0
                ));


            }
        }

        [TestMethod]
        public void ComplexDiLogUnitCircleTest () {

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedComplexMath.DiLog(1.0),
                Math.PI * Math.PI / 6.0
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedComplexMath.DiLog(ComplexMath.I),
                new Complex(-Math.PI * Math.PI / 48.0, AdvancedMath.Catalan)
             ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedComplexMath.DiLog(-1.0),
                -Math.PI * Math.PI / 12.0
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                AdvancedComplexMath.DiLog(-ComplexMath.I),
                new Complex(-Math.PI * Math.PI / 48.0, -AdvancedMath.Catalan)
             ));

        }

        [TestMethod]
        public void ComplexDiLogAgreementTest () {
            // use negative arguments b/c positive real DiLog function complex for x > 1
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.DiLog(-x), AdvancedComplexMath.DiLog(-x)
                ));
            }
        }

        [TestMethod]
        public void ComplexDiLogConjugationTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-1, 1.0E1, 10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedComplexMath.DiLog(z.Conjugate), AdvancedComplexMath.DiLog(z).Conjugate
                ));
            }
        }

        [TestMethod]
        public void ComplexDiLogBranchCutTest () {
            foreach (double z in TestUtilities.GenerateRealValues(1.0, 8.0, 10)) {

                // Imaginary part should be positive just above the cut, negative just below the cut

                Complex Dp = AdvancedComplexMath.DiLog(z + ComplexMath.I * TestUtilities.TargetPrecision / 8.0);
                Complex Dm = AdvancedComplexMath.DiLog(z - ComplexMath.I * TestUtilities.TargetPrecision / 8.0);

                Assert.IsTrue(Dp.Im > 0.0);
                Assert.IsTrue(Dm.Im < 0.0);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(Dp.Re, Dm.Re));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(Dp.Im, -Dm.Im));

            }
        }

        [TestMethod]
        public void ComplexRiemannTest () {

            Complex z = new Complex(10.0, 10.0);
            Complex r = AdvancedComplexMath.Riemann_Euler(z);
            Console.WriteLine(r);

            double r1 = AdvancedMath.RiemannZeta(z.Re);
            Console.WriteLine(r1);

        }

    }
}
