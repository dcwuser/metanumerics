using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public class AdvancedComplexMathTest {

        // Complex error function

        [TestMethod]
        public void ComplexErfSpecialCase () {
            Assert.IsTrue(AdvancedComplexMath.Erf(0.0) == 0.0);
        }

        [TestMethod]
        public void ComplexErfSymmetries () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-1, 1.0E3, 16)) {
                // for large imaginary parts, erf overflows, so don't try them
                if (z.Im * z.Im > 500.0) continue;
                //Console.WriteLine("{0} {1} {2} {3}", z, AdvancedComplexMath.Erf(z), AdvancedComplexMath.Erf(-z), AdvancedComplexMath.Erf(z.Conjugate));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedComplexMath.Erf(-z), -AdvancedComplexMath.Erf(z)
                ));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedComplexMath.Erf(z.Conjugate), AdvancedComplexMath.Erf(z).Conjugate
                ));
            }
        }

        [TestMethod]
        public void ComplexErfFresnel () {
            // don't let x get too big or we run into problem of accruacy for arguments of large trig-like functions
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-1, 1.0E2, 16)) {
                Complex z = Math.Sqrt(Math.PI) * (1.0 - ComplexMath.I) * x / 2.0;
                Complex w = (1.0 + ComplexMath.I) * AdvancedComplexMath.Erf(z) / 2.0;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(w.Re, AdvancedMath.FresnelC(x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(w.Im, AdvancedMath.FresnelS(x)));
            }
        }

        [TestMethod]
        public void ComplexErfFaddevaAgreement () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-1, 1.0E2, 16)) {
                Complex w = AdvancedComplexMath.Faddeeva(z);
                Complex erf = AdvancedComplexMath.Erf(-ComplexMath.I * z);
                Complex erfc = 1.0 - erf;
                Complex f = ComplexMath.Exp(z * z);
                Console.WriteLine("z={0} w={1} erf={2} erfc={3} f={4} w'={5} erf'={6}", z, w, erf, erfc, f, erfc / f, 1.0 - f * w);
                if (Double.IsInfinity(f.Re) || Double.IsInfinity(f.Im) || f == 0.0) continue;
                //Console.WriteLine("{0} {1}", TestUtilities.IsNearlyEqual(w, erfc / f), TestUtilities.IsNearlyEqual(erf, 1.0 - f * w));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    erf, 1.0 - f * w
                ));
            }
        }

        [TestMethod]
        public void ComplexFaddeevaSpecialCase () {
            Assert.IsTrue(AdvancedComplexMath.Faddeeva(0.0) == 1.0);
        }

        [TestMethod]
        public void ComplexFaddeevaConjugation () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4,1.0E4,50)) {
                if (z.Im * z.Im > Math.Log(Double.MaxValue / 10.0)) continue; // large imaginary values blow up
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.Faddeeva(z.Conjugate), AdvancedComplexMath.Faddeeva(-z).Conjugate));
            }
        }

        [TestMethod]
        public void ComplexFaddeevaDawson () {
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
        public void ComplexGammaConjugation () {
            // limited to 10^-2 to 10^2 to avoid overflow
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 10)) {
                Console.WriteLine("z={0} G(z*)={1} G*(z)={2}", z, AdvancedComplexMath.Gamma(z.Conjugate), AdvancedComplexMath.Gamma(z).Conjugate);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.Gamma(z.Conjugate), AdvancedComplexMath.Gamma(z).Conjugate));
            }
        }

        [TestMethod]
        public void ComplexGammaRecurrance () {
            // fails when extended to more numbers due to a loss of last 4 digits of accuracy in an extreme case; look into it
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 10)) {
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
        public void ComplexGammaInequality () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 10)) {
                Console.WriteLine(z);
                Assert.IsTrue(ComplexMath.Abs(AdvancedComplexMath.Gamma(z)) <= Math.Abs(AdvancedMath.Gamma(z.Re)));
            }
        }

        [TestMethod]
        public void ComplexGammaKnownLines () {
            // to avoid problems with trig functions of large arguments, don't let x get too big
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E2, 16)) {
                Console.WriteLine(x);
                double px = Math.PI * x;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.Gamma(new Complex(0.0, x)) * AdvancedComplexMath.Gamma(new Complex(0.0, -x)), Math.PI / x / Math.Sinh(px)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.Gamma(new Complex(0.5, x)) * AdvancedComplexMath.Gamma(new Complex(0.5, -x)), Math.PI / Math.Cosh(px)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.Gamma(new Complex(1.0, x)) * AdvancedComplexMath.Gamma(new Complex(1.0, -x)), px / Math.Sinh(px)));
            }
        }

        // Log Gamma

        [TestMethod]
        public void ComplexLogGammaConjugation () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 24)) {
                if (z.Re <= 0.0) continue;
                Console.WriteLine("z={0} lnG(z*)={1} lnG*(z)={2}", z, AdvancedComplexMath.LogGamma(z.Conjugate), AdvancedComplexMath.LogGamma(z).Conjugate);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedComplexMath.LogGamma(z.Conjugate),
                    AdvancedComplexMath.LogGamma(z).Conjugate
                ));
            }
        }

        [TestMethod]
        public void ComplexLogGammaAgreement () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-4, 1.0E4, 24)) {
                //Console.WriteLine("{0} {1}", AdvancedComplexMath.LogGamma(x), AdvancedMath.LogGamma(x));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(AdvancedComplexMath.LogGamma(x), AdvancedMath.LogGamma(x)));
            }
        }

        // periodicity in imaginary part of ln z means that LogGamma recurrence and LogGamma duplication tests fail

        [TestMethod]
        public void ComplexPsiRecurrence () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 16)) {
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
        public void ComplexPsiDuplication () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E3, 12)) {
                Console.WriteLine(z);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedComplexMath.Psi(2.0 * z),
                    AdvancedComplexMath.Psi(z) / 2.0 + AdvancedComplexMath.Psi(z + 0.5) / 2.0 + Math.Log(2.0)
                ));
            }
        }

        [TestMethod]
        public void ComplexPsiConjugation () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E3, 12)) {
                Console.WriteLine(z);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedComplexMath.Psi(z.Conjugate),
                    AdvancedComplexMath.Psi(z).Conjugate
                ));
            }
        }

        [TestMethod]
        public void ComplexPsiImaginaryParts () {
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
        public void ComplexDiLogSymmetry () {
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
        public void ComplexDiLogUnitCircle () {

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
        public void ComplexDiLogAgreement () {
            // use negative arguments b/c positive real DiLog function complex for x > 1
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.DiLog(-x), AdvancedComplexMath.DiLog(-x)
                ));
            }
        }

        [TestMethod]
        public void ComplexDiLogConjugation () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-1, 1.0E1, 8)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedComplexMath.DiLog(z.Conjugate), AdvancedComplexMath.DiLog(z).Conjugate
                ));
            }
        }

        [TestMethod]
        public void ComplexDiLogBranchCut () {
            foreach (double z in TestUtilities.GenerateRealValues(1.0, 1.0E1, 8)) {

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
        public void ComplexDiLogClausen () {

            foreach (double t in TestUtilities.GenerateUniformRealValues(0.0, 2.0 * Math.PI, 4)) {

                Complex z = AdvancedComplexMath.DiLog(ComplexMath.Exp(ComplexMath.I * t));

                Assert.IsTrue(TestUtilities.IsNearlyEqual(z.Re, Math.PI * Math.PI / 6.0 - t * (2.0 * Math.PI - t) / 4.0));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(z.Im, AdvancedMath.Clausen(t)));

            }

        }

        [TestMethod]
        public void ComplexEinAgreement () {
            foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E4, 16)) {
                Console.WriteLine("Ei {0} {1} {2}", x, AdvancedMath.IntegralEi(x), AdvancedMath.EulerGamma + Math.Log(x) - AdvancedComplexMath.Ein(-x));
                if (x < Math.Log(Double.MaxValue / 10.0)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        AdvancedMath.IntegralEi(x),
                        AdvancedMath.EulerGamma + Math.Log(x) - AdvancedComplexMath.Ein(-x)
                    ));
                }
                Console.WriteLine("E1 {0} {1} {2}", x, AdvancedMath.IntegralE(1, x) + AdvancedMath.EulerGamma + Math.Log(x), AdvancedComplexMath.Ein(x));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    AdvancedMath.IntegralE(1, x) + AdvancedMath.EulerGamma + Math.Log(x),
                    AdvancedComplexMath.Ein(x)
                ));
            }

        }

    }
}
