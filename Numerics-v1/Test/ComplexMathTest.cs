using System;
using Meta.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
namespace Test
{
    
    
    /// <summary>
    ///This is a test class for ComplexMathTest and is intended
    ///to contain all ComplexMathTest Unit Tests
    ///</summary>
    [TestClass()]
    public class ComplexMathTest {


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


        private Complex a = new Complex(3.0, 4.0);
        private Complex b = new Complex(1.5, -2.3);

        [TestMethod()]
        public void ComplexITest () {
            Assert.AreEqual<double>(ComplexMath.I.Re, 0.0);
            Assert.AreEqual<double>(ComplexMath.I.Im, 1.0);
        }

        // Square root

        [TestMethod]
        public void ComplexSqrtTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 10)) {
                Complex sz = ComplexMath.Sqrt(z);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(sz * sz,z));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Arg(z) / 2.0, ComplexMath.Arg(sz)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(Math.Sqrt(ComplexMath.Abs(z)), ComplexMath.Abs(sz)));
            }
        }

        // Trig functions

        [TestMethod()]
        public void ComplexSinCosTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 10)) {
                if (Math.Abs(z.Im) > Math.Log(Double.MaxValue / 10.0)) continue; // sin and cos blow up in imaginary part of argument gets too big 
                Complex sin = ComplexMath.Sin(z);
                Complex cos = ComplexMath.Cos(z);
                Console.WriteLine("z={0} s={1} c={2} s^2+c^2={3}", z, sin, cos, sin*sin + cos*cos);
                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(sin * sin, cos * cos, 1));
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
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 6)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Sin(z), ComplexMath.Sin(z + 2.0 * Math.PI)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Cos(z), ComplexMath.Cos(z + 2.0 * Math.PI)));
            }
        }

        [TestMethod]
        public void ComplexNegativeAnglesTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 10)) {
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
        public void ComplexSinSinhTest () {
            foreach (Complex x in TestUtilities.GenerateComplexValues(1.0E-2,1.0E2,10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Sinh(x), -ComplexMath.I * ComplexMath.Sin(ComplexMath.I * x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.I * ComplexMath.Sinh(-ComplexMath.I * x), ComplexMath.Sin(x)));
            }
        }

        [TestMethod]
        public void ComplexCosCoshTest () {
            foreach (Complex x in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Cosh(x), ComplexMath.Cos(ComplexMath.I * x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Cosh(-ComplexMath.I * x), ComplexMath.Cos(x)));
            }
        }

        [TestMethod]
        public void ComplexTanTanhTest () {
            foreach (Complex x in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Tanh(x), -ComplexMath.I * ComplexMath.Tan(ComplexMath.I * x)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Tanh(-ComplexMath.I * x), -ComplexMath.I * ComplexMath.Tan(x)));
            }
        }


        [TestMethod]
        public void ComplexSinhCoshTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 10)) {
                Complex sinh = ComplexMath.Sinh(z);
                Complex cosh = ComplexMath.Cosh(z);
                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(cosh * cosh, - sinh * sinh, 1));
            }
        }

        // Powers

        [TestMethod()]
        public void PowTest () {
            Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Pow(a, 2.0) * ComplexMath.Pow(a, 3.0), ComplexMath.Pow(a, 5.0)));
        }

        [TestMethod]
        public void ComplexPowSpecialCaseTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 10)) {
                Assert.IsTrue(ComplexMath.Pow(z, 0) == 1, "0");
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Pow(z, 1), z), String.Format("z={0} z^1={1}", z, ComplexMath.Pow(z,1)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Pow(z, -1), 1.0 / z), "-1");
            }
        }


        [TestMethod]
        public void ComplexPowOneHalf () {
            Assert.IsTrue(ComplexMath.Pow(a, 0.5) == ComplexMath.Sqrt(a));
        }

        // log and exp

        [TestMethod]
        public void ComlexLogSpecialCaseTest () {
            Assert.IsTrue(ComplexMath.Log(Math.E) == 1);
            Assert.IsTrue(ComplexMath.Log(1.0) == 0.0);
            Assert.IsTrue(ComplexMath.Log(-1.0) == ComplexMath.I * Math.PI);
            Assert.IsTrue(ComplexMath.Log(ComplexMath.I) == ComplexMath.I * Math.PI / 2.0);
            Assert.IsTrue(ComplexMath.Log(-ComplexMath.I) == -ComplexMath.I * Math.PI / 2.0);
        }

        [TestMethod]
        public void ComplexLogExpTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 10)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Exp(ComplexMath.Log(z)), z));
            }
        }

        [TestMethod]
        public void ComplexLogSumTest () {
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
        public void ComplexExpSpecialCaseTest () {
            Assert.IsTrue(ComplexMath.Exp(0) == 1);
            Assert.IsTrue(ComplexMath.Exp(1) == Math.E);
        }

        [TestMethod]
        public void ComplexExpReflectionTest () {
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
        public void ComplexAbsTest () {
            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 10)) {
                double abs = ComplexMath.Abs(z);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(abs * abs, z * z.Conjugate));
            }
        }

        [TestMethod]
        public void PowComplexExponent () {

            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-2, 1.0E2, 6)) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Pow(Math.E, z), ComplexMath.Exp(z)));
            }

        }

        [TestMethod]
        public void ComplexPowSpecialCases () {

            foreach (Complex z in TestUtilities.GenerateComplexValues(1.0E-4, 1.0E4, 10)) {

                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Pow(z, 0), 1.0));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Pow(z, 1), z));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(ComplexMath.Pow(z, 2), z * z));

            }

        }

    }
}
