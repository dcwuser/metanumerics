using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Test {

    internal class TestExtremum {

        public TestExtremum (Func<double, double> function, Interval interval, double location, double value, double curvature) {
            this.Function = function;
            this.Interval = interval;
            this.Location = location;
            this.Value = value;
            this.Curvature = curvature;
        }

        public Func<double, double> Function { get; set; }
        public Interval Interval { get; set; }

        public double Location { get; set; }
        public double Value { get; set; }
        public double Curvature { get; set; }

        public bool Agrees (LineExtremum extremum) {
            return (
                TestUtilities.IsNearlyEqual(extremum.Value, Value, 10.0 * TestUtilities.TargetPrecision) &&
                TestUtilities.IsNearlyEqual(extremum.Location, Location, 10.0 * Math.Sqrt(TestUtilities.TargetPrecision)) &&
                (Double.IsNaN(Curvature) || TestUtilities.IsNearlyEqual(extremum.Curvature, Curvature, 0.01))
            );
        }

    }
    

    [TestClass()]
    public class ExtremaTest {


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

        private TestExtremum[] testMinima = new TestExtremum[] {
            new TestExtremum(x => x * Math.Log(x), Interval.FromEndpoints(0.1, 1.0), 1.0 / Math.E, -1.0 / Math.E, Math.E),
            new TestExtremum(x => Math.Cos(x), Interval.FromEndpoints(0.0, 5.0), Math.PI, -1.0, 1.0),
            new TestExtremum(AdvancedMath.Gamma, Interval.FromEndpoints(0.0, 2.0), 1.461632144968362, 0.885603194410889, Double.NaN),
            new TestExtremum(Math.Abs, Interval.FromEndpoints(-2.0, 3.0), 0.0, 0.0, Double.NaN),
            new TestExtremum(Math.Cosh, Interval.FromEndpoints(-1.5, 2.5), 0.0, 1.0, 1.0),
            new TestExtremum(x => MoreMath.Pow(x, 4), Interval.FromEndpoints(-2.0, 1.0), 0.0, 0.0, Double.NaN)
        };

        private TestExtremum[] testMaxima = new TestExtremum[] {
            new TestExtremum(Math.Sin, Interval.FromEndpoints(0.0, 5.0), Math.PI / 2.0, 1.0, 1.0),
            new TestExtremum(AdvancedMath.Dawson, Interval.FromEndpoints(0.0, 1.0), 0.92413887300459176701, 0.54104422463518169847, Double.NaN)
        };

        [TestMethod]
        public void FindMinimaFromInterval () {
            foreach (TestExtremum testMinimum in testMinima) {
                LineExtremum minimum = FunctionMath.FindMinimum(testMinimum.Function, testMinimum.Interval);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Location, testMinimum.Location, Math.Sqrt(TestUtilities.TargetPrecision)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Value, testMinimum.Value, TestUtilities.TargetPrecision));
                //if (!Double.IsNaN(testMinimum.Curvature)) Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Curvature, testMinimum.Curvature, 0.1));
            }
        }

        [TestMethod]
        public void FindMinimaFromPoint () {
            foreach (TestExtremum testExtremum in testMinima) {
                LineExtremum extremum = FunctionMath.FindMinimum(testExtremum.Function, testExtremum.Interval.Midpoint);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(extremum.Location, testExtremum.Location, Math.Sqrt(TestUtilities.TargetPrecision)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(extremum.Value, testExtremum.Value, TestUtilities.TargetPrecision));
                //if (!Double.IsNaN(testMinimum.Curvature)) Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Curvature, testMinimum.Curvature, 0.05));
             }
        }

        [TestMethod]
        public void FindMaximaFromInterval () {
            foreach (TestExtremum testExtremum in testMaxima) {
                LineExtremum extremum = FunctionMath.FindMaximum(testExtremum.Function, testExtremum.Interval);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(extremum.Location, testExtremum.Location, Math.Sqrt(TestUtilities.TargetPrecision)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(extremum.Value, testExtremum.Value, TestUtilities.TargetPrecision));
            }
        }

        [TestMethod]
        public void FindMaximaFromPoint () {
            foreach (TestExtremum testExtremum in testMaxima) {
                LineExtremum extremum = FunctionMath.FindMaximum(testExtremum.Function, testExtremum.Interval.Midpoint);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(extremum.Location, testExtremum.Location, Math.Sqrt(TestUtilities.TargetPrecision)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(extremum.Value, testExtremum.Value, TestUtilities.TargetPrecision));
            }
        }

        [TestMethod]
        public void FindMinimumOfGamma () {

            Func<double, double> f = new Func<double, double>(AdvancedMath.Gamma);
            LineExtremum minimum = FunctionMath.FindMinimum(f, 1.5);

            Console.WriteLine(minimum.Location);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Location, 1.46163214496836234126, Math.Sqrt(TestUtilities.TargetPrecision)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Value, 0.88560319441088870028));
            Console.WriteLine(minimum.Curvature);
            // add test for known curvature 0.8569736317111708

            // test conversion to space extremum
            SpaceExtremum space = minimum;
            Assert.IsTrue(space.Dimension == 1);
            double[] location = space.Location();
            Assert.IsTrue(location.Length == 1);
            Assert.IsTrue(location[0] == minimum.Location);
            Assert.IsTrue(space.Value == minimum.Value);
            SymmetricMatrix curvature = space.Curvature();
            Assert.IsTrue(curvature.Dimension == 1);
            Assert.IsTrue(curvature[0, 0] == minimum.Curvature);

            
        }

        // Some standard minimization test functions: http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf

        [TestMethod]
        public void FindSpaceMinimumOfRosenbock () {

            int n = 4; {
            //foreach (int n in TestUtilities.GenerateIntegerValues(2, 16, 4)) {

                Func<double[], double> f = delegate(double[] x) {
                    double s = 0.0;
                    for (int i = 0; i < (n - 1); i++) {
                        s += 100.0 * MoreMath.Pow(x[i + 1] - x[i] * x[i], 2) + MoreMath.Pow(1.0 - x[i], 2);
                    }
                    return (s);
                };

                double[] x0 = new double[n];

                EvaluationSettings settings = new EvaluationSettings() { EvaluationBudget = 200000, AbsolutePrecision = TestUtilities.TargetPrecision * TestUtilities.TargetPrecision, RelativePrecision = TestUtilities.TargetPrecision };

                SpaceExtremum minimum = FunctionMath.FindMinimum(f, x0, settings);

                double[] x1 = new double[n]; for (int i = 0; i < n; i++) x1[i] = 1.0;


                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Location(), x1, Math.Sqrt(TestUtilities.TargetPrecision)));

            }

        }

        [TestMethod]
        public void MinimizeGoldsteinPrice () {

            double[] p0 = new double[] { -0.2, -0.8 };
            SpaceExtremum m = FunctionMath.FindMinimum(
                (double[] p) => {
                    double x = p[0]; double y = p[1]; return (
                        (1 + MoreMath.Pow(x + y + 1, 2) * (19 - 14 * x + 3 * x * x - 14 * y + 6 * x * y + 6 * y * y)) *
                        (30 + MoreMath.Pow(2 * x - 3 * y, 2) * (18 - 32 * x + 12 * x * x + 48 * y - 36 * x * y + 27 * y * y))
                    );
                }, p0
            );

            double[] p1 = m.Location();
            Console.WriteLine("{0} {1} {2}", p1[0], p1[1], m.Value);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(p1[0], 0.0, Math.Sqrt(TestUtilities.TargetPrecision)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(p1[1], -1.0, Math.Sqrt(TestUtilities.TargetPrecision)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(m.Value, 3.0, Math.Sqrt(TestUtilities.TargetPrecision)));
        }
    }
}
