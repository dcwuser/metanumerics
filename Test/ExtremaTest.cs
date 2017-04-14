using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
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

        public bool Agrees (Extremum extremum) {
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
                Extremum minimum = FunctionMath.FindMinimum(testMinimum.Function, testMinimum.Interval);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Location, testMinimum.Location, Math.Sqrt(TestUtilities.TargetPrecision)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Value, testMinimum.Value, TestUtilities.TargetPrecision));
                //if (!Double.IsNaN(testMinimum.Curvature)) Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Curvature, testMinimum.Curvature, 0.1));
            }
        }

        [TestMethod]
        public void FindMinimaFromPoint () {
            foreach (TestExtremum testExtremum in testMinima) {
                Extremum extremum = FunctionMath.FindMinimum(testExtremum.Function, testExtremum.Interval.Midpoint);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(extremum.Location, testExtremum.Location, Math.Sqrt(TestUtilities.TargetPrecision)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(extremum.Value, testExtremum.Value, TestUtilities.TargetPrecision));
                //if (!Double.IsNaN(testMinimum.Curvature)) Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Curvature, testMinimum.Curvature, 0.05));
             }
        }

        [TestMethod]
        public void FindMaximaFromInterval () {
            foreach (TestExtremum testExtremum in testMaxima) {
                Extremum extremum = FunctionMath.FindMaximum(testExtremum.Function, testExtremum.Interval);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(extremum.Location, testExtremum.Location, Math.Sqrt(TestUtilities.TargetPrecision)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(extremum.Value, testExtremum.Value, TestUtilities.TargetPrecision));
            }
        }

        [TestMethod]
        public void FindMaximaFromPoint () {
            foreach (TestExtremum testExtremum in testMaxima) {
                Extremum extremum = FunctionMath.FindMaximum(testExtremum.Function, testExtremum.Interval.Midpoint);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(extremum.Location, testExtremum.Location, Math.Sqrt(TestUtilities.TargetPrecision)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(extremum.Value, testExtremum.Value, TestUtilities.TargetPrecision));
            }
        }

        [TestMethod]
        public void FindMinimumOfGamma () {

            int count = 0;
            ExtremumSettings settings = new ExtremumSettings() { Listener = r => { count++; } };

            Func<double, double> f = new Func<double, double>(AdvancedMath.Gamma);
            Extremum minimum = FunctionMath.FindMinimum(f, 1.5, settings);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Location, 1.46163214496836234126, Math.Sqrt(TestUtilities.TargetPrecision)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Value, 0.88560319441088870028));
            //Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Curvature, Math.Sqrt(Math.Sqrt(0.8569736317111708))));
            Assert.IsTrue(count > 0);
            
        }

        [TestMethod]
        public void FindExtremaNegativeGamma () {

            // https://en.wikipedia.org/wiki/Particular_values_of_the_Gamma_function#Other_constants

            Tuple<double, double>[] extrema = new Tuple<double, double>[] {
                Tuple.Create(-0.5040830082644554092582693045, -3.5446436111550050891219639933),
                Tuple.Create(-1.5734984731623904587782860437, 2.3024072583396801358235820396),
                Tuple.Create(-2.6107208684441446500015377157, -0.8881363584012419200955280294),
                Tuple.Create(-3.6352933664369010978391815669, 0.2451275398343662504382300889)
            };

            foreach(Tuple<double, double> extremum in extrema) {

                // We use brackets since we know the extremum must lie between the singularities.
                // We should be able to use the actual singularities at endpoints, but this doesn't work; look into it.

                Interval bracket = Interval.FromEndpoints(Math.Floor(extremum.Item1) + 0.01, Math.Ceiling(extremum.Item1) - 0.01);
                Extremum result;
                if (extremum.Item2 < 0.0) {
                    result = FunctionMath.FindMaximum(AdvancedMath.Gamma, bracket);
                } else {
                    result = FunctionMath.FindMinimum(AdvancedMath.Gamma, bracket);
                }
                Assert.IsTrue(result.Bracket.OpenContains(extremum.Item1));
                Assert.IsTrue(result.Bracket.OpenContains(result.Location));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(result.Value, extremum.Item2));

            }

        }

        // Some standard minimization test functions: http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf

        /*
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
        */

    }
}
