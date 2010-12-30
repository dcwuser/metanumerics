using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Test {


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

        [TestMethod]
        public void FindMaximumOfDawson () {

            Func<double, double> f = delegate(double x) {
                return (-AdvancedMath.Dawson(x));
            };
            LineExtremum maximum = FunctionMath.FindMinimum(f, Interval.FromEndpoints(0.0, 2.0));

            //Assert.IsTrue(TestUtilities.IsNearlyEqual(maximum.Location, 0.92413887300459176701));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(maximum.Value, -0.54104422463518169847));
        }

        [TestMethod]
        public void FindMinimumOfCosine () {

            Func<double, double> f = delegate(double x) {
                return (Math.Cos(x));
            };
            LineExtremum minimum = FunctionMath.FindMinimum(f, Interval.FromEndpoints(1.0, 6.0));
            Console.WriteLine("X = {0}", minimum.Location);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Location, Math.PI, Math.Sqrt(TestUtilities.TargetPrecision)));
            Console.WriteLine("Y = {0}", minimum.Value);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Value, -1.0));
            Console.WriteLine("d2Y/dX2 = {0}", minimum.Curvature);
            //Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Curvature, 1.0, Math.Sqrt(TestUtilities.TargetPrecision)));

        }

        [TestMethod]
        public void FindSpaceMinimumOfRosenbock () {

            Func<double[], double> f = delegate(double[] x) {
                double xx = x[0] * x[0];
                return( 100.0 * Math.Pow(x[1] - xx, 2) + Math.Pow(1.0-x[0],2) );
            };
            SpaceExtremum minimum = FunctionMath.FindMinimum(f, new double[] { 0, 0 });

            Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Location(), new double[] { 1, 1 }));

        }
    }
}
