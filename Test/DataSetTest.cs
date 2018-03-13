using System;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

// add tests of fit parameter uncertainties: are they really distributed as claimed?
// add T and KS tests against other samples

namespace Test {
    
    
    /// <summary>
    ///This is a test class for DataSetTest and is intended
    ///to contain all DataSetTest Unit Tests
    ///</summary>
    [TestClass()]
    public class DataSetTest {


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


        [TestMethod]
        public void DataSetManipulationsTest () {

            UncertainMeasurement<double> d1 = new UncertainMeasurement<double>(3.0, new UncertainValue(2.0, 1.0));
            UncertainMeasurement<double> d2 = new UncertainMeasurement<double>(-3.0, new UncertainValue(2.0, 1.0));
            UncertainMeasurement<double> d3 = new UncertainMeasurement<double>(3.0, new UncertainValue(-2.0, 1.0));

            Assert.IsTrue(d1 != null);

            UncertainMeasurement<double>[] data = new UncertainMeasurement<double>[] { d1, d2 };
            UncertainMeasurementSample set = new UncertainMeasurementSample();
            set.Add(data);

            Assert.IsFalse(set.Contains(d3));
            Assert.IsTrue(set.Count == data.Length);
            set.Add(d3);
            Assert.IsTrue(set.Contains(d3));
            Assert.IsTrue(set.Count == data.Length + 1);
            set.Remove(d3);
            Assert.IsFalse(set.Contains(d3));
            Assert.IsTrue(set.Count == data.Length);

            set.Clear();
            Assert.IsTrue(set.Count == 0);

        }

        // create a data set at n random points on the interval r based on the function fv
        // the function fu gives the uncertainty for each point
        private UncertainMeasurementSample CreateDataSet (Interval r, Func<double, double> fv, Func<double, double> fu, int n) {
            return (CreateDataSet(r, fv, fu, n, 1));
        }


        private UncertainMeasurementSample CreateDataSet (Interval r, Func<double, double> fv, Func<double, double> fu, int n, int seed) {
            UncertainMeasurementSample set = new UncertainMeasurementSample();

            UniformDistribution xd = new UniformDistribution(r);

            Random rng = new Random(seed);
            for (int i = 0; i < n; i++) {
                double x = xd.InverseLeftProbability(rng.NextDouble());
                double ym = fv(x);
                double ys = fu(x);
                NormalDistribution yd = new NormalDistribution(ym,ys);
                double y = yd.InverseLeftProbability(rng.NextDouble());

                //Console.WriteLine("{0}, {1}", x, new UncertainValue(y, ys));
                UncertainMeasurement<double> point = new UncertainMeasurement<double>(x, y, ys);
                set.Add(point);

            }

            return (set);
        }

        private UncertainMeasurementSample CreateDataSet (double[] xs, Func<double, double> fv, Func<double, double> fu, int seed) {
            UncertainMeasurementSample set = new UncertainMeasurementSample();
            Random rng = new Random(seed);
            foreach (double x in xs) {
                double ym = fv(x);
                double ys = fu(x);
                NormalDistribution yd = new NormalDistribution(ym, ys);
                double y = yd.InverseLeftProbability(rng.NextDouble());
                UncertainMeasurement<double> point = new UncertainMeasurement<double>(x, y, ys);
                set.Add(point);
            }
            return (set);
        }

        [TestMethod]
        public void FitDataToProportionalityTest () {
            Interval r = Interval.FromEndpoints(0.0, 0.1);
            Func<double,double> fv = delegate (double x) {
                return(0.5 * x);
            };
            Func<double, double> fu = delegate(double x) {
                return (0.02);
            };
            UncertainMeasurementSample set = CreateDataSet(r, fv, fu, 20);

            // fit to proportionality
            FitResult prop = set.FitToProportionality();
            Assert.IsTrue(prop.Dimension == 1);
            Assert.IsTrue(prop.Parameter(0).ConfidenceInterval(0.95).ClosedContains(0.5));
            Assert.IsTrue(prop.GoodnessOfFit.Probability > 0.05);

            // fit to line
            FitResult line = set.FitToLine();
            Assert.IsTrue(line.Dimension == 2);

            // line's intercept should be compatible with zero and slope with proportionality constant
            Assert.IsTrue(line.Parameter(0).ConfidenceInterval(0.95).ClosedContains(0.0));
            Assert.IsTrue(line.Parameter(1).ConfidenceInterval(0.95).ClosedContains(prop.Parameter(0).Value));

            // the fit should be better, but not too much better
            Assert.IsTrue(line.GoodnessOfFit.Statistic < prop.GoodnessOfFit.Statistic);

        }

        [TestMethod]
        public void FitDataToLineTest () {
            Interval r = Interval.FromEndpoints(0.0, 10.0);
            Func<double,double> fv = delegate (double x) {
                return(2.0 * x - 1.0);
            };
            Func<double, double> fu = delegate(double x) {
                return (1.0 + x);
            };
            UncertainMeasurementSample data = CreateDataSet(r, fv, fu, 20);


            // sanity check the data set
            Assert.IsTrue(data.Count == 20);

            // fit to a line
            FitResult line = data.FitToLine();
            Assert.IsTrue(line.Dimension == 2);
            Assert.IsTrue(line.Parameter(0).ConfidenceInterval(0.95).ClosedContains(-1.0));
            Assert.IsTrue(line.Parameter(1).ConfidenceInterval(0.95).ClosedContains(2.0));
            Assert.IsTrue(line.GoodnessOfFit.Probability > 0.05);

            // correlation coefficient should be related to covariance as expected
            Assert.IsTrue(TestUtilities.IsNearlyEqual(line.CorrelationCoefficient(0,1),line.Covariance(0,1)/line.Parameter(0).Uncertainty/line.Parameter(1).Uncertainty));

            // fit to a 1st order polynomial and make sure it agrees
            FitResult poly = data.FitToPolynomial(1);
            Assert.IsTrue(poly.Dimension == 2);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(poly.Parameters, line.Parameters));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(poly.CovarianceMatrix, line.CovarianceMatrix));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(poly.GoodnessOfFit.Statistic, line.GoodnessOfFit.Statistic));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(poly.GoodnessOfFit.Probability, line.GoodnessOfFit.Probability));

            // fit to a constant; the result should be poor
            FitResult constant = data.FitToConstant();
            Assert.IsTrue(constant.GoodnessOfFit.Probability < 0.05);

        }


        [TestMethod]
        public void FitDataToLineChiSquaredTest () {

            // we want to make sure that the chi^2 values we are producing from line fits are distributed as expected

            // create a sample to hold chi^2 values
            Sample chis = new Sample();

            // define a model
            Interval r = Interval.FromEndpoints(-5.0, 5.0);
            Func<double, double> fv = delegate(double x) {
                return (1.0 - 2.0 * x);
            };
            Func<double, double> fu = delegate(double x) {
                return (1.0 + 0.5 * Math.Sin(x));
            };

            // draw 50 data sets from the model and fit year
            // store the resulting chi^2 value in the chi^2 set
            for (int i = 0; i < 50; i++) {
                UncertainMeasurementSample xs = CreateDataSet(r, fv, fu, 10, i);
                FitResult fit = xs.FitToLine();
                double chi = fit.GoodnessOfFit.Statistic;
                chis.Add(chi);
            }

            // sanity check the sample
            Assert.IsTrue(chis.Count == 50);

            // test whether the chi^2 values are distributed as expected
            ContinuousDistribution chiDistribution = new ChiSquaredDistribution(8);
            TestResult ks = chis.KolmogorovSmirnovTest(chiDistribution);
            Assert.IsTrue(ks.Probability > 0.05);

        }

        [TestMethod]
        public void FitDataToPolynomialTest () {

            Interval r = Interval.FromEndpoints(-10.0, 10.0);
            Polynomial p = Polynomial.FromCoefficients(1.0, -2.0, 3.0, -4.0, 5.0, -6.0);
            Func<double, double> fu = delegate(double x) {
                return (1.0 + 0.5 * Math.Cos(x));
            };

            UncertainMeasurementSample set = CreateDataSet(r, p.Evaluate, fu, 50);
            Assert.IsTrue(set.Count == 50);

            // fit to an appropriate polynomial
            FitResult poly = set.FitToPolynomial(5);

            // the coefficients should match
            for (int i = 0; i < poly.Dimension; i++) {
                Assert.IsTrue(poly.Parameter(i).ConfidenceInterval(0.95).ClosedContains(p.Coefficient(i)));
            }

            // the fit should be good
            Assert.IsTrue(poly.GoodnessOfFit.Probability > 0.05);

            // fit to a lower order polynomial
            FitResult low = set.FitToPolynomial(4);

            // the fit should be bad
            Assert.IsTrue(low.GoodnessOfFit.Statistic > poly.GoodnessOfFit.Statistic);
            Assert.IsTrue(low.GoodnessOfFit.Probability < 0.05);

            // fit to a higher order polynomial
            FitResult high = set.FitToPolynomial(6);

            // the higher order coefficients should be compatible with zero
            Assert.IsTrue(high.Parameter(6).ConfidenceInterval(0.95).ClosedContains(0.0));

            // the fit should be better, but not too much better
            Assert.IsTrue(high.GoodnessOfFit.Statistic < poly.GoodnessOfFit.Statistic);

        }

        [TestMethod]
        public void FitDataToPolynomialChiSquaredTest () {

            // we want to make sure that the chi^2 values we are producing from polynomial fits are distributed as expected

            // create a sample to hold chi^2 values
            Sample chis = new Sample();

            // define a model
            Interval r = Interval.FromEndpoints(-5.0,15.0);
            Func<double,double> fv = delegate(double x) {
                return(1.0*x - 2.0*x*x);
            };
            Func<double,double> fu = delegate(double x) {
                return(1.0 + 0.5 * Math.Sin(x));
            };

            // draw 50 data sets from the model and fit year
            // store the resulting chi^2 value in the chi^2 set
            for (int i = 0; i < 50; i++) {
                UncertainMeasurementSample xs = CreateDataSet(r, fv, fu, 10, i);
                FitResult fit = xs.FitToPolynomial(2);
                double chi = fit.GoodnessOfFit.Statistic;
                chis.Add(chi);
            }

            // sanity check the sample
            Assert.IsTrue(chis.Count == 50);

            // test whether the chi^2 values are distributed as expected
            ContinuousDistribution chiDistribution = new ChiSquaredDistribution(7);
            TestResult ks = chis.KolmogorovSmirnovTest(chiDistribution);
            Assert.IsTrue(ks.Probability > 0.05);

        }

        [TestMethod]
        public void FitDataToPolynomialUncertaintiesTest () {

            // make sure the reported uncertainties it fit parameters really represent their standard deviation,
            // and that the reported off-diagonal elements really represent their correlations

            double[] xs = TestUtilities.GenerateUniformRealValues(-1.0, 2.0, 10);
            Func<double, double> fv = delegate(double x) {
                return (0.0 + 1.0 * x + 2.0 * x * x);
            };
            Func<double, double> fu = delegate(double x) {
                return (0.5);
            };

            // keep track of best-fit parameters and claimed parameter covariances
            MultivariateSample sample = new MultivariateSample(3);

            // generate 50 small data sets and fit each
            FitResult[] fits = new FitResult[50];
            for (int i = 0; i < fits.Length; i++) {
                //DataSet set = CreateDataSet(Interval.FromEndpoints(-1.0,2.0), fv, fu, 10, i);
                UncertainMeasurementSample set = CreateDataSet(xs, fv, fu, 314159+i);
                /*
                foreach (DataPoint point in set) {
                    Console.WriteLine("  i={0} x={1} y={2}", i, point.X, point.Y);
                }
                */
                fits[i] = set.FitToPolynomial(2);
                /*
                for (int j = 0; j < fits[i].Dimension; j++) {
                    Console.WriteLine("p[{0}] = {1}", j, fits[i].Parameter(j));
                }
                */

                for (int j = 0; j < fits[i].Dimension; j++) {
                    sample.Add(fits[i].Parameters);
                }


            }

            // check that parameters agree
            for (int i = 0; i < 3; i++) {
                Console.WriteLine(sample.Column(i).PopulationMean);
            }
            //Assert.IsTrue(sample.PopulationMean(0).ConfidenceInterval(0.95).ClosedContains(0.0));
            //Assert.IsTrue(sample.PopulationMean(1).ConfidenceInterval(0.95).ClosedContains(1.0));
            //Assert.IsTrue(sample.PopulationMean(2).ConfidenceInterval(0.95).ClosedContains(2.0));

            /*
            // check that the claimed covariances agree with the measured covariances
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {

                    Console.WriteLine("{0},{1} {2} {3}", i, j, sample.PopulationCovariance(i, j), C[i, j]);

                    Assert.IsTrue(sample.PopulationCovariance(i, j).ConfidenceInterval(0.95).ClosedContains(C[i, j]));
                }
            }
            */

            // for each parameter, verify that the standard devition of the reported values agrees with the (average) reported uncertainty
            double[] pMeans = new double[3];
            for (int i = 0; i <= 2; i++) {
                Console.WriteLine("paramater {0}", i);
                Sample values = new Sample();
                Sample uncertainties = new Sample();
                for (int j = 0; j < fits.Length; j++) {
                    UncertainValue p = fits[j].Parameter(i);
                    values.Add(p.Value);
                    uncertainties.Add(p.Uncertainty);
                }
                pMeans[i] = values.Mean;
                Console.WriteLine("reported values mean={0}, standard deviation={1}", values.Mean, values.StandardDeviation);
                Console.WriteLine("reported uncertainties mean={0}", uncertainties.Mean);
                Assert.IsTrue(values.PopulationStandardDeviation.ConfidenceInterval(0.95).ClosedContains(uncertainties.Mean));
            }

            // for each parameter pair, verify that the covariances of the reported values agrees with the (average) reported covarance
            /*
            for (int i = 0; i <= 2; i++) {
                for (int j = 0; j <= 2; j++) {
                    // compute cov(i,j)
                    double cov = 0.0;
                    for (int k = 0; k < fits.Length; k++) {
                        cov += (fits[k].Parameter(i).Value - pMeans[i]) * (fits[k].Parameter(j).Value - pMeans[j]);
                    }
                    cov = cov / fits.Length;
                    // collect the reported covarainces
                    Sample covs = new Sample();
                    for (int k = 0; k < fits.Length; k++) {
                        covs.Add(fits[k].Covariance(i, j));
                    }
                    Console.WriteLine("cov({0},{1}) computed={2} reported={3}", i, j, cov, covs.PopulationMean);
                    // the computed covariance should agree with the (average) reported covariance
                    // problem: we need to estimate the uncertainty in our covariance, but that isn't ever computed
                    // note: covs just depend on x's, so the reported cov is the same for all fits
                    // solution: for now, just assume it scales with 1/Sqrt(N); long term, do a multivariate fit
                    //Assert.IsTrue((new UncertainValue(cov, Math.Abs(cov)* Math.Sqrt(2.0/fits.Length))).ConfidenceInterval(0.95).ClosedContains(covs.Mean));
                }
            }
            */

        }

        [TestMethod]
        public void FitDataToLinearFunctionTest () {

            // create a data set from a linear combination of sine and cosine
            Interval r = Interval.FromEndpoints(-4.0, 6.0);
            double[] c = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
            Func<double, double> fv = delegate(double x) {
                return (2.0 * Math.Cos(x) + 1.0 * Math.Sin(x));
            };
            Func<double, double> fu = delegate(double x) {
                return (0.1 + 0.1 * Math.Abs(x));
            };
            UncertainMeasurementSample set = CreateDataSet(r, fv, fu, 20, 2);

            // fit the data set to a linear combination of sine and cosine
            Func<double,double>[] fs = new Func<double,double>[]
                { delegate (double x) { return(Math.Cos(x)); }, delegate (double x) { return(Math.Sin(x)); } };
            FitResult result = set.FitToLinearFunction(fs);

            // the fit should be right right dimension
            Assert.IsTrue(result.Dimension == 2);

            // the coefficients should match
            Console.WriteLine(result.Parameter(0));
            Console.WriteLine(result.Parameter(1));
            Assert.IsTrue(result.Parameter(0).ConfidenceInterval(0.95).ClosedContains(2.0));
            Assert.IsTrue(result.Parameter(1).ConfidenceInterval(0.95).ClosedContains(1.0));

            // diagonal covarainces should match errors
            Assert.IsTrue(TestUtilities.IsNearlyEqual(Math.Sqrt(result.Covariance(0,0)), result.Parameter(0).Uncertainty));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(Math.Sqrt(result.Covariance(1,1)), result.Parameter(1).Uncertainty));

        }

        [TestMethod]
        public void FitDataToFunctionTest () {

            // create a data set from a nonlinear function
            /*
            Interval r = Interval.FromEndpoints(-3.0, 5.0);
            double[] c = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
            Function<double, double> fv = delegate(double x) {
                return (3.0 * Math.Cos(2.0 * Math.PI * x / 2.0 - 1.0));
            };
            Function<double, double> fu = delegate(double x) {
                return (0.1 + 0.1 * Math.Abs(x));
            };
            DataSet set = CreateDataSet(r, fv, fu, 20, 2);
            */

            UncertainMeasurementSample set = new UncertainMeasurementSample();
            set.Add(new UncertainMeasurement<double>(1.0, 1.0, 0.1));
            set.Add(new UncertainMeasurement<double>(2.0, 0.7, 0.1));
            set.Add(new UncertainMeasurement<double>(3.0, 0.0, 0.1));
            set.Add(new UncertainMeasurement<double>(4.0, -0.7, 0.1));
            set.Add(new UncertainMeasurement<double>(5.0, -1.0, 0.1));
            set.Add(new UncertainMeasurement<double>(6.0, -0.7, 0.1));
            set.Add(new UncertainMeasurement<double>(7.0, 0.0, 0.1));
            set.Add(new UncertainMeasurement<double>(8.0, 0.7, 0.1));
            set.Add(new UncertainMeasurement<double>(9.0, 1.0, 0.1));

            // fit it to a parameterized fit function
            /*
            Function<double[], double, double> ff = delegate(double[] p, double x) {
                return (p[0] * Math.Cos(2.0 * Math.PI / p[1] + p[2]));
            };
            */
            Func<double[], double, double> ff = delegate(double[] p, double x) {
                //Console.WriteLine("    p[0]={0}, x={1}", p[0], x);
                return (p[1] * Math.Cos(x / p[0] + p[2]));
                //return (x / p[0]);
            };
            FitResult fit = set.FitToFunction(ff, new double[] { 1.3, 1.1, 0.1 });

            Console.WriteLine(fit.Parameter(0));
            Console.WriteLine(fit.Parameter(1));
            Console.WriteLine(fit.Parameter(2));

        }

        [TestMethod]
        public void FitToFunctionLinearCompatibilityTest () {

            // specify a cubic function
            Interval r = Interval.FromEndpoints(-5.0, 5.0);
            Func<double, double> fv = delegate(double x) {
                return (1.0 + 2.0 * x);
                //return (0.0 - 1.0 * x + 2.0 * x * x - 3.0 * x * x * x);
            };
            Func<double, double> fu = delegate(double x) {
                return (1.0 + 0.5 * Math.Sin(x));
            };

            // create a data set from it
            UncertainMeasurementSample set = CreateDataSet(r, fv, fu, 30);

            // fit it to a cubic polynomial
            FitResult pFit = set.FitToLine();
            //FitResult pFit = set.FitToPolynomial(3);

            // fit it to a cubic polynomaial
            Func<double[], double, double> ff = delegate(double[] p, double x) {
                return (p[0] + p[1] * x);
                //return (p[0] + p[1] * x + p[2] * x * x + p[3] * x * x * x);
            };
            FitResult fFit = set.FitToFunction(ff, new double[] { 0, 0});

            // the fits should agree
            for (int i = 0; i < pFit.Dimension; i++) {
                Console.WriteLine("{0} ?= {1}", pFit.Parameter(i), fFit.Parameter(i));
            }
            // dimension
            Assert.IsTrue(pFit.Dimension == fFit.Dimension);
            // chi squared
            Assert.IsTrue(TestUtilities.IsNearlyEqual(pFit.GoodnessOfFit.Statistic, fFit.GoodnessOfFit.Statistic, Math.Sqrt(TestUtilities.TargetPrecision)));
            // parameters
            Assert.IsTrue(TestUtilities.IsNearlyEqual(pFit.Parameters, fFit.Parameters, Math.Sqrt(TestUtilities.TargetPrecision)));
            // covariance
            Assert.IsTrue(TestUtilities.IsNearlyEqual(pFit.CovarianceMatrix, fFit.CovarianceMatrix, Math.Sqrt(TestUtilities.TargetPrecision)));

        }

        [TestMethod]
        public void FitToFunctionPolynomialCompatibilityTest () {

            // specify a cubic function
            Interval r = Interval.FromEndpoints(-10.0, 10.0);
            Func<double, double> fv = delegate(double x) {
                return (0.0 - 1.0 * x + 2.0 * x * x - 3.0 * x * x * x);
            };
            Func<double, double> fu = delegate(double x) {
                return (1.0 + 0.5 * Math.Cos(x));
            };

            // create a data set from it
            UncertainMeasurementSample set = CreateDataSet(r, fv, fu, 60);

            // fit it to a cubic polynomial
            FitResult pFit = set.FitToPolynomial(3);

            // fit it to a cubic polynomaial
            Func<double[], double, double> ff = delegate(double[] p, double x) {
                return (p[0] + p[1] * x + p[2] * x * x + p[3] * x * x * x);
            };
            FitResult fFit = set.FitToFunction(ff, new double[] { 0, 0, 0, 0 });

            // the fits should agree
            Console.WriteLine("{0} ?= {1}", pFit.GoodnessOfFit.Statistic, fFit.GoodnessOfFit.Statistic);
            for (int i = 0; i < pFit.Dimension; i++) {
                Console.WriteLine("{0} ?= {1}", pFit.Parameter(i), fFit.Parameter(i));
                Assert.IsTrue(pFit.Parameter(i).ConfidenceInterval(0.01).ClosedContains(fFit.Parameter(i).Value));
            }
            // dimension
            Assert.IsTrue(pFit.Dimension == fFit.Dimension);
            // chi squared
            Assert.IsTrue(TestUtilities.IsNearlyEqual(pFit.GoodnessOfFit.Statistic, fFit.GoodnessOfFit.Statistic, Math.Sqrt(TestUtilities.TargetPrecision)));
            // don't demand super-high precision agreement of parameters and covariance matrix
            // parameters
            Assert.IsTrue(TestUtilities.IsNearlyEqual(pFit.Parameters, fFit.Parameters, Math.Pow(TestUtilities.TargetPrecision,0.3)));
            // covariance
            Assert.IsTrue(TestUtilities.IsNearlyEqual(pFit.CovarianceMatrix, fFit.CovarianceMatrix, Math.Pow(TestUtilities.TargetPrecision,0.3)));

        }

        [TestMethod]
        public void FitDataToLineUncertaintyTest () {

            double[] xs = TestUtilities.GenerateUniformRealValues(0.0, 10.0, 10);
            Func<double,double> fv = delegate (double x) {
                return(2.0*x - 1.0);
            };
            Func<double, double> fu = delegate(double x) {
                return (1.0+x);
            };

            MultivariateSample sample = new MultivariateSample(2);
            SymmetricMatrix covariance = new SymmetricMatrix(2);

            // create a bunch of small data sets
            for (int i = 0; i < 100; i++) {
                UncertainMeasurementSample data = CreateDataSet(xs, fv, fu, i);
                FitResult fit = data.FitToLine();

                sample.Add(fit.Parameters);
                covariance = fit.CovarianceMatrix;
                // because it depends only on the x's and sigmas, the covariance is always the same

                Console.WriteLine("cov_00 = {0}", covariance[0, 0]);
            }

            // the measured covariances should agree with the claimed covariances
            //Assert.IsTrue(sample.PopulationCovariance(0,0).ConfidenceInterval(0.95).ClosedContains(covariance[0,0]));
            //Assert.IsTrue(sample.PopulationCovariance(0,1).ConfidenceInterval(0.95).ClosedContains(covariance[0,1]));
            //Assert.IsTrue(sample.PopulationCovariance(1,0).ConfidenceInterval(0.95).ClosedContains(covariance[1,0]));
            //Assert.IsTrue(sample.PopulationCovariance(1,1).ConfidenceInterval(0.95).ClosedContains(covariance[1,1]));

        }
    
    }

}
