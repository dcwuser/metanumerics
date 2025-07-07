using System;
using System.Collections.Generic;

using TestClassAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestClassAttribute;
using TestMethodAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestMethodAttribute;
using ExpectedExceptionAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.ExpectedExceptionAttribute;
using Assert = Microsoft.VisualStudio.TestTools.UnitTesting.Assert;

using Meta.Numerics;
using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

// add tests of fit parameter uncertainties: are they really distributed as claimed?
// add T and KS tests against other samples

namespace Test {
        
    [TestClass]
    public class UncertainMeasurementSampleTests {


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
            UncertainMeasurementFitResult prop = set.FitToProportionality();
            Assert.IsTrue(prop.Parameters.Count == 1);
            Assert.IsTrue(prop.Parameters[0].Estimate.ConfidenceInterval(0.95).ClosedContains(0.5));
            Assert.IsTrue(prop.GoodnessOfFit.Probability > 0.05);

            // fit to line
            UncertainMeasurementFitResult line = set.FitToLine();
            Assert.IsTrue(line.Parameters.Count == 2);

            // line's intercept should be compatible with zero and slope with proportionality constant
            Assert.IsTrue(line.Parameters[0].Estimate.ConfidenceInterval(0.95).ClosedContains(0.0));
            Assert.IsTrue(line.Parameters[1].Estimate.ConfidenceInterval(0.95).ClosedContains(prop.Parameters[0].Estimate.Value));

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
            UncertainMeasurementFitResult line = data.FitToLine();
            Assert.IsTrue(line.Parameters.Count == 2);
            Assert.IsTrue(line.Parameters[0].Estimate.ConfidenceInterval(0.95).ClosedContains(-1.0));
            Assert.IsTrue(line.Parameters[1].Estimate.ConfidenceInterval(0.95).ClosedContains(2.0));
            Assert.IsTrue(line.GoodnessOfFit.Probability > 0.05);

            // correlation coefficient should be related to covariance as expected
            //Assert.IsTrue(TestUtilities.IsNearlyEqual(line.Parameters.CorrelationCoefficient(0,1),line.Covariance(0,1)/line.Parameter(0).Uncertainty/line.Parameter(1).Uncertainty));

            // fit to a 1st order polynomial and make sure it agrees
            UncertainMeasurementFitResult poly = data.FitToPolynomial(1);
            Assert.IsTrue(poly.Parameters.Count == 2);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(poly.Parameters.ValuesVector, line.Parameters.ValuesVector));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(poly.Parameters.CovarianceMatrix, line.Parameters.CovarianceMatrix));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(poly.GoodnessOfFit.Statistic, line.GoodnessOfFit.Statistic));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(poly.GoodnessOfFit.Probability, line.GoodnessOfFit.Probability));

            // fit to a constant; the result should be poor
            UncertainMeasurementFitResult constant = data.FitToConstant();
            Assert.IsTrue(constant.GoodnessOfFit.Probability < 0.05);

        }


        [TestMethod]
        public void FitDataToLineChiSquaredTest () {

            // we want to make sure that the chi^2 values we are producing from line fits are distributed as expected

            // create a sample to hold chi^2 values
            List<double> chis = new List<double>();

            // define a model
            Interval r = Interval.FromEndpoints(-5.0, 5.0);
            Func<double, double> fv = (double x) => 1.0 - 2.0 * x;
            Func<double, double> fu = (double x) => 1.0 + 0.5 * Math.Sin(x);

            // draw 50 data sets from the model and fit year
            // store the resulting chi^2 value in the chi^2 set
            ContinuousDistribution chiSquaredDistribution = null;
            for (int i = 0; i < 50; i++) {
                UncertainMeasurementSample xs = CreateDataSet(r, fv, fu, 10, i);
                UncertainMeasurementFitResult fit = xs.FitToLine();
                double chi = fit.GoodnessOfFit.Statistic;
                chis.Add(chi);
                chiSquaredDistribution = fit.GoodnessOfFit.Statistic.Distribution;
            }

            // sanity check the sample
            Assert.IsTrue(chis.Count == 50);

            // test whether the chi^2 values are distributed as expected
            TestResult ks = chis.KolmogorovSmirnovTest(chiSquaredDistribution);
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
            UncertainMeasurementFitResult poly = set.FitToPolynomial(5);

            // the coefficients should match
            for (int i = 0; i < poly.Parameters.Count; i++) {
                Assert.IsTrue(poly.Parameters[i].Estimate.ConfidenceInterval(0.95).ClosedContains(p.Coefficient(i)));
            }

            // the fit should be good
            Assert.IsTrue(poly.GoodnessOfFit.Probability > 0.05);

            // fit to a lower order polynomial
            UncertainMeasurementFitResult low = set.FitToPolynomial(4);

            // the fit should be bad
            Assert.IsTrue(low.GoodnessOfFit.Statistic > poly.GoodnessOfFit.Statistic);
            Assert.IsTrue(low.GoodnessOfFit.Probability < 0.05);

            // fit to a higher order polynomial
            UncertainMeasurementFitResult high = set.FitToPolynomial(6);

            // the higher order coefficients should be compatible with zero
            Assert.IsTrue(high.Parameters[6].Estimate.ConfidenceInterval(0.95).ClosedContains(0.0));

            // the fit should be better, but not too much better
            Assert.IsTrue(high.GoodnessOfFit.Statistic < poly.GoodnessOfFit.Statistic);

        }

        [TestMethod]
        public void FitDataToPolynomialChiSquaredTest () {

            // we want to make sure that the chi^2 values we are producing from polynomial fits are distributed as expected

            // create a sample to hold chi^2 values
            List<double> chis = new List<double>();

            // define a model
            Interval r = Interval.FromEndpoints(-5.0,15.0);
            Func<double,double> fv = (double x) => 1.0 * x - 2.0 * x * x;
            Func<double,double> fu = (double x) => 1.0 + 0.5 * Math.Sin(x);

            // draw 50 data sets from the model and fit year
            // store the resulting chi^2 value in the chi^2 set
            for (int i = 0; i < 50; i++) {
                UncertainMeasurementSample xs = CreateDataSet(r, fv, fu, 10, i);
                UncertainMeasurementFitResult fit = xs.FitToPolynomial(2);
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

            // make sure the reported uncertainties in fit parameters really represent their standard deviation,
            // and that the reported off-diagonal elements really represent their correlations

            double[] xs = TestUtilities.GenerateUniformRealValues(-1.0, 2.0, 10);
            Func<double, double> fv = (double x) => 0.0 + 1.0 * x + 2.0 * x * x;
            Func<double, double> fu = (double x) => 0.5;

            // keep track of best-fit parameters and claimed parameter covariances
            //MultivariateSample sample = new MultivariateSample(3);

            // generate 50 small data sets and fit each
            UncertainMeasurementFitResult[] fits = new UncertainMeasurementFitResult[50];
            for (int i = 0; i < fits.Length; i++) {
                UncertainMeasurementSample set = CreateDataSet(xs, fv, fu, 314159+i);
                fits[i] = set.FitToPolynomial(2);
                //sample.Add(fits[i].Parameters.ValuesVector);
            }

            // check that parameters agree
            //for (int i = 0; i < 3; i++) {
            //    Console.WriteLine(sample.Column(i).PopulationMean);
            //}
 
            // for each parameter, verify that the standard deviation of the reported values agrees with the (average) reported uncertainty
            double[] pMeans = new double[3];
            for (int i = 0; i <= 2; i++) {
                List<double> values = new List<double>();
                List<double> uncertainties = new List<double>();
                for (int j = 0; j < fits.Length; j++) {
                    UncertainValue p = fits[j].Parameters[i].Estimate;
                    values.Add(p.Value);
                    uncertainties.Add(p.Uncertainty);
                }
                pMeans[i] = values.Mean();
                Assert.IsTrue(values.PopulationStandardDeviation().ConfidenceInterval(0.95).Contains(uncertainties.Mean()));
            }
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
            UncertainMeasurementFitResult result = set.FitToLinearFunction(fs);

            // the fit should be right right dimension
            Assert.IsTrue(result.Parameters.Count == 2);

            // the coefficients should match
            Assert.IsTrue(result.Parameters[0].Estimate.ConfidenceInterval(0.95).ClosedContains(2.0));
            Assert.IsTrue(result.Parameters[1].Estimate.ConfidenceInterval(0.95).ClosedContains(1.0));

            // diagonal covarainces should match errors
            Assert.IsTrue(TestUtilities.IsNearlyEqual(Math.Sqrt(result.Parameters.CovarianceMatrix[0,0]), result.Parameters[0].Estimate.Uncertainty));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(Math.Sqrt(result.Parameters.CovarianceMatrix[1,1]), result.Parameters[1].Estimate.Uncertainty));

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
            UncertainMeasurementFitResult fit = set.FitToFunction(ff, new double[] { 1.3, 1.1, 0.1 });
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
            UncertainMeasurementFitResult pFit = set.FitToLine();
            //FitResult pFit = set.FitToPolynomial(3);

            // fit it to a cubic polynomial
            Func<double[], double, double> ff = delegate(double[] p, double x) {
                return (p[0] + p[1] * x);
                //return (p[0] + p[1] * x + p[2] * x * x + p[3] * x * x * x);
            };
            UncertainMeasurementFitResult fFit = set.FitToFunction(ff, new double[] { 0, 0});

            // dimension
            Assert.IsTrue(pFit.Parameters.Count == fFit.Parameters.Count);
            // chi squared
            Assert.IsTrue(TestUtilities.IsNearlyEqual(pFit.GoodnessOfFit.Statistic, fFit.GoodnessOfFit.Statistic, Math.Sqrt(TestUtilities.TargetPrecision)));
            // parameters
            Assert.IsTrue(TestUtilities.IsNearlyEqual(pFit.Parameters.ValuesVector, fFit.Parameters.ValuesVector, Math.Sqrt(TestUtilities.TargetPrecision)));
            // covariance
            Assert.IsTrue(TestUtilities.IsNearlyEqual(pFit.Parameters.CovarianceMatrix, fFit.Parameters.CovarianceMatrix, Math.Sqrt(TestUtilities.TargetPrecision)));

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
            UncertainMeasurementFitResult pFit = set.FitToPolynomial(3);

            // fit it to a cubic polynomial
            Func<double[], double, double> ff = delegate(double[] p, double x) {
                return (p[0] + p[1] * x + p[2] * x * x + p[3] * x * x * x);
            };
            UncertainMeasurementFitResult fFit = set.FitToFunction(ff, new double[] { 0, 0, 0, 0 });

            // dimension
            Assert.IsTrue(pFit.Parameters.Count == fFit.Parameters.Count);
            // chi squared
            Assert.IsTrue(TestUtilities.IsNearlyEqual(pFit.GoodnessOfFit.Statistic, fFit.GoodnessOfFit.Statistic, Math.Sqrt(TestUtilities.TargetPrecision)));
            // don't demand super-high precision agreement of parameters and covariance matrix
            // parameters
            Assert.IsTrue(TestUtilities.IsNearlyEqual(pFit.Parameters.ValuesVector, fFit.Parameters.ValuesVector, Math.Pow(TestUtilities.TargetPrecision,0.3)));
            // covariance
            Assert.IsTrue(TestUtilities.IsNearlyEqual(pFit.Parameters.CovarianceMatrix, fFit.Parameters.CovarianceMatrix, Math.Pow(TestUtilities.TargetPrecision,0.3)));

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

            //MultivariateSample sample = new MultivariateSample(2);
            SymmetricMatrix covariance = new SymmetricMatrix(2);

            // create a bunch of small data sets
            for (int i = 0; i < 100; i++) {
                UncertainMeasurementSample data = CreateDataSet(xs, fv, fu, i);
                UncertainMeasurementFitResult fit = data.FitToLine();

                //sample.Add(fit.Parameters.ValuesVector);
                covariance = fit.Parameters.CovarianceMatrix;
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
