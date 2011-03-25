using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace Test {

    [TestClass]
    public class BivariateSampleTest {

        [TestMethod]
        public void BivariateSampleManipulations () {

            BivariateSample s = new BivariateSample();
            s.Add(1.0, 3.0);
            s.Add(2.0, 2.0);
            s.Add(3.0, 1.0);
            Assert.IsTrue(s.Count == 3);

            Assert.IsTrue(s.Remove(2.0, 2.0));
            Assert.IsTrue(s.Count == 2);
            Assert.IsFalse(s.Remove(2.0, 2.0));

            Assert.IsTrue(s.Contains(1.0, 3.0));
            Assert.IsFalse(s.Contains(3.0, 3.0));

            s.Clear();
            Assert.IsTrue(s.Count == 0);

        }


        [TestMethod]
        public void BivariateSampleCopy () {

            // test independency of copy

            BivariateSample sample1 = new BivariateSample();
            sample1.Add(1.0, 2.0);

            BivariateSample sample2 = sample1.Copy();
            sample2.Add(3.0, 4.0);

            Assert.IsTrue(sample1.Count == 1);
            Assert.IsTrue(sample2.Count == 2);

        }

        [TestMethod]
        public void PearsonRDistribution () {

            Random rng = new Random(1);

            // pick some underlying distributions for the sample variables, which must be normal but can have any parameters
            NormalDistribution xDistribution = new NormalDistribution(1, 2);
            NormalDistribution yDistribution = new NormalDistribution(3, 4);

            // try this for several sample sizes, all low so that we see the difference from the normal distribution
            // n = 3 maxima at ends; n = 4 uniform; n = 5 semi-circular "mound"; n = 6 parabolic "mound"
            foreach (int n in new int[] { 3, 4, 5, 8 }) {
                Console.WriteLine("n={0}", n);

                // find r values
                Sample rSample = new Sample();
                for (int i = 0; i < 50; i++) {

                    // to get each r value, construct a bivariate sample of the given size with no cross-correlation
                    BivariateSample xySample = new BivariateSample();
                    for (int j = 0; j < n; j++) {
                        xySample.Add(xDistribution.GetRandomValue(rng), yDistribution.GetRandomValue(rng));
                    }
                    double r = xySample.PearsonRTest().Statistic;
                    rSample.Add(r);

                }

                // check whether r is distributed as expected
                TestResult result = rSample.KolmogorovSmirnovTest(new PearsonRDistribution(n));
                Console.WriteLine("P={0}", result.LeftProbability);
                Assert.IsTrue(result.LeftProbability < 0.95);
            }


        }

        [TestMethod]
        public void LinearLogisticRegression () {

            // do a set of logistic regression fits
            // make sure not only that the fit parameters are what they should be, but that their variances/covariances are as returned

            Random rng = new Random(314159);

            // define logistic parameters
            double a0 = 1.0; double b0 = -1.0 / 2.0;
            //double a0 = -0.5; double b0 = 2.0;

            // keep track of sample of returned a and b fit parameters
            BivariateSample ps = new BivariateSample();

            // also keep track of returned covariance estimates
            // since these vary slightly from fit to fit, we will average them
            double caa = 0.0;
            double cbb = 0.0;
            double cab = 0.0;

            // do 50 fits
            for (int k = 0; k < 50; k++) {

                // generate a synthetic data set
                BivariateSample s = new BivariateSample();
                for (int i = 0; i < 50; i++) {
                    double x = 2.0 * rng.NextDouble() - 1.0;
                    double ez = Math.Exp(a0 + b0 * x);
                    double P = ez / (1.0 + ez);
                    if (rng.NextDouble() < P) {
                        s.Add(x, 1.0);
                    } else {
                        s.Add(x, 0.0);
                    }
                }

                // do the regression
                FitResult r = s.LinearLogisticRegression();

                // record best fit parameters
                double a = r.Parameter(0).Value;
                double b = r.Parameter(1).Value;
                ps.Add(a, b);

                // record estimated covariances
                caa += r.Covariance(0, 0);
                cbb += r.Covariance(1, 1);
                cab += r.Covariance(0, 1); 

            }

            caa /= ps.Count;
            cbb /= ps.Count;
            cab /= ps.Count;

            // check that mean parameter estimates are what they should be: the underlying population parameters
            Assert.IsTrue(ps.X.PopulationMean.ConfidenceInterval(0.95).ClosedContains(a0));
            Assert.IsTrue(ps.Y.PopulationMean.ConfidenceInterval(0.95).ClosedContains(b0));

            // check that parameter covarainces are what they should be: the reported covariance estimates
            Assert.IsTrue(ps.X.PopulationVariance.ConfidenceInterval(0.95).ClosedContains(caa));
            Assert.IsTrue(ps.Y.PopulationVariance.ConfidenceInterval(0.95).ClosedContains(cbb));
            Assert.IsTrue(ps.PopulationCovariance.ConfidenceInterval(0.95).ClosedContains(cab));

        }

        [TestMethod]
        public void BivariateLinearRegression () {

            // do a set of logistic regression fits
            // make sure not only that the fit parameters are what they should be, but that their variances/covariances are as returned

            Random rng = new Random(314159);

            // define logistic parameters
            double a0 = 2.0; double b0 = -1.0;

            // keep track of sample of returned a and b fit parameters
            BivariateSample ps = new BivariateSample();

            // also keep track of returned covariance estimates
            // since these vary slightly from fit to fit, we will average them
            double caa = 0.0;
            double cbb = 0.0;
            double cab = 0.0;

            // also keep track of test statistics
            Sample fs = new Sample();

            // do 50 fits
            for (int k = 0; k < 50; k++) {

                // we should be able to draw x's from any distribution; noise should be drawn from a normal distribution
                //Distribution xd = new ExponentialDistribution();
                Distribution xd = new LogisticDistribution();
                Distribution nd = new NormalDistribution(0.0, 2.0);

                // generate a synthetic data set
                BivariateSample s = new BivariateSample();
                for (int i = 0; i < 50; i++) {
                    double x = xd.GetRandomValue(rng);
                    double y = a0 + b0 * x + nd.GetRandomValue(rng);
                    s.Add(x, y);
                }

                // do the regression
                FitResult r = s.LinearRegression();

                // record best fit parameters
                double a = r.Parameter(0).Value;
                double b = r.Parameter(1).Value;
                ps.Add(a, b);

                // record estimated covariances
                caa += r.Covariance(0, 0);
                cbb += r.Covariance(1, 1);
                cab += r.Covariance(0, 1);

                // record the fit statistic
                fs.Add(r.GoodnessOfFit.Statistic);
                Console.WriteLine("F={0}", r.GoodnessOfFit.Statistic);

            }

            caa /= ps.Count;
            cbb /= ps.Count;
            cab /= ps.Count;

            // check that mean parameter estimates are what they should be: the underlying population parameters
            Assert.IsTrue(ps.X.PopulationMean.ConfidenceInterval(0.95).ClosedContains(a0));
            Assert.IsTrue(ps.Y.PopulationMean.ConfidenceInterval(0.95).ClosedContains(b0));

            Console.WriteLine("{0} {1}", caa, ps.X.PopulationVariance);
            Console.WriteLine("{0} {1}", cbb, ps.Y.PopulationVariance);

            // check that parameter covarainces are what they should be: the reported covariance estimates
            Assert.IsTrue(ps.X.PopulationVariance.ConfidenceInterval(0.95).ClosedContains(caa));
            Assert.IsTrue(ps.Y.PopulationVariance.ConfidenceInterval(0.95).ClosedContains(cbb));
            Assert.IsTrue(ps.PopulationCovariance.ConfidenceInterval(0.95).ClosedContains(cab));

            // check that F is distributed as it should be
            Console.WriteLine(fs.KolmogorovSmirnovTest(new FisherDistribution(2, 48)).LeftProbability);

        }

        [TestMethod]
        public void BivariateLinearRegressionGoodnessOfFitDistribution () {

            // create uncorrelated x and y values
            // the distribution of F-test statistics returned by linear fits should follow the expected F-distribution

            Random rng = new Random(987654321);
            NormalDistribution xd = new NormalDistribution(1.0, 2.0);
            NormalDistribution yd = new NormalDistribution(-3.0, 4.0);

            Sample fs = new Sample();

            for (int i = 0; i < 127; i++) {
                BivariateSample xys = new BivariateSample();
                for (int j = 0; j < 7; j++) {
                    xys.Add(xd.GetRandomValue(rng), yd.GetRandomValue(rng));
                }
                double f = xys.LinearRegression().GoodnessOfFit.Statistic;
                fs.Add(f);
            }

            Distribution fd = new FisherDistribution(1, 5);
            Console.WriteLine("{0} v. {1}", fs.PopulationMean, fd.Mean);
            TestResult t = fs.KolmogorovSmirnovTest(fd);
            Console.WriteLine(t.LeftProbability);
            Assert.IsTrue(t.LeftProbability < 0.95);

        }

    }
}
