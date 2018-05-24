using System;
using System.Collections.Generic;

using TestClassAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestClassAttribute;
using TestMethodAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestMethodAttribute;
using ExpectedExceptionAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.ExpectedExceptionAttribute;
using Assert = Microsoft.VisualStudio.TestTools.UnitTesting.Assert;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace Test {
    
    [TestClass]
    public class DistributionTest {

        // We see a reliabile failure of KS or Kuiper tests for Beta distribution with small parameters, e.g. Beta(0.01,0.01).
        // This appears to occur whether we use inverse CDF or x/(x+y) to generate beta deviates.
        // Perhaps it indicates a problem with P computation for beta in this region?

        private static List<ContinuousDistribution> CreateDistributions () {

            List<ContinuousDistribution> distributions = new List<ContinuousDistribution>( new ContinuousDistribution[] {
                new NoncentralChiSquaredDistribution(2, 3.0),
                new CauchyDistribution(1.0, 2.0),
                new UniformDistribution(Interval.FromEndpoints(-2.0,1.0)), new UniformDistribution(Interval.FromEndpoints(7.0, 9.0)),
                new NormalDistribution(3.0,2.0),
                new ExponentialDistribution(2.0),
                new ChiSquaredDistribution(3),
                new StudentDistribution(5),
                new LognormalDistribution(0.2,0.4),
                new WeibullDistribution(2.0, 3.0),
                new LogisticDistribution(-4.0,5.0),
                new FisherDistribution(4.0, 7.0),
                new KuiperDistribution(),
                new KolmogorovDistribution(),
                new TriangularDistribution(1.0,2.0,4.0),
                new BetaDistribution(0.5, 2.0),
                new ParetoDistribution(1.0, 3.0),
                new WaldDistribution(3.0, 1.0),
                new PearsonRDistribution(7),
                new GammaDistribution(0.8), new GammaDistribution(3.0, 5.0), new GammaDistribution(96.2),
                new GumbelDistribution(1.2, 2.3),
                new LaplaceDistribution(4.5, 6.0),
                new ChiDistribution(1), new ChiDistribution(4),
                new RayleighDistribution(3.0),
                new FrechetDistribution(2.9, 4.0)
            });

            // Add some distributions that come from tests.
            Sample small = TestUtilities.CreateSample(distributions[0], 7);
            Sample large = TestUtilities.CreateSample(distributions[1], 127);
            distributions.Add(small.KolmogorovSmirnovTest(distributions[2]).Statistic.Distribution);
            distributions.Add(large.KolmogorovSmirnovTest(distributions[3]).Statistic.Distribution);
            //distributions.Add(small.KuiperTest(distributions[4]).Distribution);
            //distributions.Add(large.KuiperTest(distributions[5]).Distribution);
            //distributions.Add(Sample.MannWhitneyTest(small, large).Distribution);

            BivariateSample two = new BivariateSample();
            two.Add(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 }, new double[] { 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 });
            //distributions.Add(two.SpearmanRhoTest().Distribution);
            //distributions.Add(two.KendallTauTest().Distribution);

            return (distributions);
        }

        private static IEnumerable<ContinuousDistribution> distributions = CreateDistributions();

        public static ContinuousDistribution[] GetDistributions () { return (new ContinuousDistribution[] {
            new CauchyDistribution(1.0, 2.0),
            new UniformDistribution(Interval.FromEndpoints(-2.0,1.0)), new UniformDistribution(Interval.FromEndpoints(7.0, 9.0)),
            new NormalDistribution(3.0,2.0),
            new ExponentialDistribution(2.0),
            new ChiSquaredDistribution(3),
            new StudentDistribution(5),
            new LognormalDistribution(0.2,0.4),
            new WeibullDistribution(2.0, 3.0),
            new LogisticDistribution(-4.0,5.0),
            new FisherDistribution(4.0, 7.0),
            new KuiperDistribution(),
            new KolmogorovDistribution(),
            new TriangularDistribution(1.0,2.0,4.0),
            new BetaDistribution(0.5, 2.0),
            new ParetoDistribution(1.0, 3.0),
            new WaldDistribution(3.0, 1.0),
            new PearsonRDistribution(7),
            new GammaDistribution(0.8), new GammaDistribution(3.0, 5.0), new GammaDistribution(96.2),
            new GumbelDistribution(1.2, 2.3)            });
        }

        private double[] probabilities = new double[] {
            0.00001, 0.01, 0.05, 0.1, 1.0 / 3.0, 1.0 / 2.0, 2.0 / 3.0, 0.9, 0.95, 0.99999
        };

        [TestMethod]
        public void DistributionMedian () {
            foreach (ContinuousDistribution distribution in distributions) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(distribution.Median, distribution.InverseLeftProbability(0.5)));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(distribution.Median, distribution.InverseRightProbability(0.5)));
            }
        }


        [TestMethod]
        public void DistributionSkewness () {
            foreach (ContinuousDistribution distribution in distributions) {
                Console.WriteLine(distribution.GetType().FullName);
                if (!Double.IsNaN(distribution.Skewness)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        distribution.Skewness, distribution.CentralMoment(3) / Math.Pow(distribution.CentralMoment(2), 3.0 / 2.0)
                    ));
                }
            }
        }

        [TestMethod]
        public void DistributionCentralInequality () {
            foreach (ContinuousDistribution distribution in distributions) {
                double mean = distribution.Mean;
                if (Double.IsNaN(mean)) continue;
                double median = distribution.Median;
                double sigma = distribution.StandardDeviation;
                Assert.IsTrue(Math.Abs(mean - median) <= sigma);
            }
        }

        [TestMethod]
        public void DistributionMonotonicity () {
            foreach (ContinuousDistribution distribution in distributions) {
                for (int i = 0; i < (probabilities.Length - 1); i++) {
                    Assert.IsTrue(distribution.InverseLeftProbability(probabilities[i]) < distribution.InverseLeftProbability(probabilities[i+1]));
                    Assert.IsTrue(distribution.InverseRightProbability(probabilities[i]) > distribution.InverseRightProbability(probabilities[i + 1]));
                }
            }
        }

        [TestMethod]
        public void DistributionProbability () {
            foreach (ContinuousDistribution distribution in distributions) {
                // some of these x's will be outside range,
                // but that should just produce zero probability values
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 5, 3)) {
                    DistributionProbabilityTestHelper(distribution, x);
                    DistributionProbabilityTestHelper(distribution, -x);
                }
            }
        }

        private void DistributionProbabilityTestHelper (ContinuousDistribution distribution, double x) {
            double P = distribution.LeftProbability(x);
            double Q = distribution.RightProbability(x);
            Assert.IsTrue((0.0 <= P) && (P <= 1.0));
            Assert.IsTrue((0.0 <= Q) && (Q <= 1.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(P + Q, 1.0));
            double p = distribution.ProbabilityDensity(x);
            Assert.IsTrue(p >= 0.0);
            double h = distribution.Hazard(x);
            if (p > 0.0 && Q > 0.0) Assert.IsTrue(TestUtilities.IsNearlyEqual(h, p / Q));
        }

        [TestMethod]
        public void DistributionUnitarityIntegral () {
            foreach (ContinuousDistribution distribution in distributions) {
                // Gamma distribution has a power law divergence at zero, which our integrator is documented not to handle (look for lower precision)
                if ((distribution is GammaDistribution) && (((GammaDistribution) distribution).Shape < 1.0)) continue;
                Console.WriteLine(distribution.GetType().Name);
                double M0 = FunctionMath.Integrate(distribution.ProbabilityDensity, distribution.Support);
                Console.WriteLine("  1 {0}", M0);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(M0, 1.0));
            }
        }

        [TestMethod]
        public void DistributionMeanIntegral () {
            foreach (ContinuousDistribution distribution in distributions) {
                double mean = distribution.Mean;
                if (Double.IsNaN(mean) || Double.IsInfinity(mean)) continue;
                Func<double, double> f = delegate(double x) {
                    return (distribution.ProbabilityDensity(x) * x);
                };
                double M1 = FunctionMath.Integrate(f, distribution.Support);
                Console.WriteLine("{0} {1} {2}", distribution.GetType().Name, mean, M1);
                if (mean == 0.0) {
                    Assert.IsTrue(Math.Abs(M1) <= TestUtilities.TargetPrecision);
                } else {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(M1, distribution.Mean));
                }
            }
        }

        // test variance
        [TestMethod]
        public void DistributionVarianceIntegral () {
            foreach (ContinuousDistribution distribution in distributions) {
                double v = distribution.Variance;
                // Skip distributions with infinite variance
                if (Double.IsNaN(v) || Double.IsInfinity(v)) continue;
                // Determine target precision, which must be reduced for some cases where an integral singularity means that cannot achieve full precision
                double e = TestUtilities.TargetPrecision;
                GammaDistribution gammaDistribution = distribution as GammaDistribution;
                if ((gammaDistribution != null) && (gammaDistribution.Shape < 1.0)) e = Math.Sqrt(e);
                Func<double, double> f = delegate(double x) {
                    double z = x - distribution.Mean;
                    return (distribution.ProbabilityDensity(x) * z * z);
                };
                double C2 = FunctionMath.Integrate(f, distribution.Support, new IntegrationSettings() { RelativePrecision = e, AbsolutePrecision = 0.0 }).Value;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(C2, distribution.Variance, e));
            }
        }

        // test higher raw moments

        [TestMethod]
        public void DistributionRawMomentIntegral () {
            foreach (ContinuousDistribution distribution in distributions) {
                foreach (int r in TestUtilities.GenerateIntegerValues(2, 32, 8)) {
                    double M = distribution.RawMoment(r);
                    if (Double.IsInfinity(M) || Double.IsNaN(M)) continue; // don't try to do a non-convergent integral
                    Func<double, double> f = delegate(double x) {
                        return (distribution.ProbabilityDensity(x) * Math.Pow(x, r));
                    };
                    try {
                        IntegrationResult MI = FunctionMath.Integrate(f, distribution.Support);
                        Console.WriteLine("{0} {1} {2} {3}", distribution.GetType().Name, r, M, MI);
                        Assert.IsTrue(MI.Estimate.ConfidenceInterval(0.99).ClosedContains(M));
                        /*
                        if (M == 0.0) {
                            Assert.IsTrue(Math.Abs(MI) < TestUtilities.TargetPrecision);
                        } else {
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(M, MI));
                        }
                        */
                    } catch (NonconvergenceException) {
                        Console.WriteLine("{0} {1} {2} NC", distribution.GetType().Name, r, M);
                    }
                }
            }
        }

        // test higher central moments
        [TestMethod]
        public void DistributionCentralMomentIntegral () {
            foreach (ContinuousDistribution distribution in distributions) {
                foreach (int n in TestUtilities.GenerateIntegerValues(2, 24, 8)) {

                    // get the predicted central moment
                    double C = distribution.CentralMoment(n);

                    // don't try to integrate infinite moments
                    if (Double.IsInfinity(C) || Double.IsNaN(C)) continue;

                    if (C == 0.0) continue;

                    IntegrationSettings settings = new IntegrationSettings();
                    if (C == 0.0) {
                        // if moment is zero, use absolute precision
                        settings.AbsolutePrecision = TestUtilities.TargetPrecision;
                        settings.RelativePrecision = 0.0;
                    } else {
                        // if moment in non-zero, use relative precision
                        settings.AbsolutePrecision = 0.0;
                        settings.RelativePrecision = TestUtilities.TargetPrecision;
                    }

                    // do the integral
                    double m = distribution.Mean;
                    Func<double, double> f = delegate(double x) {
                        return (distribution.ProbabilityDensity(x) * MoreMath.Pow(x - m, n));
                    };
                    try {
                        double CI = FunctionMath.Integrate(f, distribution.Support, settings).Value;
                        Console.WriteLine("{0} {1} {2} {3}", distribution.GetType().Name, n, C, CI);
                        if (C == 0.0) {
                            Assert.IsTrue(Math.Abs(CI) < TestUtilities.TargetPrecision);
                        } else {
                            double e = TestUtilities.TargetPrecision;
                            // reduce required precision, because some distributions (e.g. Kolmogorov, Weibull)
                            // have no analytic expressions for central moments, which must therefore be
                            // determined via raw moments and are thus subject to cancelation error
                            // can we revisit this later?
                            if (distribution is WeibullDistribution) e = Math.Sqrt(Math.Sqrt(e));
                            if (distribution is KolmogorovDistribution) e = Math.Sqrt(e);
                            if (distribution is KuiperDistribution) e = Math.Sqrt(Math.Sqrt(e));
                            if (distribution is TriangularDistribution) e = Math.Sqrt(e);
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(C, CI, e));
                        }
                    } catch (NonconvergenceException) {
                        Console.WriteLine("{0} {1} {2} {3}", distribution.GetType().Name, n, C, "NC");
                        // deal with these later; they are integration problems, not distribution problems
                    }
                }
            }
        }

        // test P values
        [TestMethod]
        public void DistributionProbabilityIntegral () {
            Random rng = new Random(1);

            // if integral is very small, we still want to get it very accurately
            IntegrationSettings settings = new IntegrationSettings() { AbsolutePrecision = 0.0 };

            foreach (ContinuousDistribution distribution in distributions) {

                for (int i = 0; i < 4; i++) {
                    double x;
                    if (Double.IsNegativeInfinity(distribution.Support.LeftEndpoint) && Double.IsPositiveInfinity(distribution.Support.RightEndpoint)) {
                        // pick an exponentially distributed random point with a random sign
                        double y = rng.NextDouble();
                        x = - Math.Log(y);
                        if (rng.NextDouble() < 0.5) x = -x;
                    } else if (Double.IsPositiveInfinity(distribution.Support.RightEndpoint)) {
                        // pick an exponentially distributed random point
                        double y = rng.NextDouble();
                        x = distribution.Support.LeftEndpoint - Math.Log(y);
                    } else {
                        // pick a random point within the support
                        x = distribution.Support.LeftEndpoint + rng.NextDouble() * distribution.Support.Width;
                    }
                    double P = FunctionMath.Integrate(distribution.ProbabilityDensity, Interval.FromEndpoints(distribution.Support.LeftEndpoint, x), settings).Value;
                    double Q = FunctionMath.Integrate(distribution.ProbabilityDensity, Interval.FromEndpoints(x, distribution.Support.RightEndpoint), settings).Value;
                    if (!TestUtilities.IsNearlyEqual(P + Q, 1.0)) {
                        // the numerical integral for the triangular distribution can be inaccurate, because
                        // its locally low-polynomial behavior fools the integration routine into thinking it need
                        // not integrate as much near the inflection point as it must; this is a problem
                        // of the integration routine (or arguably the integral), not the triangular distribution,
                        // so skip it here
                        continue;
                    }

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(P, distribution.LeftProbability(x)));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(Q, distribution.RightProbability(x)));
                }

            }
        }

        [TestMethod]
        public void DistributionInvalidProbabilityInput () {

            foreach (ContinuousDistribution distribution in distributions) {

                try {
                    distribution.InverseLeftProbability(1.1);
                    Assert.IsTrue(false);
                } catch (ArgumentOutOfRangeException) {
                    Assert.IsTrue(true);
                }

                try {
                    distribution.InverseLeftProbability(-0.1);
                    Assert.IsTrue(false);
                } catch (ArgumentOutOfRangeException) {
                    Assert.IsTrue(true);
                }

            }

        }

        [TestMethod]
        public void StudentFromNormal () {

            // make sure Student t is consistent with its definition

            // we are going to take a sample that we expect to be t-distributed
            Sample tSample = new Sample();

            // begin with an underlying normal distribution
            ContinuousDistribution xDistribution = new NormalDistribution();

            // compute a bunch of t statistics from the distribution
            Random rng = new Random(314159);
            for (int i = 0; i < 10000; i++) {
                double p = xDistribution.GetRandomValue(rng);
                double q = 0.0;
                for (int j = 0; j < 5; j++) {
                    double x = xDistribution.GetRandomValue(rng);
                    q += x * x;
                }
                q = q / 5;

                double t = p / Math.Sqrt(q);
                tSample.Add(t);
            }

            ContinuousDistribution tDistribution = new StudentDistribution(5);
            TestResult tResult = tSample.KolmogorovSmirnovTest(tDistribution);
            Assert.IsTrue(tResult.Probability > 0.05);

            // Distribution should be demonstrably non-normal
            ContinuousDistribution zDistribution = new NormalDistribution(tDistribution.Mean, tDistribution.StandardDeviation);
            TestResult zResult = tSample.KolmogorovSmirnovTest(zDistribution);
            Assert.IsTrue(zResult.Probability < 0.05);
        }

        [TestMethod]
        public void FisherFromChiSquared () {

            // we will need a RNG
            Random rng = new Random(314159);

            int n1 = 1;
            int n2 = 2;

            // define chi squared distributions
            ContinuousDistribution d1 = new ChiSquaredDistribution(n1);
            ContinuousDistribution d2 = new ChiSquaredDistribution(n2);

            // create a sample of chi-squared variates
            Sample s = new Sample();
            for (int i = 0; i < 250; i++) {
                double x1 = d1.GetRandomValue(rng);
                double x2 = d2.GetRandomValue(rng);
                double x = (x1/n1) / (x2/n2);
                s.Add(x);
            }

            // it should match a Fisher distribution with the appropriate parameters
            ContinuousDistribution f0 = new FisherDistribution(n1, n2);
            TestResult t0 = s.KuiperTest(f0);
            Assert.IsTrue(t0.Probability > 0.05);

            // it should be distinguished from a Fisher distribution with different parameters
            ContinuousDistribution f1 = new FisherDistribution(n1 + 1, n2);
            TestResult t1 = s.KuiperTest(f1);
            Assert.IsTrue(t1.Probability < 0.05);

        }

        [TestMethod]
        public void UniformOrderStatistics () {

            // Check that the order statistics of the uniform distribution are distributed as expected.

            Random rng = new Random(1);
            UniformDistribution u = new UniformDistribution();

            Sample maxima = new Sample();
            Sample minima = new Sample();

            for (int i = 0; i < 100; i++) {

                double maximum = 0.0;
                double minimum = 1.0;
                for (int j = 0; j < 4; j++) {
                    double value = u.GetRandomValue(rng);
                    if (value > maximum) maximum = value;
                    if (value < minimum) minimum = value;
                }

                maxima.Add(maximum);
                minima.Add(minimum);

            }

            // maxima should be distributed according to Beta(n,1)
            TestResult maxTest = maxima.KolmogorovSmirnovTest(new BetaDistribution(4, 1));
            Assert.IsTrue(maxTest.Probability > 0.05);

            // minima should be distributed according to Beta(1,n)
            TestResult minTest = minima.KolmogorovSmirnovTest(new BetaDistribution(1, 4));
            Assert.IsTrue(minTest.Probability > 0.05);


        }

        [TestMethod]
        public void FisherInversion () {

            // x ~ Fisher(a,b) => 1/x ~ Fisher(b,a)

            FisherDistribution f = new FisherDistribution(2.3, 5.6);
            FisherDistribution fi = new FisherDistribution(f.DenominatorDegreesOfFreedom, f.NumeratorDegreesOfFreedom);

            Random rng = new Random(1);
            for (int i = 0; i < 10; i++) {

                double x = f.GetRandomValue(rng);
                double xi = 1.0 / x;

                // LeftProbability <-> RightProbability because as x increases, 1/x decreases
                Assert.IsTrue(TestUtilities.IsNearlyEqual(f.LeftProbability(x), fi.RightProbability(xi)));

            }

        }

        [TestMethod]
        public void GammaFromExponential () {

            // test that x_1 + x_2 + ... + x_n ~ Gamma(n) when z ~ Exponential()

            Random rng = new Random(1);
            ExponentialDistribution eDistribution = new ExponentialDistribution();

            // pick some low values of n so distribution is not simply normal
            foreach (int n in new int[] { 2, 3, 4, 5 }) {
                Sample gSample = new Sample();
                for (int i = 0; i < 100; i++) {

                    double sum = 0.0;
                    for (int j = 0; j < n; j++) {
                        sum += eDistribution.GetRandomValue(rng);
                    }
                    gSample.Add(sum);

                }

                GammaDistribution gDistribution = new GammaDistribution(n);
                TestResult result = gSample.KolmogorovSmirnovTest(gDistribution);
                Assert.IsTrue(result.Probability > 0.05);

            }

        }

        [TestMethod]
        public void InverseGaussianSummation () {

            // X_i ~ IG(\mu,\lambda) \rightarrow \sum_{i=0}^{n} X_i ~ IG(n \mu, n^2 \lambda)
            
            Random rng = new Random(1);
            WaldDistribution d0 = new WaldDistribution(1.0, 2.0);
            List<double> s = new List<double>();
            for (int i = 0; i < 64; i++) {
                s.Add(d0.GetRandomValue(rng) + d0.GetRandomValue(rng) + d0.GetRandomValue(rng));
            }
            WaldDistribution d1 = new WaldDistribution(3.0 * 1.0, 9.0 * 2.0);
            TestResult r = s.KolmogorovSmirnovTest(d1);
            Assert.IsTrue(r.Probability > 0.05);
        }

        [TestMethod]
        public void CauchyStudentAgreement () {

            StudentDistribution S = new StudentDistribution(1);
            CauchyDistribution C = new CauchyDistribution();

            // don't compare moments directly, because NaN != NaN

            foreach (double P in probabilities) {
                double xS = S.InverseLeftProbability(P);
                double xC = C.InverseLeftProbability(P);
                Console.WriteLine("{0} {1} {2}", P, xS, xC);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(xS, xC));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(S.ProbabilityDensity(xS), C.ProbabilityDensity(xC)));
            }
            
        }

        [TestMethod]
        public void CauchyFWHM () {
            // Check that FWHM really is the full-width at half-maximum.
            CauchyDistribution D = new CauchyDistribution(1.0, 2.0);
            double p = D.ProbabilityDensity(D.Median);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(D.ProbabilityDensity(D.Median - D.FullWithAtHalfMaximum / 2.0), p / 2.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(D.ProbabilityDensity(D.Median + D.FullWithAtHalfMaximum / 2.0), p / 2.0));
        }

        [TestMethod]
        public void DistributionRandomDeviates () {
            // Check that random deviates generated match distribution
            foreach (ContinuousDistribution distribution in distributions) {
                Console.WriteLine(distribution.GetType().Name);
                Sample s = TestUtilities.CreateSample(distribution, 128);
                TestResult r = s.KolmogorovSmirnovTest(distribution);
                Assert.IsTrue(r.Probability > 0.01);

            }
        }

        [TestMethod]
        public void OutsideDistributionSupport () {
            foreach (ContinuousDistribution distribution in distributions) {
                Interval support = distribution.Support;
                if (support.LeftEndpoint > Double.NegativeInfinity) {
                    Assert.IsTrue(distribution.ProbabilityDensity(support.LeftEndpoint - 1.0) == 0.0);
                    Assert.IsTrue(distribution.LeftProbability(support.LeftEndpoint - 1.0) == 0.0);
                    Assert.IsTrue(distribution.RightProbability(support.LeftEndpoint - 1.0) == 1.0);
                }
                if (support.RightEndpoint < Double.PositiveInfinity) {
                    Assert.IsTrue(distribution.ProbabilityDensity(support.RightEndpoint + 1.0) == 0.0);
                    Assert.IsTrue(distribution.LeftProbability(support.RightEndpoint + 1.0) == 1.0);
                    Assert.IsTrue(distribution.RightProbability(support.RightEndpoint + 1.0) == 0.0);
                }
            }
        }

        [TestMethod]
        public void DistributionBase () {

            // test that implementations on base Distribution classes function and agree with overridden implementations

            ContinuousDistribution d = new TestDistribution();
            ContinuousDistribution t = new TriangularDistribution(0.0, 1.0, 1.0);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(d.Mean, t.Mean));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(d.StandardDeviation, t.StandardDeviation));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(d.Skewness, t.Skewness));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(d.Median, t.Median));

        }

        private double[] cs = new double[] { 0.02, 0.2, 2.0, 20.0, 200.0 };

        [TestMethod]
        public void BetaInversion () {

            // test that beta distribution is accurately inverted over a wide range of a, b, P

            foreach (double a in TestUtilities.GenerateRealValues(0.01, 100.0, 8)) {
                foreach (double b in cs) {

                    BetaDistribution B = new BetaDistribution(a, b);

                    foreach (double P in probabilities) {
                        Console.WriteLine("a={0} b={1} P={2}", a, b, P);
                        double x = B.InverseLeftProbability(P);
                        Console.WriteLine("  x={0} P(x)={1}", x, B.LeftProbability(x));
                        // we would like to test that P(x) = P, but floating point limitations prevent us from meeting this standard
                        // P(x) changes so fast at extremes that sometimes even the minimal change in x causes a change
                        // in P(x) larger than our target precision; so instead we test that our routine gets us
                        // as close as it can, by checking that P(x-e) < P < P(x+e)
                        double Px = B.LeftProbability(x);
                        double Pxm = B.LeftProbability(Math.Min(0.0, x * (1.0 - TestUtilities.TargetPrecision)));
                        double Pxp = B.LeftProbability(Math.Max(x * (1.0 + TestUtilities.TargetPrecision), 1.0));
                        Assert.IsTrue((Pxm <= P) && (P <= Pxp));
                    }
                }
            }

        }

        [TestMethod]
        public void TransformedBetaMoments () {

            // For x ~ B(a,b), kth cumulant of \ln\left(\frac{x}{1-x}\right) = \psi^{(k-1)}(a) \pm \psi^{(k-1)}(b)

            BetaDistribution D = new BetaDistribution(2.0, 3.0);

            double m = D.ExpectationValue(x => Math.Log(x / (1.0 - x)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(m, AdvancedMath.Psi(D.Alpha) - AdvancedMath.Psi(D.Beta)));

            double c2 = D.ExpectationValue(x => MoreMath.Pow(Math.Log(x / (1.0 - x)) - m, 2));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(c2, AdvancedMath.Psi(1, D.Alpha) + AdvancedMath.Psi(1, D.Beta)));

            double c3 = D.ExpectationValue(x => MoreMath.Pow(Math.Log(x / (1.0 - x)) - m, 3));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(c3, AdvancedMath.Psi(2, D.Alpha) - AdvancedMath.Psi(2, D.Beta)));

            double c4 = D.ExpectationValue(x => MoreMath.Pow(Math.Log(x / (1.0 - x)) - m, 4));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(c4 - 3.0 * c2 * c2, AdvancedMath.Psi(3, D.Alpha) + AdvancedMath.Psi(3, D.Beta)));

        }

        [TestMethod]
        public void DistributionMomentTranslations () {

            int n = 4;

            foreach (ContinuousDistribution distribution in distributions) {

                if (Double.IsNaN(distribution.Mean)) continue;
                Console.Write(distribution.GetType().Name);

                // Convert central moments to raw moments
                double[] centralInputs = new double[n];
                for (int k = 0; k < n; k++) centralInputs[k] = distribution.CentralMoment(k);
                double[] rawOutputs = MomentMath.CentralToRaw(distribution.Mean, centralInputs);
                Assert.IsTrue(rawOutputs.Length == n);
                for (int k = 0; k < n; k++) Assert.IsTrue(TestUtilities.IsNearlyEqual(rawOutputs[k], distribution.RawMoment(k)));

                // Convert cumulants to central moments
                double[] cumulantInputs = new double[n];
                for (int k = 0; k < n; k++) cumulantInputs[k] = distribution.Cumulant(k);
                double[] centralOutputs = MomentMath.CumulantToCentral(cumulantInputs);
                Assert.IsTrue(centralOutputs.Length == n);
                for (int k = 0; k < n; k++) Assert.IsTrue(TestUtilities.IsNearlyEqual(centralOutputs[k], distribution.CentralMoment(k)));
            }

        }

#if FUTURE
        [TestMethod]
        public void GammaFitTest () {

            double k0 = 0.05;
            double t0 = 2.0;
            GammaDistribution G = new GammaDistribution(k0, t0);
            Random rng = new Random(1);

            Sample S = new Sample();
            for (int i = 0; i < 1000; i++) {
                S.Add(G.GetRandomValue(rng));
            }

            Console.WriteLine("mean sample={0} distribution={1}", S.Mean, G.Mean);
            Console.WriteLine("variance sample={0} distribution={1}", S.Variance, G.Variance);

            Console.WriteLine("moment estimate k={0} t={1}", S.Mean * S.Mean / S.Variance, S.Variance / S.Mean);

            double q = 0;
            foreach (double x in S) {
                q += Math.Log(x / S.Mean);
            }
            q /= S.Count;
            q = -q;
            Console.WriteLine("q = {0}, Log(k)-Psi(k)={1}", q, Math.Log(k0) - AdvancedMath.Psi(k0));

            double ke0 = S.Mean * S.Mean / S.Variance;
            Console.WriteLine("ke0={0} Log(ke0)-Psi(ke0)={1}", ke0, Math.Log(ke0) - AdvancedMath.Psi(ke0));
            double ke1 = 1.0 / (2.0 * q - 2.0 / 3.0 * q * q + 4.0 / 9.0 * q * q * q - 14.0 / 135.0 * q * q * q * q);
            Console.WriteLine("ke1={0} Log(ke1)-Psi(ke1)={1}", ke1, Math.Log(ke1) - AdvancedMath.Psi(ke1));
            double ke2 = 1.0 / (q - AdvancedMath.EulerGamma);
            Console.WriteLine("ke2={0} Log(ke2)-Psi(ke2)={1}", ke2, Math.Log(ke2) - AdvancedMath.Psi(ke2));

        }

        private double PsiDeficit (double x) {
            if (x > 4.0) {
                return ((1.0 / 2.0) / x + (1.0 / 12.0) / (x * x) - (1.0 / 120.0) * (x * x * x * x) + (1.0 / 252.0) * (x * x * x * x * x * x));
            } else {
                return (Math.Log(x) - AdvancedMath.Psi(x));
            }
        }

#endif

    }

    // This is a very simple distribution we define in order to test the implementations of method on the
    // Distribution base class. Most of these methods are not invoked by our "real" distribution classes,
    // because we override them with better implementations.

    public class TestDistribution : ContinuousDistribution {

        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, 1.0));
            }
        }

        public override double ProbabilityDensity (double x) {
            if (x < 0.0) {
                return (0.0);
            } else if (x > 1.0) {
                return (1.0);
            } else {
                return (2.0 * x);
            }
        }

    }
}
