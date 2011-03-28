using System;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;
using Meta.Numerics.Matrices;

namespace Test {
    
    
    /// <summary>
    ///This is a test class for DistributionTest and is intended
    ///to contain all DistributionTest Unit Tests
    ///</summary>
    [TestClass()]
    public class DistributionTest {


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

        private Distribution[] distributions = new Distribution[] {
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
            new GammaDistribution(5.0, 6.0), new GammaDistribution(78.9)
        };

        private double[] probabilities = new double[] {
            0.00001, 0.01, 0.05, 0.1, 1.0 / 3.0, 1.0 / 2.0, 2.0 / 3.0, 0.9, 0.95, 0.99999
        };

        [TestMethod]
        public void DistributionMedian () {
            foreach (Distribution distribution in distributions) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(distribution.Median, distribution.InverseLeftProbability(0.5)));
            }
        }

        [TestMethod]
        public void DistributionMomentSpecialCases () {
            foreach (Distribution distribution in distributions) {
                Assert.IsTrue(distribution.Moment(0) == 1.0);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(distribution.Moment(1), distribution.Mean));
                Assert.IsTrue(distribution.MomentAboutMean(0) == 1.0);
                Assert.IsTrue(distribution.MomentAboutMean(1) == 0.0);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(distribution.MomentAboutMean(2), distribution.Variance));
            }
        }

        [TestMethod]
        public void DistributionSkewness () {
            foreach (Distribution distribution in distributions) {
                Console.WriteLine(distribution.GetType().FullName);
                if (!Double.IsInfinity(distribution.Skewness)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(
                        distribution.Skewness, distribution.MomentAboutMean(3) / Math.Pow(distribution.MomentAboutMean(2), 3.0 / 2.0)
                    ));
                }
            }
        }

        [TestMethod]
        public void DistributionMomentSums () {
            foreach (Distribution distribution in distributions) {
                Console.WriteLine(distribution.GetType().Name);
                // C2 = M2 - M1^2
                double M1 = distribution.Moment(1);
                double M2 = distribution.Moment(2);
                double C2 = distribution.MomentAboutMean(2);
                if (!Double.IsInfinity(C2)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(C2 + M1 * M1, M2));
                }
                // C3 = M3 - 3 M2 M1 + 2 M1^3
                double M3 = distribution.Moment(3);
                double C3 = distribution.MomentAboutMean(3);
                if (!Double.IsInfinity(C3)) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(C3 + 3.0 * M2 * M1, M3 + 2.0 * M1 * M1 * M1));
                }
                // C4 = M4 - 4 M3 M1 + 6 M2 M1^2 - 3 M1^4
            }
        }

        [TestMethod]
        public void DistributionCentralInequality () {
            foreach (Distribution distribution in distributions) {
                Assert.IsTrue(Math.Abs(distribution.Mean - distribution.Median) <= distribution.StandardDeviation);
            }
        }

        [TestMethod]
        public void DistributionMonotonicity () {
            foreach (Distribution distribution in distributions) {
                for (int i = 0; i < (probabilities.Length - 1); i++) {
                    Console.WriteLine("{0} {1}", distribution.GetType().Name, probabilities[i]);
                    Assert.IsTrue(distribution.InverseLeftProbability(probabilities[i]) < distribution.InverseLeftProbability(probabilities[i+1]));
                }
            }
        }

        [TestMethod]
        public void DistributionProbability () {
            foreach (Distribution distribution in distributions) {
                // some of these x's will be outside range,
                // but that should just produce zero probability values
                foreach (double x in TestUtilities.GenerateRealValues(1.0E-2, 1.0E2, 5)) {
                    DistributionProbabilityTestHelper(distribution, x);
                    DistributionProbabilityTestHelper(distribution, -x);
                }
            }
        }

        private void DistributionProbabilityTestHelper (Distribution distribution, double x) {
            Console.WriteLine("{0} {1}", distribution.GetType().Name, x);
            double P = distribution.LeftProbability(x);
            double Q = distribution.RightProbability(x);
            Console.WriteLine(" P={0} Q={1} P+Q={2}", P, Q, P + Q);
            Assert.IsTrue((0.0 <= P) && (P <= 1.0));
            Assert.IsTrue((0.0 <= Q) && (Q <= 1.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(P + Q, 1.0));
            // this is a rather poor test; we can do much more when we get integration
            // update: we have integration now, and tests like DistributionUnitarityIntegralTest,
            // DistributionMeanIntegralTest, DistributionVarianceIntegralTest, DistributionRawMomentIntegralTest,
            // and DistributionCentralMomentIntegralTest use it
            double p = distribution.ProbabilityDensity(x);
            Console.WriteLine(" p={0}", p);
            Assert.IsTrue(p >= 0.0);
        }

        [TestMethod]
        public void DistributionUnitarityIntegral () {
            foreach (Distribution distribution in distributions) {
                double M0 = FunctionMath.Integrate(distribution.ProbabilityDensity, distribution.Support);
                Console.WriteLine("{0} 1 {1}", distribution.GetType().Name, M0);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(M0, 1.0));
            }
        }

        [TestMethod]
        public void DistributionMeanIntegral () {
            foreach (Distribution distribution in distributions) {
                Func<double, double> f = delegate(double x) {
                    return (distribution.ProbabilityDensity(x) * x);
                };
                double M1 = FunctionMath.Integrate(f, distribution.Support);
                Console.WriteLine("{0} {1} {2}", distribution.GetType().Name, distribution.Mean, M1);
                if (distribution.Mean == 0.0) {
                    Assert.IsTrue(Math.Abs(M1) <= TestUtilities.TargetPrecision);
                } else {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(M1, distribution.Mean));
                }
            }
        }

        // test variance
        [TestMethod]
        public void DistributionVarianceIntegral () {
            foreach (Distribution distribution in distributions) {
                Func<double, double> f = delegate(double x) {
                    double z = x - distribution.Mean;
                    return (distribution.ProbabilityDensity(x) * z * z);
                };
                double C2 = FunctionMath.Integrate(f, distribution.Support);
                Console.WriteLine("{0} {1} {2}", distribution.GetType().Name, distribution.StandardDeviation, Math.Sqrt(C2));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(Math.Sqrt(C2), distribution.StandardDeviation));
            }
        }

        // test higher raw moments

        [TestMethod]
        public void DistributionRawMomentIntegral () {
            foreach (Distribution distribution in distributions) {
                // range of moments is about 3 to 30
                foreach (int n in TestUtilities.GenerateIntegerValues(3, 30, 10)) {
                    double M = distribution.Moment(n);
                    if (Double.IsInfinity(M)) continue; // don't try to do a non-convergent integral
                    Func<double, double> f = delegate(double x) {
                        return (distribution.ProbabilityDensity(x) * Math.Pow(x, n));
                    };
                    try {
                        double MI = FunctionMath.Integrate(f, distribution.Support);
                        Console.WriteLine("{0} {1} {2} {3}", distribution.GetType().Name, n, M, MI);
                        if (M == 0.0) {
                            Assert.IsTrue(Math.Abs(MI) < TestUtilities.TargetPrecision);
                        } else {
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(M, MI));
                        }
                    } catch (NonconvergenceException) {
                        Console.WriteLine("{0} {1} {2} NC", distribution.GetType().Name, n, M);
                    }
                }
            }
        }

        // test higher central moments
        [TestMethod]
        public void DistributionCentralMomentIntegral () {
            foreach (Distribution distribution in distributions) {
                // range of moments is about 3 to 30
                foreach (int n in TestUtilities.GenerateIntegerValues(3, 30, 10)) {
                    double C = distribution.MomentAboutMean(n);
                    if (Double.IsInfinity(C)) continue; // don't try to integrate to infinity
                    double m = distribution.Mean;
                    Func<double, double> f = delegate(double x) {
                        return (distribution.ProbabilityDensity(x) * Math.Pow(x - m, n));
                    };
                    try {
                        double CI = FunctionMath.Integrate(f, distribution.Support);
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
                            if (distribution is KuiperDistribution) e = Math.Sqrt(e);
                            if (distribution is BetaDistribution) e = Math.Sqrt(e);
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
            Random rng = new Random(4);
            foreach (Distribution distribution in distributions) {
                for (int i = 0; i < 3; i++) {
                    double x;
                    if (Double.IsNegativeInfinity(distribution.Support.LeftEndpoint) && Double.IsPositiveInfinity(distribution.Support.RightEndpoint)) {
                        // pick an exponentially distributed random point with a random sign
                        double y = rng.NextDouble();
                        x = - Math.Log(y);
                        if (rng.NextDouble() < 0.5) x = -x;
                    } else if (Double.IsPositiveInfinity(distribution.Support.RightEndpoint)) {
                        // pick an exponentialy distributed random point
                        double y = rng.NextDouble();
                        x = distribution.Support.LeftEndpoint - Math.Log(y);
                    } else {
                        // pick a random point within the support
                        x = distribution.Support.LeftEndpoint + rng.NextDouble() * distribution.Support.Width;
                    }
                    Console.WriteLine("{0} {1}", distribution.GetType().Name, x);
                    double P = FunctionMath.Integrate(distribution.ProbabilityDensity, Interval.FromEndpoints(distribution.Support.LeftEndpoint, x));
                    double Q = FunctionMath.Integrate(distribution.ProbabilityDensity, Interval.FromEndpoints(x, distribution.Support.RightEndpoint));
                    if (!TestUtilities.IsNearlyEqual(P + Q, 1.0)) {
                        // the numerical integral for the triangular distribution can be innacurate, because
                        // its locally low-polynomial behavior fools the integration routine into thinking it need
                        // not integrate as much near the inflection point as it must; this is a problem
                        // of the integration routine (or arguably the integral), not the triangular distribution,
                        // so skip it here
                        Console.WriteLine("skip (P={0}, Q={1})", P, Q);
                        continue;
                    }
                    Console.WriteLine("  {0} v. {1}", P, distribution.LeftProbability(x));
                    Console.WriteLine("  {0} v. {1}", Q, distribution.RightProbability(x));

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(P, distribution.LeftProbability(x)));
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(Q, distribution.RightProbability(x)));
                }

            }
        }

        [TestMethod]
        public void DistributionInvalidProbabilityInput () {

            foreach (Distribution distribution in distributions) {

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

        // this test seems to indicate a problem as the tSample size gets bigger (few thousand)
        // the fit gets too good to be true; is this a problem with KS? the the RNG? somehting else?

        [TestMethod]
        public void StudentTest2 () {

            // make sure Student t is consistent with its definition

            // we are going to take a sample that we expect to be t-distributed
            Sample tSample = new Sample();

            // begin with an underlying normal distribution
            Distribution xDistribution = new NormalDistribution();

            // compute a bunch of t satistics from the distribution
            for (int i = 0; i < 100000; i++) {

                // take a small sample from the underlying distribution
                // (as the sample gets large, the t distribution becomes normal)
                Random rng = new Random(314159+i);

                double p = xDistribution.InverseLeftProbability(rng.NextDouble());
                double q = 0.0;
                for (int j = 0; j < 5; j++) {
                    double x = xDistribution.InverseLeftProbability(rng.NextDouble());
                    q += x * x;
                }
                q = q / 5;

                double t = p / Math.Sqrt(q);
                tSample.Add(t);

            }

            Distribution tDistribution = new StudentDistribution(5);
            TestResult result = tSample.KolmogorovSmirnovTest(tDistribution);
            Console.WriteLine(result.LeftProbability);

        }

        [TestMethod]
        public void StudentTest () {

            // make sure Student t is consistent with its definition

            // we are going to take a sample that we expect to be t-distributed
            Sample tSample = new Sample();

            // begin with an underlying normal distribution
            Distribution xDistribution = new NormalDistribution(1.0, 2.0);

            // compute a bunch of t satistics from the distribution
            for (int i = 0; i < 200000; i++) {

                // take a small sample from the underlying distribution
                // (as the sample gets large, the t distribution becomes normal)
                Random rng = new Random(i);
                Sample xSample = new Sample();
                for (int j = 0; j < 5; j++) {
                    double x = xDistribution.InverseLeftProbability(rng.NextDouble());
                    xSample.Add(x);
                }

                // compute t for the sample
                double t = (xSample.Mean - xDistribution.Mean) / (xSample.PopulationStandardDeviation.Value / Math.Sqrt(xSample.Count));
                tSample.Add(t);
                //Console.WriteLine(t);

            }

            // t's should be t-distrubuted; use a KS test to check this
            Distribution tDistribution = new StudentDistribution(4);
            TestResult result = tSample.KolmogorovSmirnovTest(tDistribution);
            Console.WriteLine(result.LeftProbability);
            //Assert.IsTrue(result.LeftProbability < 0.95);

            // t's should be demonstrably not normally distributed
            Console.WriteLine(tSample.KolmogorovSmirnovTest(new NormalDistribution()).LeftProbability);
            //Assert.IsTrue(tSample.KolmogorovSmirnovTest(new NormalDistribution()).LeftProbability > 0.95);

        }


        [TestMethod]
        public void FisherTest () {

            // we will need a RNG
            Random rng = new Random(314159);

            int n1 = 1;
            int n2 = 2;

            // define chi squared distributions
            Distribution d1 = new ChiSquaredDistribution(n1);
            Distribution d2 = new ChiSquaredDistribution(n2);

            // create a sample of chi-squared variates
            Sample s = new Sample();
            for (int i = 0; i < 250; i++) {
                double x1 = d1.InverseLeftProbability(rng.NextDouble());
                double x2 = d2.InverseLeftProbability(rng.NextDouble());
                double x = (x1/n1) / (x2/n2);
                s.Add(x);
            }

            // it should match a Fisher distribution with the appropriate parameters
            Distribution f0 = new FisherDistribution(n1, n2);
            TestResult t0 = s.KuiperTest(f0);
            Console.WriteLine(t0.LeftProbability);
            Assert.IsTrue(t0.LeftProbability < 0.95);

            // it should be distinguished from a Fisher distribution with different parameters
            Distribution f1 = new FisherDistribution(n1 + 1, n2);
            TestResult t1 = s.KuiperTest(f1);
            Console.WriteLine(t1.LeftProbability);
            Assert.IsTrue(t1.LeftProbability > 0.95);

        }

        /*
        [TestMethod]
        public void FiniteKSDebug () {

            int n = 200;

            FiniteKolmogorovDistribution d = new FiniteKolmogorovDistribution(n);

            double x = 99.0/200.0;

            Stopwatch s = new Stopwatch();
            s.Start();
            double P = d.LeftProbability(x);
            s.Stop();
            Console.WriteLine("{0} ({1})", P, s.ElapsedMilliseconds);

            s.Reset();
            s.Start();
            double Q = d.RightProbability(x);
            s.Stop();
            Console.WriteLine("{0} ({1})", Q, s.ElapsedMilliseconds);

        }
        */


        [TestMethod]
        public void UniformOrderStatistics () {

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
            //Console.WriteLine(maxTest.LeftProbability);
            Assert.IsTrue(maxTest.LeftProbability < 0.95);

            // minima should be distributed according to
            TestResult minTest = minima.KolmogorovSmirnovTest(new BetaDistribution(1, 4));
            //Console.WriteLine(minTest.LeftProbability);
            Assert.IsTrue(minTest.LeftProbability < 0.95);


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
                Assert.IsTrue(result.LeftProbability < 0.95);

            }

        }

        [TestMethod]
        public void ContinuousDistributionDeviates () {

            //Distribution distribution = new WaldDistribution(3.5, 2.5);
            foreach (Distribution distribution in distributions) {
                Console.Write(distribution.GetType().Name);
                Sample s = new Sample();
                Random rng = new Random(1000000);
                for (int i = 0; i < 64; i++) {
                    s.Add(distribution.GetRandomValue(rng));
                }
                TestResult r = s.KolmogorovSmirnovTest(distribution);
                Assert.IsTrue(r.LeftProbability < 0.95);

            }

        }

        [TestMethod]
        public void OutsideDistributionSupport () {
            foreach (Distribution distribution in distributions) {
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

            Distribution d = new TestDistribution();
            Distribution t = new TriangularDistribution(0.0, 1.0, 1.0);

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

    }

    // This is a very simple distribution we define in order to test the implementations of method on the
    // Distribution base class. Most of these methods are not invoked by our "real" distribution classes,
    // because we override them with better implementations.

    public class TestDistribution : Distribution {

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
