using System;
using System.Collections.Generic;

using TestClassAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestClassAttribute;
using TestMethodAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestMethodAttribute;
using ExpectedExceptionAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.ExpectedExceptionAttribute;
using Assert = Microsoft.VisualStudio.TestTools.UnitTesting.Assert;

using Meta.Numerics;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;


namespace Test {
    
    [TestClass]
    public class DiscreteDistributionTest {

        private DiscreteDistribution[] distributions = GetDistributions();

        public static DiscreteDistribution[] GetDistributions () {
            return (new DiscreteDistribution[] {
                new BernoulliDistribution(0.1),
                new BinomialDistribution(0.2, 30), new BinomialDistribution(0.4, 5),
                new PoissonDistribution(0.54), new PoissonDistribution(5.4), new PoissonDistribution(540.0),
                new DiscreteUniformDistribution(5, 11),
                new GeometricDistribution(0.6),
                new NegativeBinomialDistribution(7.8, 0.4),
                new HypergeometricDistribution(9, 3, 5),
                new DiscreteTestDistribution()
            });
        }

        [TestMethod]
        public void DiscreteDistributionUnitarity () {
            foreach (DiscreteDistribution distribution in distributions) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    distribution.ExpectationValue(delegate (int k) { return (1.0); }), 1.0
                ));
            }
        }

        [TestMethod]
        public void DiscreteDistributionMean () {
            foreach (DiscreteDistribution distribution in distributions) {
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    distribution.ExpectationValue(delegate(int k) { return (k); }), distribution.Mean
                ));
            }
        }

        [TestMethod]
        public void DiscreteDistributionVariance () {
            foreach (DiscreteDistribution distribution in distributions) {
                double m = distribution.Mean;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(
                    distribution.ExpectationValue(delegate(int x) { return (Math.Pow(x-m, 2)); }), distribution.Variance
                ));
            }
        }

        [TestMethod]
        public void DiscreteDistributionSupport () {

            foreach (DiscreteDistribution distribution in distributions) {

                // To left of support
                if (distribution.Support.LeftEndpoint > Int32.MinValue) {
                    int k0 = distribution.Support.LeftEndpoint - 1;
                    Assert.IsTrue(distribution.LeftInclusiveProbability(k0) == 0.0);
                    Assert.IsTrue(distribution.LeftExclusiveProbability(k0) == 0.0);
                    Assert.IsTrue(distribution.RightExclusiveProbability(k0) == 1.0);
                    Assert.IsTrue(distribution.ProbabilityMass(k0) == 0.0);
                }

                // At left end of support
                int k1 = distribution.Support.LeftEndpoint;
                Assert.IsTrue(distribution.LeftInclusiveProbability(k1) == distribution.ProbabilityMass(k1));
                Assert.IsTrue(distribution.LeftExclusiveProbability(k1) == 0.0);
                Assert.IsTrue(distribution.InverseLeftProbability(0.0) == k1);

                // At right end of support
                int k2 = distribution.Support.RightEndpoint;
                Assert.IsTrue(distribution.LeftInclusiveProbability(k2) == 1.0);
                //Assert.IsTrue(distribution.RightExclusiveProbability(k2) == 0.0);

                // To right of support
                if (distribution.Support.RightEndpoint < Int32.MaxValue) {
                    int k3 = distribution.Support.RightEndpoint + 1;
                    Assert.IsTrue(distribution.LeftInclusiveProbability(k3) == 1.0);
                    Assert.IsTrue(distribution.LeftExclusiveProbability(k3) == 1.0);
                    Assert.IsTrue(distribution.RightExclusiveProbability(k3) == 0.0);
                }

            }
        }

        [TestMethod]
        public void DiscreteDistributionProbabilityAxioms () {

            foreach (DiscreteDistribution distribution in distributions) {

                // some of these values will be outside the support, but that's fine, our results should still be consistent with probability axioms
                foreach (int k in TestUtilities.GenerateUniformIntegerValues(-10, +100, 8)) {

                    double DP = distribution.ProbabilityMass(k);
                    Assert.IsTrue(0.0 <= DP && DP <= 1.0);

                    double P = distribution.LeftInclusiveProbability(k);
                    Assert.IsTrue(0.0 <= P && P <= 1.0);

                    double Q = distribution.RightExclusiveProbability(k);
                    Assert.IsTrue(0.0 <= Q && Q <= 1.0);

                    Assert.IsTrue(TestUtilities.IsNearlyEqual(P + Q, 1.0));
                }

            }

        }

        [TestMethod]
        public void DiscreteDistributionInverseCDF () {

            Random rng = new Random(1);
            for (int i = 0; i < 10; i++) {

                double P = rng.NextDouble();

                foreach (DiscreteDistribution distribution in distributions) {
                    int k = distribution.InverseLeftProbability(P);
                    Assert.IsTrue(distribution.LeftExclusiveProbability(k) < P);
                    Assert.IsTrue(P <= distribution.LeftInclusiveProbability(k));
                }


            }

        }

        [TestMethod]
        public void BinomialNegativeBinomialRelation () {

            int k = 2;

            int r = 3;
            double p = 0.4;
            NegativeBinomialDistribution nb = new NegativeBinomialDistribution(r, p);

            int n = r + k;
            BinomialDistribution b = new BinomialDistribution(p, n);

            double nbP = nb.LeftInclusiveProbability(k);
            double bP = b.LeftInclusiveProbability(k);


        }

        [TestMethod]
        public void HypergeometricDistributionSymmetry () {
            foreach (int nPopulation in TestUtilities.GenerateIntegerValues(1, 1000, 4)) {
                foreach (int nSuccess in TestUtilities.GenerateIntegerValues(1, nPopulation, 2)) {
                    foreach (int nDraws in TestUtilities.GenerateIntegerValues(1, nPopulation, 2)) {
                        HypergeometricDistribution d = new HypergeometricDistribution(nPopulation, nSuccess, nDraws);

                        HypergeometricDistribution d1 = new HypergeometricDistribution(nPopulation, nPopulation - nSuccess, nDraws);
                        HypergeometricDistribution d2 = new HypergeometricDistribution(nPopulation, nSuccess, nPopulation - nDraws);
                        HypergeometricDistribution d3 = new HypergeometricDistribution(nPopulation, nDraws, nSuccess);

                        int kMin = d.Support.LeftEndpoint;
                        int kMax = d.Support.RightEndpoint;
                        foreach (int k in TestUtilities.GenerateUniformIntegerValues(kMin, kMin, 4)) {
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(d.ProbabilityMass(k), d1.ProbabilityMass(nDraws - k)));
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(d.ProbabilityMass(k), d2.ProbabilityMass(nSuccess - k)));
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(d.ProbabilityMass(k), d3.ProbabilityMass(k)));
                        }
                    }
                }
            }
        }

        /*
        [TestMethod]
        public void DiscreteContinuousAgreement () {

            DiscreteDistribution dd = new BinomialDistribution(0.6, 7);
            Distribution cd = new DiscreteAsContinuousDistribution(dd);

            Assert.IsTrue(cd.Mean == dd.Mean);
            Assert.IsTrue(cd.StandardDeviation == dd.StandardDeviation);
            Assert.IsTrue(cd.Variance == dd.Variance);
            Assert.IsTrue(cd.Skewness == dd.Skewness);
            Assert.IsTrue(cd.Moment(5) == dd.Moment(5));
            Assert.IsTrue(cd.MomentAboutMean(5) == dd.MomentAboutMean(5));

            // this should cause an interval conversion
            Assert.IsTrue(TestUtilities.IsNearlyEqual(cd.Support.LeftEndpoint, dd.Minimum));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(cd.Support.RightEndpoint, dd.Maximum));

            //Assert.IsTrue(cd.LeftProbability(4.5) == dd.LeftProbability(4));
            //Assert.IsTrue(cd.RightProbability(4.5) == dd.RightProbability(4));

            // Switch LeftProbablity for discrete distributions to be exclusive.
            // This is already the case for internal distributions used for exact null distributions, but not for public discrete distributions.

        }
        */

        [TestMethod]
        public void OutsideDiscreteDistributionSupport () {
            foreach (DiscreteDistribution distribution in distributions) {
                int min = distribution.Support.LeftEndpoint;
                int max = distribution.Support.RightEndpoint;
                if (min > Int32.MinValue) {
                    Assert.IsTrue(distribution.ProbabilityMass(min - 1) == 0.0);
                    Assert.IsTrue(distribution.LeftInclusiveProbability(min - 1) == 0.0);
                    Assert.IsTrue(distribution.RightExclusiveProbability(min - 1) == 1.0);
                    Assert.IsTrue(distribution.LeftExclusiveProbability(min) == 0.0);
                }
                if (distribution.Support.RightEndpoint < Int32.MaxValue) {
                    Assert.IsTrue(distribution.ProbabilityMass(max + 1) == 0.0);
                    Assert.IsTrue(distribution.LeftInclusiveProbability(max + 1) == 1.0);
                    Assert.IsTrue(distribution.RightExclusiveProbability(max) == 0.0);
                }
            }
        }

        [TestMethod]
        public void PoissonBug () {

            PoissonDistribution pd = new PoissonDistribution(0.5);
            double x = pd.InverseLeftProbability(0.7716);
            Console.WriteLine(x);

        }

        [TestMethod]
        public void DiscreteDistributionChisquaredTest () {

            Random rng = new Random(271828);
            foreach (DiscreteDistribution distribution in distributions) {

                // High mean Poisson has at most ~2 counts in all bins
                // Make chi square test re-bin before trying it here
                if (distribution is PoissonDistribution) continue;

                List<int> sample = new List<int>();
                for (int i = 0; i < 128; i++) sample.Add(distribution.GetRandomValue(rng));

                TestResult chi2 = sample.ChiSquaredTest(distribution);
                Console.WriteLine($"{distribution.GetType().Name} {((ChiSquaredDistribution) chi2.Statistic.Distribution).DegreesOfFreedom}");
                Assert.IsTrue(chi2.Probability > 0.01);
            }

        }

        [TestMethod]
        public void DiscreteDistributionRandomValues () {

            // Test of histogram to the distribution that created it should find compatibility.

            foreach (DiscreteDistribution distribution in distributions) {

                int max = distribution.Support.RightEndpoint;
                if (max < 128) {
                    max = max + 1;
                } else {
                    max = (int) Math.Round(distribution.Mean + 2.0 * distribution.StandardDeviation);
                }

                Histogram h = new Histogram(max);

                Random rng = new Random(314159265);
                for (int i = 0; i < 1024; i++) h.Add(distribution.GetRandomValue(rng));
                TestResult result = h.ChiSquaredTest(distribution);
                Assert.IsTrue(result.Probability > 0.01);

            }

        }

        [TestMethod]
        public void BernoulliDistributionFit ()
        {
            double p = 0.25;

            int n0 = 0;
            int n1 = 0;
            Random rng = new Random(3);
            for (int i = 0; i < 100; i++) {
                if (rng.NextDouble() < p) { n1++; } else { n0++; }
            }
            BernoulliFitResult fit = BernoulliDistribution.FitToSample(n0, n1);
            Assert.IsTrue(fit.P.ConfidenceInterval(0.95).ClosedContains(p));
            Assert.IsTrue(fit.GoodnessOfFit.Probability > 0.05);

            Assert.IsTrue(fit.Parameters.Count == 1);
            Assert.IsTrue(fit.Parameters[0].Estimate == fit.P);
        }

        [TestMethod]
        public void DiscreteDistributionBase () {

            DiscreteDistribution D = new DiscreteTestDistribution();

            double M0 = D.ExpectationValue(k => 1.0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(M0, 1.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(M0, D.RawMoment(0)));

            double M1 = D.ExpectationValue(k => k);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(M1, D.Mean));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(M1, D.RawMoment(1)));

            double C2 = D.ExpectationValue(k => MoreMath.Sqr(k - M1));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(C2, D.Variance));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(C2, D.CentralMoment(2)));

            Assert.IsTrue(D.InverseLeftProbability(D.LeftInclusiveProbability(2)) == 2);

        }

    }

    // a minimal implementation to test base methods on abstract DiscreteDistribution class

    public class DiscreteTestDistribution : DiscreteDistribution {

        public override DiscreteInterval Support {
            get { return DiscreteInterval.FromEndpoints(1, 3); }
        }

        public override double ProbabilityMass (int k) {
            switch (k) {
                case 1:
                    return (1.0 / 6.0);
                case 2:
                    return (2.0 / 6.0);
                case 3:
                    return (3.0 / 6.0);
                default:
                    return (0.0);
            }
        }

    }

}
