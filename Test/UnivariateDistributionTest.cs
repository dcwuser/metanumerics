using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Statistics.Distributions;

namespace Test {

    [TestClass]
    public class UnivariateDistributionTest {

        public UnivariateDistribution[] Distributions;

        [TestInitialize]
        public void Initialize () {
            ContinuousDistribution[] continuous = DistributionTest.GetDistributions();
            DiscreteDistribution[] discrete = DiscreteDistributionTest.GetDistributions();
            Distributions = new UnivariateDistribution[continuous.Length + discrete.Length];
            Array.Copy(continuous, 0, Distributions, 0, continuous.Length);
            Array.Copy(discrete, 0, Distributions, continuous.Length, discrete.Length);
        }

        [TestMethod]
        public void UnivariateDistributionRawMomentSpecialCases () {
            foreach (UnivariateDistribution distribution in Distributions) {
                Console.WriteLine(distribution.GetType().Name);
                Assert.IsTrue(distribution.RawMoment(0) == 1.0);
                if (Double.IsNaN(distribution.Mean)) {
                    Assert.IsTrue(Double.IsNaN(distribution.RawMoment(1)));
                } else {
                    // near-equality added for Binomial distribution
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(distribution.RawMoment(1), distribution.Mean));
                }
            }
        }

        [TestMethod]
        public void UnivariateDistributionCentralMomentSpecialCases () {
            foreach (UnivariateDistribution distribution in Distributions) {
                Assert.IsTrue(distribution.CentralMoment(0) == 1.0);
                if (!Double.IsNaN(distribution.Mean)) Assert.IsTrue(distribution.CentralMoment(1) == 0.0);
                if (Double.IsNaN(distribution.Variance)) {
                    Assert.IsTrue(Double.IsNaN(distribution.CentralMoment(2)));
                } else {
                    // near equality rather than perfect equality added for Weibull distribution
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(distribution.CentralMoment(2), distribution.Variance));
                }
            }
        }


        [TestMethod]
        public void UnivariantDistributionCumulantSpecialCases () {
            foreach (UnivariateDistribution distribution in Distributions) {
                Console.WriteLine(distribution.GetType().Name);
                Assert.IsTrue(distribution.Cumulant(0) == 0.0);
                if (Double.IsNaN(distribution.Mean)) {
                    Assert.IsTrue(Double.IsNaN(distribution.Cumulant(1)));
                } else {
                    Assert.IsTrue(distribution.Cumulant(1) == distribution.Mean);
                }
                if (Double.IsNaN(distribution.Variance)) {
                    Assert.IsTrue(Double.IsNaN(distribution.Cumulant(2)));
                } else {
                    // Near-equality added for Weibull
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(distribution.Cumulant(2), distribution.Variance));
                }
            }
        }

        [TestMethod]
        public void UnivariateDistributionSkewness () {
            foreach (UnivariateDistribution distribution in Distributions) {
                Console.WriteLine(distribution.GetType().Name);
                if (Double.IsNaN(distribution.ExcessKurtosis)) continue;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(distribution.Skewness, distribution.CentralMoment(3) / Math.Pow(distribution.CentralMoment(2), 3.0 / 2.0)));
            }
        }

        [TestMethod]
        public void UnivariateDistributionExcessKurtosis () {
            foreach (UnivariateDistribution distribution in Distributions) {
                Console.WriteLine(distribution.GetType().Name);
                if (Double.IsNaN(distribution.ExcessKurtosis)) continue;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(distribution.ExcessKurtosis, distribution.Cumulant(4) / MoreMath.Sqr(distribution.Cumulant(2))));
            }
        }

        [TestMethod]
        public void UnivariateDistributionMomentSums () {
            foreach (UnivariateDistribution distribution in Distributions) {
                Console.WriteLine(distribution.GetType().Name);
                // C2 = M2 - M1^2
                double M1 = distribution.RawMoment(1);
                double M2 = distribution.RawMoment(2);
                double C2 = distribution.CentralMoment(2);
                if (Double.IsNaN(M1) || Double.IsNaN(M2)) continue;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(C2 + M1 * M1, M2));
                // C3 = M3 - 3 M2 M1 + 2 M1^3
                double M3 = distribution.RawMoment(3);
                double C3 = distribution.CentralMoment(3);
                if (Double.IsNaN(M3)) continue;
                Assert.IsTrue(TestUtilities.IsNearlyEqual(C3 + 3.0 * M2 * M1, M3 + 2.0 * M1 * M1 * M1));
                // C4 = M4 - 4 M3 M1 + 6 M2 M1^2 - 3 M1^4
            }
        }

        [TestMethod]
        public void CumulantToCentralAndRaw () {

            int n = 8;

            foreach (UnivariateDistribution d in Distributions) {
                Console.WriteLine(d.GetType().Name);

                // Problems with NaN/Infinity
                if (d is CauchyDistribution) continue;
                if (d is StudentDistribution) continue;
                if (d is FisherDistribution) continue;
                if (d is ParetoDistribution) continue;

                // Real problems
                if (d is KolmogorovDistribution) continue;
                if (d is KuiperDistribution) continue;

                // From cumulants to central and raw moments

                double[] inK = new double[n];
                for (int r = 0; r < n; r++) inK[r] = d.Cumulant(r);

                double[] outC = MomentMath.CumulantToCentral(inK);
                for (int r = 0; r < n; r++) Console.WriteLine("r={0} K={1} -> C={2} v C={3}", r, inK[r], outC[r], d.CentralMoment(r));
                for (int r = 0; r < n; r++) Assert.IsTrue(TestUtilities.IsNearlyEqual(outC[r], d.CentralMoment(r)));

                double[] outM = MomentMath.CumulantToRaw(inK);
                for (int r = 0; r < n; r++) Console.WriteLine("r={0} K={1} -> M={2} v M={3}", r, inK[r], outM[r], d.RawMoment(r));
                for (int r = 0; r < n; r++) Assert.IsTrue(TestUtilities.IsNearlyEqual(outM[r], d.RawMoment(r)));

            }
        }

    }
}
