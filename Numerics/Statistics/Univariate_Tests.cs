using System;
using System.Collections.Generic;

using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {
    public static partial class Univariate {

        /// <summary>
        /// Tests whether the sample is compatible with the given distribution.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <param name="distribution">The distribution.</param>
        /// <returns>The test result. The test statistic is the D statistic and the probability is the chance of
        /// obtaining such a large value of D under the assumption that the sample is drawn from the given distribution.</returns>
        /// <remarks>
        /// <para>The null hypothesis of the Kolmogorov-Smirnov (KS) test is that the sample is drawn from the given continuous distribution.
        /// The test statsitic D is the maximum deviation of the sample's empirical distribution function (EDF) from
        /// the distribution's cumulative distribution function (CDF). A high value of the test statistic, corresponding
        /// to a low right tail probability, indicates that the sample distribution disagrees with the given distribution
        /// to a degree unlikely to arise from statistical fluctuations.</para>
        /// <para>For small sample sizes, we compute the null distribution of D exactly. For large sample sizes, we use an accurate
        /// asympototic approximation. Therefore it is safe to use this method for all sample sizes.</para>
        /// <para>A variant of this test, <see cref="Sample.KolmogorovSmirnovTest(Sample, Sample)"/>, allows you to non-parametrically
        /// test whether two samples are drawn from the same underlying distribution, without having to specify that distribution.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="distribution"/> is <see langword="null"/>.</exception>
        /// <exception cref="InsufficientDataException">There is no data in the sample.</exception>
        /// <seealso cref="KolmogorovDistribution"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"/>
        public static TestResult KolmogorovSmirnovTest (this IReadOnlyList<double> sample, ContinuousDistribution distribution) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (distribution == null) throw new ArgumentNullException(nameof(distribution));

            int n = sample.Count;
            if (n < 1) throw new InsufficientDataException();

            // Compute the D statistic, which is the maximum of the D+ and D- statistics.
            double DP, DM;
            ComputeDStatistics(sample, distribution, out DP, out DM);
            double D = Math.Max(DP, DM);

            ContinuousDistribution DDistribution;
            if (n < 32) {
                DDistribution = new TransformedDistribution(new KolmogorovExactDistribution(n), 0.0, 1.0 / n);
            } else {
                DDistribution = new TransformedDistribution(new KolmogorovAsymptoticDistribution(n), 0.0, 1.0 / Math.Sqrt(n));
            }
            return (new TestResult("D", D, TestType.RightTailed, DDistribution));

        }


        /// <summary>
        /// Tests whether the sample is compatible with the given distribution.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <param name="distribution">The distribution.</param>
        /// <returns>The test result. The test statistic is the V statistic and the chance to obtain such a large
        /// value of V under the assumption that the sample is drawn from the given distribution.</returns>
        /// <remarks>
        /// <para>Like the Kolmogorov-Smirnov test (<see cref="Univariate.KolmogorovSmirnovTest(IReadOnlyList{Double},ContinuousDistribution)"/>),
        /// Kuiper's test compares the EDF of the sample to the CDF of the given distribution.</para>
        /// <para>For small sample sizes, we compute the null distribution of V exactly. For large sample sizes, we use an accurate
        /// asympototic approximation. Therefore it is safe to use this method for all sample sizes.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="distribution"/> is <see langword="null"/>.</exception>
        /// <exception cref="InsufficientDataException">There is no data in the sample.</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Kuiper%27s_test"/>
        public static TestResult KuiperTest (this IReadOnlyList<double> sample, ContinuousDistribution distribution) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (distribution == null) throw new ArgumentNullException(nameof(distribution));

            int n = sample.Count;
            if (n < 1) throw new InsufficientDataException();

            // Compute the V statistic, which is the sum of the D+ and D- statistics.
            double DP, DM;
            ComputeDStatistics(sample, distribution, out DP, out DM);
            double V = DP + DM;

            ContinuousDistribution VDistribution;
            if (n < 32) {
                VDistribution = new TransformedDistribution(new KuiperExactDistribution(n), 0.0, 1.0 / n);
            } else {
                VDistribution = new TransformedDistribution(new KuiperAsymptoticDistribution(n), 0.0, 1.0 / Math.Sqrt(n));
            }
            return (new TestResult("V", V, TestType.RightTailed, VDistribution));

        }

        private static void ComputeDStatistics (IReadOnlyList<double> sample, ContinuousDistribution distribution, out double D1, out double D2) {

            List<int> order = GetSortOrder(sample);

            D1 = 0.0;
            D2 = 0.0;

            double P1 = 0.0;
            for (int i = 0; i < sample.Count; i++) {

                // Compute the theoretical CDF value at the data-point.
                double x = sample[order[i]];
                double P = distribution.LeftProbability(x);

                // compute the experimental CDF value at x+dx
                double P2 = (i + 1.0) / sample.Count;

                // Compare the the theoretical and experimental CDF values at x-dx and x+dx;
                // Update D if it is the largest seperation yet observed.
                double DU = P2 - P;
                if (DU > D1) D1 = DU;
                double DL = P - P1;
                if (DL > D2) D2 = DL;

                // Remember the experimental CDF value at x+dx;
                // it will be the experimental CDF value at x-dx for the next step
                P1 = P2;

            }

        }

    }
}
