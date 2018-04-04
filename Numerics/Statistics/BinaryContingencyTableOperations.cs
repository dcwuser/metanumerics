using System;
using System.Diagnostics;

using Meta.Numerics.Functions;
using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Exposes properties which are only defined for a 2 X 2 contingency table.
    /// </summary>
    public sealed class BinaryContingencyTableOperations {

        internal BinaryContingencyTableOperations (int[,] counts) {
            Debug.Assert(counts != null);
            Debug.Assert(counts.GetLength(0) == 2);
            Debug.Assert(counts.GetLength(1) == 2);
            this.counts = counts;
        }

        private int[,] counts;

        /// <summary>
        /// Gets an estimate of the odds ratio in the underlying population.
        /// </summary>
        /// <remarks>
        /// <para>For entries in the first row, the odds of landing in the first column are given by N[0,0] / N[0,1].
        /// For entries in the second row, the odds of landing in the first column are given by N[1,0] / N[1,1].
        /// The odds ratio is the ratio of these two odds. An odds ratio significantly different from 1 indicates a
        /// correlation between row and column values.</para>
        /// <para>Note that the odds ratio is inverted under the exchange of rows or columns.</para>
        /// </remarks>
        /// <seealso href="https://en.wikipedia.org/wiki/Odds_ratio"/>
        /// <seealso cref="LogOddsRatio"/>
        public UncertainValue OddsRatio {
            get {
                double r = ((double) counts[0, 0] / counts[0, 1]) / ((double) counts[1, 0] / counts[1, 1]);
                double dlnr = Math.Sqrt(1.0 / counts[0, 0] + 1.0 / counts[0, 1] + 1.0 / counts[1, 0] + 1.0 / counts[1, 1]);
                return (new UncertainValue(r, r * dlnr));
            }
        }

        /// <summary>
        /// Gets an estimate of the log of the odds ratio in the underlying population.
        /// </summary>
        /// <remarks>
        /// <para>A log odds ratio significantly different from 0 indicates a correlation between row
        /// and column values.</para>
        /// <para>The sign of the log odds ratio is reversed under the exchange of rows or columns.</para>
        /// </remarks>
        /// <seealso cref="OddsRatio"/>
        /// <seealso href="https://en.wikipedia.org/wiki/Odds_ratio"/>
        public UncertainValue LogOddsRatio {
            get {
                double lnr = Math.Log(counts[0, 0]) - Math.Log(counts[0, 1]) - Math.Log(counts[1, 0]) + Math.Log(counts[1, 1]);
                double dlnr = Math.Sqrt(1.0 / counts[0, 0] + 1.0 / counts[0, 1] + 1.0 / counts[1, 0] + 1.0 / counts[1, 1]);
                return (new UncertainValue(lnr, dlnr));
            }
        }

        /// <summary>
        /// Gets the phi coefficient.
        /// </summary>
        /// <seealso href="https://en.wikipedia.org/wiki/Phi_coefficient"/>
        public double Phi {
            get {
                double p = counts[1, 1] * counts[0, 0] - counts[1, 0] * counts[0, 1];
                double q = Math.Sqrt((counts[0, 0] + counts[0, 1]) * (counts[1, 0] + counts[1, 1]) *
                    (counts[0, 0] + counts[1, 0]) * (counts[0, 1] + counts[1, 1]));
                return (p / q);
            }
        }
        // Add non-population values
        // Add confusion matrix façade?

        /// <summary>
        /// Performs a McNemar test.
        /// </summary>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>McNemar's test is appropriate to a rather specialized circumstance, so if you aren't sure this is the
        /// right test, it probably isn't.</para>
        /// <para>The circumstance to which McNemar's test applies is paired binary measurements. Each member of
        /// a set is measured twice, once before and once after some treatment, and each measurement has a binary
        /// outcome. McNemar's test can be used to determine whether the treatment had any systematic effect. The null
        /// hypothesis is that it does not, i.e. that any discrepancies between the first and second measurements
        /// were random.</para>
        /// </remarks>
        /// <see href="http://en.wikipedia.org/wiki/McNemar%27s_test"/>
        public TestResult McNemarTest () {

            double t = MoreMath.Sqr(counts[0, 1] - counts[1, 0]) / (counts[0, 1] + counts[1, 0]);

            return (new TestResult("t", t, TestType.RightTailed, new ChiSquaredDistribution(1)));

            // Adapt this to: (1) support one-sidedness (our squaring (b-c) currently prevents this)
            // and (2) be exact for small counts.
        }


        /// <summary>
        /// Performs a Fisher exact test.
        /// </summary>
        /// <returns>The results of the test. The test statistic is the summed probability of all tables exhibiting equal or stronger correlations,
        /// and its likelihood under the null hypothesis is the (left) probability to obtain a smaller value. Note that, in this case, the test
        /// statistic itself is the likelihood.</returns>
        /// <remarks><para>The Fisher exact test tests for correlations between row and column entries. It is a robust, non-parametric test,
        /// which, unlike the &#x3C7;<sup>2</sup> test (see <see cref="ContingencyTable{R,C}.PearsonChiSquaredTest"/>), can safely be used for tables
        /// with small, even zero-valued, entries.</para>
        /// <para>The Fisher test computes, under the null hypothesis of no correlation, the exact probability of all 2 X 2 tables with the
        /// same row and column totals as the given table. It then sums the probabilities of all tables that are as or less probable than
        /// the given table. In this way it determines the total probability of obtaining a 2 X 2 table which is at least as improbable
        /// as the given one.</para>
        /// <para>The test is two-sided, i.e. when considering less probable tables it does not distinguish between tables exhibiting
        /// the same and the opposite correlation as the given one.</para>
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Fisher_exact_test"/>
        public TestResult FisherExactTest () {

            // Store row and column totals
            int[] R = new int[] { counts[0, 0] + counts[0, 1], counts[1, 0] + counts[1, 1] };
            int[] C = new int[] { counts[0, 0] + counts[1, 0], counts[0, 1] + counts[1, 1] };
            int N = counts[0, 0] + counts[0, 1] + counts[1, 0] + counts[1, 1];
            Debug.Assert(R[0] + R[1] == N);
            Debug.Assert(C[0] + C[1] == N);

            // Compute the critical probability of the measured matrix given the marginal totals
            // and the assumption of no association.
            double x =
                AdvancedIntegerMath.LogFactorial(R[0]) +
                AdvancedIntegerMath.LogFactorial(R[1]) +
                AdvancedIntegerMath.LogFactorial(C[0]) +
                AdvancedIntegerMath.LogFactorial(C[1]) -
                AdvancedIntegerMath.LogFactorial(N);
            double lnPc = x -
                AdvancedIntegerMath.LogFactorial(counts[0, 0]) -
                AdvancedIntegerMath.LogFactorial(counts[0, 1]) -
                AdvancedIntegerMath.LogFactorial(counts[1, 0]) -
                AdvancedIntegerMath.LogFactorial(counts[1, 1]);

            // Compute all possible 2 X 2 matrices with these row and column totals.
            // Find the total probability of getting a matrix as or less probable than the measured one.
            double P = 0.0;
            int[,] test = new int[2, 2];
            int min = Math.Max(C[0] - R[1], 0);
            int max = Math.Min(R[0], C[0]);

            // If the relevant hypergeometric distribution over i is symmetric, then lnP for the i on
            // the other side is exactly equal to lnPc. But floating point noise means that it might come 
            // in ever so slightly higher, causing us to incorrectly exclude that corresponding bin. This
            // caused a bug. So handle the symmetric case specially.
            if ((R[0] == R[1]) || (C[0] == C[1])) {
                double mid = (min + max) / 2.0;
                if (counts[0, 0] == mid) {
                    P = 1.0;
                } else {
                    int imax = (counts[0, 0] < mid) ? counts[0, 0] : min + (max - counts[0, 0]);
                    for (int i = min; i <= imax; i++) {
                        test[0, 0] = i;
                        test[0, 1] = R[0] - i;
                        test[1, 0] = C[0] - i;
                        test[1, 1] = R[1] - test[1, 0];
                        double lnP = x -
                            AdvancedIntegerMath.LogFactorial(test[0, 0]) -
                            AdvancedIntegerMath.LogFactorial(test[0, 1]) -
                            AdvancedIntegerMath.LogFactorial(test[1, 0]) -
                            AdvancedIntegerMath.LogFactorial(test[1, 1]);
                        P += Math.Exp(lnP);
                    }
                    P *= 2.0;
                    Debug.Assert(P <= 1.0);
                }
            } else {
                for (int i = min; i <= max; i++) {
                    test[0, 0] = i;
                    test[0, 1] = R[0] - i;
                    test[1, 0] = C[0] - i;
                    test[1, 1] = R[1] - test[1, 0];
                    double lnP = x -
                        AdvancedIntegerMath.LogFactorial(test[0, 0]) -
                        AdvancedIntegerMath.LogFactorial(test[0, 1]) -
                        AdvancedIntegerMath.LogFactorial(test[1, 0]) -
                        AdvancedIntegerMath.LogFactorial(test[1, 1]);
                    if (lnP <= lnPc) P += Math.Exp(lnP);
                }
            }

            return (new TestResult("P", P, TestType.LeftTailed, new UniformDistribution(Interval.FromEndpoints(0.0, 1.0))));
        }

    }

}