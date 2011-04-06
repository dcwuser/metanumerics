using System;

using Meta.Numerics.Functions;
using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a 2 X 2 contingency table.
    /// </summary>
    /// <remarks>
    /// <para>Binary contingency tables are the most common kind of contingency table. Any experiment with a treatment and
    /// a control group, and binary measured outcome, can be represented by a binary contingency table.</para>
    /// </remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/Contingency_table" />
    public class BinaryContingencyTable : ContingencyTable {

        /// <summary>
        /// Initializes a new binary contingency table.
        /// </summary>
        public BinaryContingencyTable () : base(2, 2) {
        }

        /// <summary>
        /// Initializes a new binary contingency table with the given entries.
        /// </summary>
        /// <param name="data">A two-dimensional matrix of table entries.</param>
        public BinaryContingencyTable (int[,] data) : base(data) {
            if ((data.GetLength(0) != 2) || (data.GetLength(1) != 2)) throw new DimensionMismatchException();
        }

        // for row 0, the odds of landing in column 0 vs. column1 are N[0,0]/N[0,1]
        // for row 1, the odds of landing in column 0 vs. column1 are N[1,0]/N[1,1]
        // the ratio of these odds are (N[0,0]/N[0,1])/(N[1,0]/N[1,1])
        // if this ratio is 1, your odds of landing in either column are not changed depending on your row
        // if this ratio is incompatible with 1, the odds of landing in either column depend on your row

        /// <summary>
        /// Computes the odds ratio of table.
        /// </summary>
        /// <remarks><para>For entries in the first row, the odds of landing in the first column are given by N[0,0] / N[0,1].
        /// For entries in the second row, the odds of landing in the first column are given by N[1,0] / N[1,1]. The odds
        /// ratio is the ratio of these two odds. An odds ratio significantly different from 1 indicates a correlation between
        /// row and column values.</para></remarks>
        /// <seealso cref="LogOddsRatio"/>
        public UncertainValue OddsRatio {
            get {
                double r = (1.0 * this[0, 0] / this[0, 1]) / (1.0 * this[1, 0] / this[1, 1]);
                double dlnr = Math.Sqrt(1.0 / this[0, 0] + 1.0 / this[0, 1] + 1.0 / this[1, 0] + 1.0 / this[1, 1]);
                return (new UncertainValue(r, r * dlnr));
            }
        }

        /// <summary>
        /// Computes the log of the odds ratio.
        /// </summary>
        /// <seealso cref="OddsRatio"/>
        public UncertainValue LogOddsRatio {
            get {
                double lnr = Math.Log(this[0, 0]) - Math.Log(this[0, 1]) - Math.Log(this[1, 0]) + Math.Log(this[1, 1]);
                double dlnr = Math.Sqrt(1.0 / this[0, 0] + 1.0 / this[0, 1] + 1.0 / this[1, 0] + 1.0 / this[1, 1]);
                return (new UncertainValue(lnr, dlnr));
            }
        }

        // the chi2 test fails when the expected values of some entries are small, because the entry values are not normally distributed
        // in the 2X2 case, the probably of any given entry value can be computed; it is given by
        //        N_R1! N_R2! N_C1! N_C2!
        // P = ----------------------------
        //      N_11! N_12! N_21! N_22! N!
        // the exact test computes this probability for the actual table and for all other tables having the same row and column totals
        // the test statistic is the fraction of tables that have a lower probability than the actual table; it is uniformly distributed
        // on [0,1]. if this fraction is very small, the actual table is particularly unlikely under the null hypothesis of no correlation

        /// <summary>
        /// Performs a Fisher exact test.
        /// </summary>
        /// <returns>The results of the test. The test statistic is the summed probability of all tables exhibiting equal or stronger correlations,
        /// and its likelyhood under the null hypothesis is the (left) probability to obtain a smaller value. Note that, in this case, the test
        /// statistic itself is the likelyhood.</returns>
        /// <remarks><para>The Fisher exact test tests for correlations between row and column entries. It is a robust, non-parametric test,
        /// which, unlike the &#x3C7;<sup>2</sup> test (see <see cref="ContingencyTable.PearsonChiSquaredTest"/>), can safely be used for tables
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
            // store row and column totals
            int[] R = this.RowTotals;
            int[] C = this.ColumnTotals;
            int N = this.Total;

            // compute the critical probability
            double x =
                AdvancedIntegerMath.LogFactorial(R[0]) +
                AdvancedIntegerMath.LogFactorial(R[1]) +
                AdvancedIntegerMath.LogFactorial(C[0]) +
                AdvancedIntegerMath.LogFactorial(C[1]) -
                AdvancedIntegerMath.LogFactorial(N);
            double lnPc = x -
                AdvancedIntegerMath.LogFactorial(this[0, 0]) -
                AdvancedIntegerMath.LogFactorial(this[0, 1]) -
                AdvancedIntegerMath.LogFactorial(this[1, 0]) -
                AdvancedIntegerMath.LogFactorial(this[1, 1]);

            // compute all possible 2 X 2 matrices with these row and column totals
            // compute the total probability of getting a matrix as or less probable than the measured one
            double P = 0.0;
            int[,] test = new int[2, 2];
            int min = Math.Max(C[0] - R[1], 0);
            int max = Math.Min(R[0], C[0]);
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

            // return the result
            return (new TestResult(P, new UniformDistribution(Interval.FromEndpoints(0.0, 1.0))));
        }

    }


}