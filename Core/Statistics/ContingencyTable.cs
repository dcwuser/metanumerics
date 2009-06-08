using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a contingency table.
    /// </summary>
    /// <remarks><para>Imagine a controlled experiment in which each data point consists of two categories. For example, a
    /// measurement of a driver's car and his nationality, or a person's sex and employment status. Such experiments are
    /// typically undertaken to determine what correlation exists between the categories.</para></remarks>
    /// <seealso cref="BinaryContingencyTable"/>
    public class ContingencyTable {

        private int[,] data;

        /// <summary>
        /// Instantiates a new contingency table.
        /// </summary>
        /// <param name="rows"></param>
        /// <param name="columns"></param>
        public ContingencyTable (int rows, int columns) {
            if (rows < 1) throw new ArgumentOutOfRangeException("rows");
            if (columns < 1) throw new ArgumentOutOfRangeException("columns");
            data = new int[rows,columns];
        }

        /// <summary>
        /// Instantiates a new contingency table with the given data set.
        /// </summary>
        /// <param name="data">A (zero-based) matrix of contingency table entries.</param>
        public ContingencyTable (int[,] data) {
            this.data = new int[data.GetLength(0), data.GetLength(1)];
            for (int r = 0; r < data.GetLength(0); r++) {
                for (int c = 0; c < data.GetLength(1); c++) {
                    if (data[r, c] < 0) throw new ArgumentOutOfRangeException("data");
                    this.data[r, c] = data[r, c];
                }
            }
        }

        /// <summary>
        /// Gets or sets the count in the specified entry.
        /// </summary>
        /// <param name="r">The entry row number.</param>
        /// <param name="c">The entry column number.</param>
        /// <returns>The count in the specified entry.</returns>
        public int this[int r, int c] {
            get {
                if ((r < 0) || (r > data.GetLength(0))) throw new ArgumentOutOfRangeException("r");
                if ((c < 0) || (c > data.GetLength(1))) throw new ArgumentOutOfRangeException("c");
                return (data[r,c]);
            }
            set {
                if ((r < 0) || (r > data.GetLength(0))) throw new ArgumentOutOfRangeException("r");
                if ((c < 0) || (c > data.GetLength(1))) throw new ArgumentOutOfRangeException("c");
                if (value < 0) throw new InvalidOperationException();
                data[r,c] = value;
            }
        }

        /// <summary>
        /// Gets the number of rows in the table.
        /// </summary>
        public int RowCount {
            get {
                return (data.GetLength(0));
            }
        }

        /// <summary>
        /// Gets the number of columns in the table.
        /// </summary>
        public int ColumnCount {
            get {
                return (data.GetLength(1));
            }
        }


        /// <summary>
        /// Gets the total counts in a row.
        /// </summary>
        /// <param name="r">The row number.</param>
        /// <returns>The sum of counts in all entries in the row.</returns>
        public int RowTotal (int r) {
            if ((r < 0) || (r > data.GetLength(0))) throw new ArgumentOutOfRangeException("r");
            int R = 0;
            for (int c = 0; c < data.GetLength(1); c++) {
                R += data[r,c];
            }
            return (R);
        }

        /// <summary>
        /// Gets the total counts in a column.
        /// </summary>
        /// <param name="c">The column number.</param>
        /// <returns>The sum of counts in all entries in the column.</returns>
        public int ColumnTotal (int c) {
            if ((c < 0) || (c > data.GetLength(1))) throw new ArgumentOutOfRangeException("c");
            int C = 0;
            for (int r = 0; r < data.GetLength(0); r++) {
                C += data[r,c];
            }
            return (C);
        }

        internal int[] RowTotals {
            get {
                int[] R = new int[data.GetLength(0)];
                for (int r = 0; r < R.Length; r++) {
                    R[r] = RowTotal(r);
                }
                return (R);
            }
        }

        internal int[] ColumnTotals {
            get {
                int[] C = new int[data.GetLength(1)];
                for (int c = 0; c < C.Length; c++) {
                    C[c] = ColumnTotal(c);
                }
                return (C);
            }
        }

        /// <summary>
        /// Gets the total counts in the table.
        /// </summary>
        public int Total {
            get {
                int N = 0;
                for (int r = 0; r < data.GetLength(0); r++) {
                    for (int c = 0; c < data.GetLength(1); c++) {
                        N += data[r, c];
                    }
                }
                return (N);
            }
        }

        // suppose there are no correlations present; then the expected number of events in a given table entry is just the fraction
        // of the total events proportional to the corresponding row and column totals: N_RC = N*(N_R/N)*(N_C/N) = N_R*N_C/N
        // if this number is large for all entries, then the entries should be distributed normally with mean N_RC and variance Sqrt(N_RC)
        // we can test this hypothesis by measuring the departure of the actual entries from the expected entries using a chi^2 test

        /// <summary>
        /// Performs a Pearson &#x3C7;<sup>2</sup> test for correlation in the table.
        /// </summary>
        /// <returns>The result of the test. The test statistic is &#x3C7;<sup>2</sup> and its likelyhood under the null hypothesis is
        /// the (right) probability to obtain a value as large or larger.</returns>
        /// <remarks><para>The Pearson Pearson &#x3C7;<sup>2</sup> test tests for correlation between the row and column values. If
        /// row and column values are uncorrelated, then the expected number of counts in a table entry is simply proportional to the
        /// totals for its row and column. If that number is large for all entries, then the central limit theorem suggests that the
        /// actual number of counts will be distributed normally with mean equal to the expected value and standard deviation equal
        /// to its square root. The &#x3C7;<sup>2</sup> statistic measures the departure of the actual table from this expectation
        /// in the uncorrelated case, and under this null hypothesis its distribution is known. Having calculated &#x3C7;<sup>2</sup>,
        /// then, we can compute just how unlikely it was to obtain a value as large or larger than the one obtained.</para>
        /// <para>In cases where either the actual or expected counts for some entries are small or zero, the assumptions of the
        /// Pearson &#x3C7;<sup>2</sup> test are violated and it should not be used. For 2 X 2 experiments, the
        /// <see cref="BinaryContingencyTable.FisherExactTest" /> is a viable alternative in these cases.</para></remarks>
        /// <seealso cref="ChiSquaredDistribution"/>
        public TestResult PearsonChiSquaredTest () {

            // get the totals needed to compute expected entries
            int[] R = RowTotals;
            int[] C = ColumnTotals;
            int N = Total;

            // compute chi squared
            double chi2 = 0.0;
            for (int r = 0; r < R.Length; r++) {
                for (int c = 0; c < C.Length; c++) {
                    double n = ((double)R[r]) * ((double)C[c]) / N;
                    double z = data[r, c] - n;
                    chi2 += z * z / n;
                }
            }

            // the degrees of freedom is just the number of entries minus the number of constraints: R*C - R - (C - 1)
            int nu = (R.Length - 1) * (C.Length - 1);

            // return the test result
            return (new TestResult(chi2, new ChiSquaredDistribution(nu)));

        }

    }

    /// <summary>
    /// Represents a 2 X 2 contingency table.
    /// </summary>
    public class BinaryContingencyTable : ContingencyTable {

        /// <summary>
        /// Initializes a new binary contingency table.
        /// </summary>
        public BinaryContingencyTable () : base(2,2) {
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
        /// the same and the opposite correlation as the given one.</para></remarks>
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
