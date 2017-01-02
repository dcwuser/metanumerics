using System;
using System.Collections.Generic;

using Meta.Numerics.Functions;
using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a contingency table.
    /// </summary>
    /// <remarks><para>Imagine a controlled experiment in which each data point consists of two category values. For example, a
    /// measurement of a driver's car and his nationality, or a person's sex and employment status. Such experiments are
    /// typically undertaken to determine what correlation exists between the categories.</para></remarks>
    /// <seealso cref="BinaryContingencyTable"/>
    /// <seealso href="http://en.wikipedia.org/wiki/Contingency_table" />
    public class ContingencyTable {

        private string name;
        //private NameCollection rowNames;
        //private NameCollection columnNames;

        private int rows, cols;
        private int[,] counts;
        private int[] rCounts;
        private int[] cCounts;
        private int tCounts;

        /// <summary>
        /// Initializes a new contingency table.
        /// </summary>
        /// <param name="rows">The number of rows, which must be positive.</param>
        /// <param name="columns">The number of columns, which must be positive.</param>
        public ContingencyTable (int rows, int columns) {
            if (rows < 1) throw new ArgumentOutOfRangeException("rows");
            if (columns < 1) throw new ArgumentOutOfRangeException("columns");
            this.rows = rows;
            this.cols = columns;
            this.counts = new int[rows,columns];
            this.rCounts = new int[rows];
            this.cCounts = new int[columns];
            //rowNames = new NameCollection(rows);
            //columnNames = new NameCollection(columns);
        }

        /// <summary>
        /// Initializes a new contingency table with the given data set.
        /// </summary>
        /// <param name="data">A (zero-based) matrix of contingency table entries.</param>
        public ContingencyTable (int[,] data) {
            if (data == null) throw new ArgumentNullException("data");
            this.rows = data.GetLength(0);
            this.cols = data.GetLength(1);
            this.counts = new int[rows, cols];
            this.rCounts = new int[rows];
            this.cCounts = new int[cols];
            for (int r = 0; r < rows; r++) {
                for (int c = 0; c < cols; c++) {
                    if (data[r, c] < 0) throw new ArgumentOutOfRangeException("data");
                    this[r, c] = data[r, c];
                }
            }
            //rowNames = new NameCollection(rows);
            //columnNames = new NameCollection(cols);
        }

        /// <summary>
        /// Gets or sets the name of the table.
        /// </summary>
        public string Name {
            get {
                return (name);
            }
            set {
                name = value;
            }
        }

        /// <summary>
        /// Gets or sets the count in the specified cell.
        /// </summary>
        /// <param name="r">The (zero-based) row index.</param>
        /// <param name="c">The (zero-based) column index.</param>
        /// <returns>The count in the specified cell.</returns>
        public int this[int r, int c] {
            get {
                if ((r < 0) || (r >= rows)) throw new ArgumentOutOfRangeException("r");
                if ((c < 0) || (c >= cols)) throw new ArgumentOutOfRangeException("c");
                return (counts[r,c]);
            }
            set {
                if ((r < 0) || (r >= rows)) throw new ArgumentOutOfRangeException("r");
                if ((c < 0) || (c >= cols)) throw new ArgumentOutOfRangeException("c");
                if (value < 0) throw new InvalidOperationException();
                rCounts[r] += (value - counts[r, c]);
                cCounts[c] += (value - counts[r, c]);
                tCounts += (value - counts[r, c]);
                counts[r,c] = value;
            }
        }

        

        /// <summary>
        /// Increments the count in the specifed cel.
        /// </summary>
        /// <param name="r">The (zero-based) row index.</param>
        /// <param name="c">The (zero-based) column index.</param>
        public void Increment (int r, int c) {
            if ((r < 0) || (r >= rows)) throw new ArgumentOutOfRangeException("r");
            if ((c < 0) || (c >= cols)) throw new ArgumentOutOfRangeException("c");
            rCounts[r]++;
            cCounts[c]++;
            tCounts++;
            counts[r, c]++;
        }

        /// <summary>
        /// Gets the number of rows in the table.
        /// </summary>
        public virtual int RowCount {
            get {
                return (rows);
            }
        }

        /// <summary>
        /// Gets the number of columns in the table.
        /// </summary>
        public virtual int ColumnCount {
            get {
                return (cols);
            }
        }


        /// <summary>
        /// Gets the total counts in the specified row.
        /// </summary>
        /// <param name="r">The (zero-based) row index.</param>
        /// <returns>The sum of all counts in the specified row.</returns>
        public virtual int RowTotal (int r) {
            if ((r < 0) || (r >= rows)) throw new ArgumentOutOfRangeException("r");
            return (rCounts[r]);
        }



        /// <summary>
        /// Gets the total counts in the specified column.
        /// </summary>
        /// <param name="c">The (zero-based) column index.</param>
        /// <returns>The sum of all counts in the specified column.</returns>
        public virtual int ColumnTotal (int c) {
            if ((c < 0) || (c >= cols)) throw new ArgumentOutOfRangeException("c");
            return (cCounts[c]);
        }

        internal int[] RowTotals {
            get {
                int[] R = new int[counts.GetLength(0)];
                for (int r = 0; r < R.Length; r++) {
                    R[r] = RowTotal(r);
                }
                return (R);
            }
        }

        internal int[] ColumnTotals {
            get {
                int[] C = new int[counts.GetLength(1)];
                for (int c = 0; c < C.Length; c++) {
                    C[c] = ColumnTotal(c);
                }
                return (C);
            }
        }

        /// <summary>
        /// Gets the total counts in the table.
        /// </summary>
        public virtual int Total {
            get {
                return (tCounts);
            }
        }

        // row and column probabilites are the fraction of total counts that a given row/column total represents
        // conditional probablities are the fraction of row/column totals that a cell represents
        // to get the uncertainties in these probabilities, use
        //   (df)^2 = (df/da)^2 (da)^2 + (df/db)^2 (db)^2 + ...
        // and use (dn)^2 = n for each cell count


        /// <summary>
        /// Estimates the probability of the given cell in the underlying population. 
        /// </summary>
        /// <param name="r">The (zero-based) row index.</param>
        /// <param name="c">The (zero-based) column index.</param>
        /// <returns>The probability for a random event to occur inthe given cell.</returns>
        public virtual UncertainValue Probability (int r, int c) {
            if ((r < 0) || (r >= rows)) throw new ArgumentOutOfRangeException("r");
            if ((c < 0) || (c >= cols)) throw new ArgumentOutOfRangeException("c");
            double p = ((double) counts[r, c]) / tCounts;
            double dp = Math.Sqrt(p * (tCounts - counts[r, c])) / tCounts;
            return (new UncertainValue(p, dp));
        }

        /// <summary>
        /// Estimates the marginal probility of the given row in the underlying population.
        /// </summary>
        /// <param name="r">The (zero-based) row index.</param>
        /// <returns>The probability for a random event to occur in the given row.</returns>
        public virtual UncertainValue ProbabilityOfRow (int r) {
            if ((r < 0) || (r >= rows)) throw new ArgumentOutOfRangeException("r");
            double p = ((double) rCounts[r]) / tCounts;
            double dp = Math.Sqrt(p * (Total - rCounts[r])) / tCounts;
            return (new UncertainValue(p, dp));
        }

        /// <summary>
        /// Estimates the marginal probility of the given column in the underlying population.
        /// </summary>
        /// <param name="c">The (zero-based) column index.</param>
        /// <returns>The probability for a random envent to occur in the given column.</returns>
        public virtual UncertainValue ProbabilityOfColumn (int c) {
            if ((c < 0) || (c >= cols)) throw new ArgumentOutOfRangeException("c");
            double p = ((double) cCounts[c]) / tCounts;
            double dp = Math.Sqrt(p * (Total - cCounts[c])) / tCounts;
            return (new UncertainValue(p, dp));
        }

        /// <summary>
        /// Estimates the probability of the given row, conditional on the given column.
        /// </summary>
        /// <param name="r">The (zero-based) row index.</param>
        /// <param name="c">The (zero-based) column index.</param>
        /// <returns>The probability for a random event to occur in the given row, if it appears in the given column.</returns>
        public virtual UncertainValue ProbabilityOfRowConditionalOnColumn (int r, int c) {
            if ((r < 0) || (r >= rows)) throw new ArgumentOutOfRangeException("r");
            if ((c < 0) || (c >= cols)) throw new ArgumentOutOfRangeException("c");
            double p = ((double) counts[r, c]) / cCounts[c];
            double dp = Math.Sqrt(p * (cCounts[c] - counts[r, c])) / cCounts[c];
            return (new UncertainValue(p, dp));
        }

        /// <summary>
        /// Estimates the probability of the given row, conditional on the given column.
        /// </summary>
        /// <param name="r">The (zero-based) row index.</param>
        /// <param name="c">The (zero-based) column index.</param>
        /// <returns>The probability for a random event to occur in the given row, if it appears in the given column.</returns>
        public virtual UncertainValue ProbabilityOfColumnConditionalOnRow (int c, int r) {
            if ((c < 0) || (c >= cols)) throw new ArgumentOutOfRangeException("c");
            if ((r < 0) || (r >= rows)) throw new ArgumentOutOfRangeException("r");
            double p = ((double) counts[r, c]) / rCounts[r];
            double dp = Math.Sqrt(p * (rCounts[r] - counts[r, c])) / rCounts[r];
            return (new UncertainValue(p, dp));
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
        public virtual TestResult PearsonChiSquaredTest () {

            // compute chi squared
            double chi2 = 0.0;
            for (int r = 0; r < rows; r++) {
                for (int c = 0; c < cols; c++) {
                    double n = ((double) rCounts[r]) * cCounts[c] / tCounts;
                    double z = counts[r, c] - n;
                    chi2 += z * z / n;
                }
            }

            // the degrees of freedom is just the number of entries minus the number of constraints: R*C - R - (C - 1)
            int nu = (rows - 1) * (cols - 1);

            // return the test result
            return (new TestResult("ChiSquared", chi2, TestType.RightTailed, new ChiSquaredDistribution(nu)));

        }

#if FUTURE

        // name-based access methods

        public NameCollection RowNames {
            get {
                return (rowNames);
            }
        }

        public NameCollection ColumnNames {
            get {
                return (columnNames);
            }
        }

        public int this[string rName, string cName] {
            get {
                if (rName == null) throw new ArgumentNullException("rName");
                if (cName == null) throw new ArgumentNullException("cName");
                int r = rowNames.GetIndexForName(rName);
                if (r < 0) throw new InvalidOperationException();
                int c = columnNames.GetIndexForName(cName);
                if (c < 0) throw new InvalidOperationException();
                return (counts[r, c]);
            }
            set {
                if (rName == null) throw new ArgumentNullException("rName");
                if (cName == null) throw new ArgumentNullException("cName");
                int r = rowNames.GetIndexForName(rName);
                if (r < 0) throw new InvalidOperationException();
                int c = columnNames.GetIndexForName(cName);
                if (c < 0) throw new InvalidOperationException();
                if (value < 0) throw new ArgumentOutOfRangeException("value");
                counts[r, c] = value;
            }
        }

        /// <summary>
        /// Gets the total counts in the row with the given name.
        /// </summary>
        /// <param name="rName">The name of the row.</param>
        /// <returns>The sum of counts in the specified row.</returns>
        public int RowTotal (string rName) {
            if (rName == null) throw new ArgumentNullException("rName");
            int r = rowNames.GetIndexForName(rName);
            return (RowTotal(r));
        }

#endif

    }

#if FUTURE

    public class NameCollection {

        internal NameCollection (int capacity) {
            names = new string[capacity];
        }

        private string[] names;

        public int Count {
            get {
                return (names.Length);
            }
        }

        public string this [int index] {
            get {
                if ((index < 0) || (index >= Count)) throw new ArgumentOutOfRangeException("index");
                return (names[index]);
            }
            set {
                if ((index < 0) || (index >= Count)) throw new ArgumentOutOfRangeException("index");
                if ((value != null) && (GetIndexForName(value) >= 0)) throw new InvalidOperationException();
                names[index] = value;
            }
        }

        internal int GetIndexForName (string name) {
            int index = Array.IndexOf(names, name);
            return (index);
        }

    }

#endif

}
