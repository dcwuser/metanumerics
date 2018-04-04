using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a contingency table.
    /// </summary>
    /// <typeparam name="R">The type of the row labels.</typeparam>
    /// <typeparam name="C">The type of the column labels.</typeparam>
    /// <remmarks>
    /// <para>Contingency tables summarize a bivariate set of categorical data by giving the total number
    /// of occurrences of every unique combination of the two variables. In many cases, this summary
    /// data is all that is needed to perform the required analysis.</para>
    /// <para>Given two columns of categorical data, you can construct a contingency table using the
    /// <see cref="Bivariate.Crosstabs{R, C}(IReadOnlyList{R}, IReadOnlyList{C})"/> method. You can
    /// also construct a contingency table explicitly and then use <see cref="this[R, C]"/> or
    /// <see cref="Increment(R, C)"/> to set the count in each cell.</para>
    /// </remmarks>
    /// <seealso href="http://en.wikipedia.org/wiki/Contingency_table" />
    public class ContingencyTable<R, C> {

        /// <summary>
        /// Initializes a new contingency table.
        /// </summary>
        /// <param name="rowValues">The low labels.</param>
        /// <param name="columnValues">The column labels.</param>
        public ContingencyTable (IReadOnlyCollection<R> rowValues, IReadOnlyCollection<C> columnValues) {

            if (rowValues == null) throw new ArgumentNullException(nameof(rowValues));
            if (columnValues == null) throw new ArgumentNullException(nameof(columnValues));

            rowMap = new NullableDictionary<R, int>();
            foreach (R rowValue in rowValues) rowMap.Add(rowValue, rowMap.Count);

            columnMap = new NullableDictionary<C, int>();
            foreach (C columnValue in columnValues) columnMap.Add(columnValue, columnMap.Count);

            counts = new int[rowMap.Count, columnMap.Count];
            rCounts = new int[rowMap.Count];
            cCounts = new int[columnMap.Count];
            tCounts = 0;
        }

        internal ContingencyTable (NullableDictionary<R, int> rowMap, NullableDictionary<C, int> columnMap, int[,] counts, int[] rCounts, int[] cCounts, int tCounts) {
            this.rowMap = rowMap;
            this.columnMap = columnMap;
            this.counts = counts;
            this.rCounts = rCounts;
            this.cCounts = cCounts;
            this.tCounts = tCounts;
        }

        private readonly NullableDictionary<R, int> rowMap;
        private readonly NullableDictionary<C, int> columnMap;

        private readonly int[,] counts;
        private readonly int[] rCounts;
        private readonly int[] cCounts;
        private int tCounts;

        /// <summary>
        /// Gets the number of events that occurred in the specified cell.
        /// </summary>
        /// <param name="row">The row value.</param>
        /// <param name="column">The column value.</param>
        /// <returns></returns>
        public virtual int this[R row, C column] {
            get {
                int rowIndex = rowMap[row];
                int columnIndex = columnMap[column];
                return (counts[rowIndex, columnIndex]);
            }
            set {
                if (value < 0) throw new ArgumentOutOfRangeException(nameof(value));
                int rowIndex = rowMap[row];
                int columnIndex = columnMap[column];
                int delta = value - counts[rowIndex, columnIndex];
                counts[rowIndex, columnIndex] = value;
                rCounts[rowIndex] += delta;
                Debug.Assert(rCounts[rowIndex] >= 0);
                cCounts[columnIndex] += delta;
                Debug.Assert(cCounts[columnIndex] >= 0);
                tCounts += delta;
                Debug.Assert(tCounts >= 0);
            }
        }

        /// <summary>
        /// Increases the count in the specified cell by one.
        /// </summary>
        /// <param name="row">The row value.</param>
        /// <param name="column">The column value.</param>
        public virtual void Increment (R row, C column) {
            int rowIndex = rowMap[row];
            int columnIndex = columnMap[column];
            counts[rowIndex, columnIndex]++;
            rCounts[rowIndex]++;
            cCounts[columnIndex]++;
            tCounts++;
        }

        /// <summary>
        /// Decreases the count in the specified cell by one.
        /// </summary>
        /// <param name="row">The row value.</param>
        /// <param name="column">The column value.</param>
        public virtual void Decrement (R row, C column) {
            int rowIndex = rowMap[row];
            int columnIndex = columnMap[column];
            if (counts[rowIndex, columnIndex] < 1) throw new InvalidOperationException();  
            counts[rowIndex, columnIndex]--;
            rCounts[rowIndex]--;
            Debug.Assert(rCounts[rowIndex] >= 0);
            cCounts[columnIndex]--;
            Debug.Assert(cCounts[columnIndex] >= 0);
            tCounts--;
            Debug.Assert(tCounts >= 0);
        }

        /// <summary>
        /// Gets the total number of events that occurred in the specified row.
        /// </summary>
        /// <param name="row">The row value.</param>
        /// <returns>The total number of counts counts in all cells in the specified row.</returns>
        public virtual int RowTotal (R row) {
            int rowIndex = rowMap[row];
            return (rCounts[rowIndex]);
        }

        /// <summary>
        /// Gets the total number of events that occurred in the specified column.
        /// </summary>
        /// <param name="column">The column value.</param>
        /// <returns>The total number of counts counts in all cells in the specified column.</returns>
        public virtual int ColumnTotal (C column) {
            int columnIndex = columnMap[column];
            return (cCounts[columnIndex]);
        }

        /// <summary>
        /// Gets the total number of events recorded in the table.
        /// </summary>
        public virtual int Total {
            get {
                return (tCounts);
            }
        }

        /// <summary>
        /// Gets the row values.
        /// </summary>
        public virtual IReadOnlyCollection<R> Rows {
            get {
                return (rowMap.Keys);
            }
        }

        /// <summary>
        /// Gets the column values.
        /// </summary>
        public virtual IReadOnlyCollection<C> Columns {
            get {
                return (columnMap.Keys);
            }
        }

        // row and column probabilities are the fraction of total counts that a given row/column total represents
        // conditional probabilities are the fraction of row/column totals that a cell represents
        // to get the uncertainties in these probabilities, use
        //   (df)^2 = (df/da)^2 (da)^2 + (df/db)^2 (db)^2 + ...
        // and use (dn)^2 = n for each cell count

        /// <summary>
        /// Estimates the probability of the given cell in the underlying population. 
        /// </summary>
        /// <param name="row">The row value.</param>
        /// <param name="column">The column value.</param>
        /// <returns>An estimate, with uncertainty, of the probability for a random event to occur in the given cell.</returns>
        public virtual UncertainValue ProbabilityOf (R row, C column) {
            int r = rowMap[row];
            int c = columnMap[column];
            double p = ((double) counts[r, c]) / tCounts;
            Debug.Assert((0.0 <= p) && (p <= 1.0));
            double dp = Math.Sqrt(p * (tCounts - counts[r, c])) / tCounts;
            return (new UncertainValue(p, dp));
        }

        /// <summary>
        /// Estimates the probability of the given row in the underlying population.
        /// </summary>
        /// <param name="row">The row value.</param>
        /// <returns>An estimate, with uncertainty, of the probability for a random event to occur in the given row.</returns>
        public virtual UncertainValue ProbabilityOfRow (R row) {
            int r = rowMap[row];
            double p = ((double) rCounts[r]) / tCounts;
            Debug.Assert((0.0 <= p) && (p <= 1.0));
            double dp = Math.Sqrt(p * (tCounts - rCounts[r])) / tCounts;
            return (new UncertainValue(p, dp));
        }

        /// <summary>
        /// Estimates the probability of the given row, conditional on the given column, in the underlying population.
        /// </summary>
        /// <param name="row">The row value.</param>
        /// <param name="column">The column value.</param>
        /// <returns>An estimate, with uncertainty, of the probability for a random event to occur in the given row, if it occurs in the given column.</returns>
        public virtual UncertainValue ProbabilityOfRowConditionalOnColumn (R row, C column) {
            int r = rowMap[row];
            int c = columnMap[column];
            double p = ((double) counts[r, c]) / cCounts[c];
            Debug.Assert((0.0 <= p) && (p <= 1.0));
            double dp = Math.Sqrt(p * (cCounts[c] - counts[r, c])) / cCounts[c];
            return (new UncertainValue(p, dp));
        }

        /// <summary>
        /// Estimates the probability of the given column in the underlying population.
        /// </summary>
        /// <param name="column">The column value.</param>
        /// <returns>An estimate, with uncertainty, of the probability for a random event to occur in the given column.</returns>
        public virtual UncertainValue ProbabilityOfColumn (C column) {
            int c = columnMap[column];
            double p = ((double) cCounts[c]) / tCounts;
            Debug.Assert((0.0 <= p) && (p <= 1.0));
            double dp = Math.Sqrt(p * (tCounts - cCounts[c])) / tCounts;
            return (new UncertainValue(p, dp));
        }

        /// <summary>
        /// Estimates the probability of the given column, conditional on the given row, in the underlying population.
        /// </summary>
        /// <param name="column">The column value.</param>
        /// <param name="row">The row value.</param>
        /// <returns>An estimate, with uncertainty, of the probability for a random event to occur in the given column, if it occurs in the given row.</returns>
        public virtual UncertainValue ProbabilityOfColumnConditionalOnRow (C column, R row) {
            int r = rowMap[row];
            int c = columnMap[column];
            double p = ((double) counts[r, c]) / rCounts[r];
            Debug.Assert((0.0 <= p) && (p <= 1.0));
            double dp = Math.Sqrt(p * (rCounts[r] - counts[r, c])) / rCounts[r];
            return (new UncertainValue(p, dp));
        }

        /// <summary>
        /// Performs a Pearson &#x3C7;<sup>2</sup> test for correlation in the table.
        /// </summary>
        /// <returns>The result of the test. The test statistic is &#x3C7;<sup>2</sup> and its likelihood under the null hypothesis is
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
        /// <see cref="BinaryContingencyTableOperations.FisherExactTest" /> is a viable alternative in these cases.</para></remarks>
        /// <seealso cref="ChiSquaredDistribution"/>
        public virtual TestResult PearsonChiSquaredTest () {

            // suppose there are no correlations present; then the expected number of events in a given table entry is just the fraction
            // of the total events proportional to the corresponding row and column totals: N_RC = N * (N_R / N) * (N_C / N) = N_R * N_C / N
            // if this number is large for all entries, then the entries should be distributed normally with mean N_RC and variance \sqrt(N_RC)
            // we can test this hypothesis by measuring the departure of the actual entries from the expected entries using a chi^2 test

            // compute chi squared
            double chi2 = 0.0;
            for (int r = 0; r < rCounts.Length; r++) {
                for (int c = 0; c < cCounts.Length; c++) {
                    double n = ((double) rCounts[r]) * cCounts[c] / tCounts;
                    double z = counts[r, c] - n;
                    chi2 += z * z / n;
                }
            }

            // the degrees of freedom is just the number of entries minus the number of constraints: R*C - R - (C - 1)
            int nu = (rCounts.Length - 1) * (cCounts.Length - 1);

            // return the test result
            return (new TestResult("χ²", chi2, TestType.RightTailed, new ChiSquaredDistribution(nu)));

        }

        /// <summary>
        /// Gets a façade that exposes properties defined only for 2 X 2 contingency tables.
        /// </summary>
        /// <exception cref="InvalidOperationException">The contingency table is not 2 X 2.</exception>
        public BinaryContingencyTableOperations Binary {
            get {
                if (rCounts.Length != 2 || cCounts.Length != 2) throw new InvalidOperationException();
                return (new BinaryContingencyTableOperations(counts));
            }
        }

    }

    /// <summary>
    /// Represents a contingency table without row and column labels.
    /// </summary>
    /// <remarks>
    /// <para>If your row and column variables represent known values, your code will be made clearer by
    /// using the <see cref="ContingencyTable{R, C}"/> class.</para>
    /// </remarks>
    public class ContingencyTable : ContingencyTable<int, int> {

        /// <summary>
        /// Initializes a new contingency table.
        /// </summary>
        /// <param name="rows">The number of rows, which must be positive.</param>
        /// <param name="columns">The number of columns, which must be positive.</param>
        public ContingencyTable (int rows, int columns) : base(GetIntegers(rows), GetIntegers(columns)) {

        }

        private static int[] GetIntegers (int n) {
            int[] array = new int[n];
            for (int i = 0; i < n; i++) array[i] = i;
            return (array);
        }

        /// <summary>
        /// Initializes a new contingency table with the given data set.
        /// </summary>
        /// <param name="data">A (zero-based) matrix of contingency table entries.</param>
        public ContingencyTable (int[,] data) : base(GetIntegers(data.GetLength(0)), GetIntegers(data.GetLength(1)))  {
            if (data == null) throw new ArgumentNullException(nameof(data));
            for (int r = 0; r < data.GetLength(0); r++) {
                for (int c = 0; c < data.GetLength(1); c++) {
                    if (data[r, c] < 0) throw new ArgumentOutOfRangeException(nameof(data));
                    this[r, c] = data[r, c];
                }
            }
        }

    }



}
