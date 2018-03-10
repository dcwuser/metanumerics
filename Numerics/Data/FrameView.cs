using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

using Meta.Numerics.Statistics;

namespace Meta.Numerics.Data
{

    /// <summary>
    /// Specifies the order of sorting.
    /// </summary>
    public enum SortOrder {

        /// <summary>
        /// From smallest to largest.
        /// </summary>
        Ascending,

        /// <summary>
        /// From largest to smallest.
        /// </summary>
        Descending
    }


    /// <summary>
    /// A read-only view of a table of data.
    /// </summary>
    public class FrameView {

        internal FrameView(List<NamedList> columns, List<int> map) {
            this.columns = columns;
            this.columnMap = new Dictionary<string, int>();
            for (int columnIndex = 0; columnIndex < columns.Count; columnIndex++) {
                this.columnMap[columns[columnIndex].Name] = columnIndex;
            }
            this.map = map;
        }


        internal readonly List<NamedList> columns;
        internal readonly Dictionary<string, int> columnMap;
        internal readonly List<int> map;

        /// <summary>
        /// Gets the index of a given column.
        /// </summary>
        /// <param name="columnName">The name of the column.</param>
        /// <returns>The (zero-based) index of the column, or -1 if there is no such column in the table.</returns>
        public int GetColumnIndex (string columnName) {
            if (columnName == null) throw new ArgumentNullException(nameof(columnName));
            int columnIndex;
            if (this.columnMap.TryGetValue(columnName, out columnIndex)) {
                return (columnIndex);
            } else {
                return (-1);
            }
        }

        /// <summary>
        /// Gets the value in the given cell.
        /// </summary>
        /// <param name="rowIndex">The (zero-based) row index of the cell.</param>
        /// <param name="columnIndex">The (zero-based) column index of the cell.</param>
        /// <returns></returns>
        public object this [int rowIndex, int columnIndex] {
            get {
                return (columns[columnIndex].GetItem(map[rowIndex]));
            }
        }

        /// <summary>
        /// Gets the rows of the table.
        /// </summary>
        public FrameRowCollection Rows {
            get {
                return (new FrameRowCollection(this));
            }
        }

        /// <summary>
        /// Gets the columns of the table.
        /// </summary>
        public FrameColumnCollection Columns {
            get {
                return (new FrameColumnCollection(this));
            }
        }

        /// <summary>
        /// Get the column with the given name.
        /// </summary>
        /// <param name="columnName">The name of the column.</param>
        /// <returns>The requested column.</returns>
        public FrameColumn this [string columnName] {
            get {
                if (columnName == null) throw new ArgumentNullException(nameof(columnName));
                int c = GetColumnIndex(columnName);
                return (new FrameColumn(this, c));
            }
        }

        /// <summary>
        /// Constructs a new view, containing only the given columns.
        /// </summary>
        /// <param name="columnNames">The names of the columns to include.</param>
        /// <returns>A view containing only the given columns.</returns>
        public FrameView Select(params string[] columnNames) {
            return (Select((ICollection<string>)columnNames));
        }

        /// <summary>
        /// Constructs a new view, containg only the given columns.
        /// </summary>
        /// <param name="columnNames">The names of the columns to include.</param>
        /// <returns>A new view containing only the given columns.</returns>
        public FrameView Select (ICollection<string> columnNames) {
            if (columnNames == null) throw new ArgumentNullException(nameof(columnNames));

            List<NamedList> newColumns = new List<NamedList>();
            foreach (string columnName in columnNames) {
                int columnIndex = GetColumnIndex(columnName);
                newColumns.Add(columns[columnIndex]);
            }
            return (new FrameView(newColumns, map));
        }

        /// <summary>
        /// Constructs a new view that does not contain the given columns.
        /// </summary>
        /// <param name="columnNames">The names of the columns to discard.</param>
        /// <returns>A new view that does not contain the given columns.</returns>
        public FrameView Discard(params string[] columnNames) {
            if (columnNames == null) throw new ArgumentNullException(nameof(columnNames));

            HashSet<string> names = new HashSet<string>(columnNames);

            List<NamedList> newColumns = new List<NamedList>();
            foreach(NamedList column in columns) {
                if (!names.Contains(column.Name)) newColumns.Add(column);
            }
            return (new FrameView(newColumns, map));
        }

        /// <summary>
        /// Sorts the rows by the values in the given column.
        /// </summary>
        /// <param name="columnName">The name of the column to sort by.</param>
        /// <returns>A new view, with rows sorted by the values in the given column.</returns>
        public FrameView OrderBy(string columnName) {
            return (OrderBy(columnName, SortOrder.Ascending));
        }

        /// <summary>
        /// Sort the rows by the values in the given column in the given direction.
        /// </summary>
        /// <param name="columnName">The name of the column to sort by.</param>
        /// <param name="order">The direction of the ordering.</param>
        /// <returns>A new view, with rows sorted by the values in the given column.</returns>
        /// <remarks>
        /// <para>The storage type of the column must implement <see cref="IComparable"/>.</para>
        /// </remarks>
        public FrameView OrderBy(string columnName, SortOrder order) {
            if (columnName == null) throw new ArgumentNullException(nameof(columnName));

            int columnIndex = GetColumnIndex(columnName);
            NamedList column = columns[columnIndex];
            List<int> newMap = new List<int>(map);
            if (order == SortOrder.Ascending) {
                newMap.Sort((i, j) => NullableComparer((IComparable) column.GetItem(i), (IComparable) column.GetItem(j)));
            } else {
                newMap.Sort((i, j) => NullableComparer((IComparable) column.GetItem(j), (IComparable) column.GetItem(i)));
            }
            return (new FrameView(this.columns, newMap));
        }

        // This comparer is able to deal with null values.

        private static int NullableComparer (IComparable a, IComparable b) {
            if (a == null) {
                if (b == null) {
                    return (0);
                } else {
                    return (-1);
                }
            } else {
                if (b == null) {
                    return (+1);
                } else {
                    return (a.CompareTo(b));
                }
            }
        }

        /// <summary>
        /// Sorts all rows by the given function of the given column.
        /// </summary>
        /// <typeparam name="T">The type of values being compared.</typeparam>
        /// <param name="columnName">The column to sort on.</param>
        /// <param name="comparer">A comparison function of <paramref name="columnName"/> values.</param>
        /// <returns>A view of the data sorted by the function of the given column.</returns>
        public FrameView OrderBy<T>(string columnName, Comparison<T> comparer) {
            if (columnName == null) throw new ArgumentNullException(nameof(columnName));
            if (comparer == null) throw new ArgumentNullException(nameof(comparer));

            int columnIndex = GetColumnIndex(columnName);
            IReadOnlyList<T> column = (IReadOnlyList<T>)columns[columnIndex];

            List<int> newMap = new List<int>(map);
            newMap.Sort((i, j) => comparer(column[i], column[j]));

            return (new FrameView(this.columns, newMap));
        }

        /// <summary>
        /// Sorts all rows by the given function.
        /// </summary>
        /// <param name="comparer">A row comparison function.</param>
        /// <returns>A new view, with rows sorted by the given function.</returns>
        public FrameView OrderBy(Comparison<FrameRow> comparer) {
            if (comparer == null) throw new ArgumentNullException(nameof(comparer));

            List<int> newMap = new List<int>(map);
            newMap.Sort((i, j) => comparer(new FrameRow(this, i), new FrameRow(this, j)));
            return (new FrameView(this.columns, newMap));
        }

        /// <summary>
        /// Selects rows matching the given criteria.
        /// </summary>
        /// <param name="selector">A function which accepts or rejects rows.</param>
        /// <returns>A new view, containing only those rows which were accepted by the <paramref name="selector"/> function.</returns>
        public FrameView Where (Func<FrameRow, bool> selector) {
            if (selector == null) throw new ArgumentNullException(nameof(selector));

            List<int> newMap = new List<int>();
            for (int i = 0; i < map.Count; i++) {
                FrameRow row = new FrameRow(this, i);
                bool include = selector(row);
                if (include) newMap.Add(map[i]);
            }
            return (new FrameView(this.columns, newMap));
        }

        /// <summary>
        /// Selects rows in which the value of the given column matches the given criteria.
        /// </summary>
        /// <typeparam name="T">The type of the row data.</typeparam>
        /// <param name="columnName">The name of the column.</param>
        /// <param name="selector">A function which accepts or rejects values.</param>
        /// <returns>A new view, containing only those rows in which the <paramref name="columnName"/> value was accepted by
        /// the <paramref name="selector"/>.</returns>
        public FrameView Where<T> (string columnName, Func<T, bool> selector) {
            if (columnName == null) throw new ArgumentNullException(nameof(columnName));
            if (selector == null) throw new ArgumentNullException(nameof(selector));

            int columnIndex = GetColumnIndex(columnName);
            NamedList<T> column = (NamedList<T>) columns[columnIndex];

            List<int> newMap = new List<int>();
            for (int i = 0; i < map.Count; i++) {
                T value = column[map[i]];
                bool include = selector(value);
                if (include) newMap.Add(map[i]);
            }

            return (new FrameView(this.columns, newMap));
        }

        /// <summary>
        /// Removes any rows that have null values in the named columns.
        /// </summary>
        /// <param name="columnNames">The columns to check for null values.</param>
        /// <returns>A new view, without the rows that have null values in the named columns.</returns>
        /// <remarks>
        /// <para>If no column names are given, rows with nulls in any column are removed.</para>
        /// </remarks>
        public FrameView WhereNotNull(params string[] columnNames)
        {
            if (columnNames == null) throw new ArgumentNullException(nameof(columnNames));

            // Convert the list of names into a list of indexes to check.
            // If no names are given, check all columns.
            // There is no need to check non-nullable columns, so leave them out.
            List<int> columnIndexes = new List<int>();
            if (columnNames.Length == 0) {
                for (int columnIndex = 0; columnIndex < columns.Count; columnIndex++) {
                    if (columns[columnIndex].IsNullable) columnIndexes.Add(columnIndex);
                }
            } else {
                foreach (string columnName in columnNames) {
                    int columnIndex = GetColumnIndex(columnName);
                    if (columns[columnIndex].IsNullable) columnIndexes.Add(columnIndex);
                }
            }

            // Go through each row, checking for nulls in the noted columns.
            // If no nulls are found, add the row to the new map.
            List<int> newMap = new List<int>();
            foreach (int rowIndex in map)
            {
                bool insert = true;
                foreach (int columnIndex in columnIndexes)
                {
                    object value = columns[columnIndex].GetItem(rowIndex);
                    if (value == null)
                    {
                        insert = false;
                        break;
                    }
                }
                if (insert) newMap.Add(rowIndex);
            }

            return (new FrameView(columns, newMap));
        }

        /// <summary>
        /// Groups the data by the values in the given column, and computes the given aggregate quantity for each group.
        /// </summary>
        /// <typeparam name="T">The type of the aggregate output.</typeparam>
        /// <param name="groupByColumnName">The name of the column to group by.</param>
        /// <param name="aggregator">A function that computes the aggregate quantity.</param>
        /// <param name="aggregateColumnName">The name of the column for the aggregate output.</param>
        /// <returns>A new data frame containing the requested aggregate values for each group.</returns>
        public FrameTable GroupBy<T>(string groupByColumnName, Func<FrameView, T> aggregator, string aggregateColumnName)
        {
            if (groupByColumnName == null) throw new ArgumentNullException(nameof(groupByColumnName));
            if (aggregator == null) throw new ArgumentNullException(nameof(aggregator));
            if (aggregateColumnName == null) throw new ArgumentNullException(nameof(aggregateColumnName));

            int groupByColumnIndex = GetColumnIndex(groupByColumnName);
            NamedList groupByColumn = columns[groupByColumnIndex];

            // Create a lookup table that maps group column values to the indexes of rows with that value.
            // We have to be a little tricky to deal with null values, because null is not allowed to be
            // a dictionary key. We use a special internal null signifier object to get around this problem.
            //Dictionary<object, List<int>> groups = new Dictionary<object, List<int>>();
            NullableDictionary<object, List<int>> groups = new NullableDictionary<object, List<int>>();
            for (int r = 0; r < map.Count; r++)
            {
                int index = map[r];
                object value = groupByColumn.GetItem(index);
                //if (value == null) value = NullSignifier.Value;
                List<int> members;
                if (!groups.TryGetValue(value, out members))
                {
                    members = new List<int>();
                    groups.Add(value, members);
                }
                members.Add(index);
            }

            // Form destination columns based on group aggregates.
            NamedList groupsColumn = NamedList.Create(groupByColumnName, groupByColumn.StorageType);
            NamedList<T> aggregateColumn = new NamedList<T>(aggregateColumnName);
            foreach (KeyValuePair<object, List<int>> group in groups)
            {
                FrameView values = new FrameView(this.columns, group.Value);
                T aggregateValue = aggregator(values);
                aggregateColumn.AddItem(aggregateValue);

                object groupKey = group.Key;
                //if (groupKey == NullSignifier.Value) groupKey = null;
                groupsColumn.AddItem(groupKey);

            }

            FrameTable result = new FrameTable(groupsColumn, aggregateColumn);
            return (result);
        }

        /// <summary>
        /// Groups the data by the values in the given column, and computes the given aggregate quantities for each group.
        /// </summary>
        /// <param name="groupByColumnName">The name of the column to group by.</param>
        /// <param name="aggregateColumns">The definitions of the aggregate columns to compute.</param>
        /// <returns>A new data frame containing the requested aggregate values for each group.</returns>
        public FrameTable GroupBy (string groupByColumnName, params AggregateDefinition[] aggregateColumns) {

            if (groupByColumnName == null) throw new ArgumentNullException(nameof(groupByColumnName));
            if (aggregateColumns == null) throw new ArgumentNullException(nameof(aggregateColumns));

            int groupByColumnIndex = GetColumnIndex(groupByColumnName);
            NamedList groupByColumn = columns[groupByColumnIndex];

            // Create a lookup table that maps group column values to the indexes of rows with that value.
            // We have to be a little tricky to deal with null values, because null is not allowed to be
            // a dictionary key. We use a special internal null signifier object to get around this problem.
            //Dictionary<object, List<int>> groups = new Dictionary<object, List<int>>();
            NullableDictionary<object, List<int>> groups = new NullableDictionary<object, List<int>>();
            for (int r = 0; r < map.Count; r++) {
                int index = map[r];
                object value = groupByColumn.GetItem(index);
                //if (value == null) value = NullSignifier.Value;
                List<int> members;
                if (!groups.TryGetValue(value, out members)) {
                    members = new List<int>();
                    groups.Add(value, members);
                }
                members.Add(index);
            }

            // Set up result columns
            List<NamedList> resultColumns = new List<NamedList>();
            NamedList groupsColumn = NamedList.Create(groupByColumnName, groupByColumn.StorageType);
            resultColumns.Add(groupsColumn);
            foreach (AggregateDefinition aggregateColumn in aggregateColumns) {
                if (aggregateColumn == null) throw new ArgumentNullException(nameof(aggregateColumn));
                NamedList outputColumn = NamedList.Create(aggregateColumn.Name, aggregateColumn.StorageType);
                resultColumns.Add(outputColumn);
            }

            // Form destination columns based on group aggregates.
            foreach (KeyValuePair<object, List<int>> group in groups) {

                // Add the value for the grouping column
                object groupValue = group.Key;
                //if (groupValue == NullSignifier.Value) groupValue = null;
                groupsColumn.AddItem(groupValue);

                // Create a filtered view of the data
                FrameView inputView = new FrameView(this.columns, group.Value);

                // Apply each aggregator to the filtered view to get the values to add to the aggregate columns.
                for (int i = 0; i < aggregateColumns.Length; i++) {
                    AggregateDefinition aggregateColumn = aggregateColumns[i];
                    object outputValue = aggregateColumn.ApplyAggregator(inputView);
                    resultColumns[i + 1].AddItem(outputValue);
                }
            }

            FrameTable result = new FrameTable(resultColumns);
            return (result);
        }

        /*
        public DataFrame Pivot(string rowNamesColumnName, string columnNamesColumnName, string valuesColumnName)
        {
            if (rowNamesColumnName == null) throw new ArgumentNullException(nameof(rowNamesColumnName));
            if (columnNamesColumnName == null) throw new ArgumentNullException(nameof(columnNamesColumnName));
            if (valuesColumnName == null) throw new ArgumentNullException(nameof(valuesColumnName));

            int rowNamesColumnIndex = GetColumnIndex(rowNamesColumnName);
            int columnNamesColumnIndex = GetColumnIndex(columnNamesColumnName);
            int valuesColumnIndex = GetColumnIndex(valuesColumnName);

            DataList rowNamesColumn = columns[rowNamesColumnIndex];
            DataList columnNamesColumn = columns[columnNamesColumnIndex];
            DataList valuesColumn = columns[valuesColumnIndex];

            // Create a list of row names, column names, and a dictionary of corresponding values
            HashSet<object> rowNames = new HashSet<object>();
            HashSet<object> columnNames = new HashSet<object>();
            Dictionary<Tuple<object, object>, object> values = new Dictionary<Tuple<object, object>, object>();
            for (int r = 0; r < map.Count; r++)
            {
                int i = map[r];
                object rowName = rowNamesColumn.GetItem(i);
                object columnName = columnNamesColumn.GetItem(i);
                object vValue = valuesColumn.GetItem(i);
                rowNames.Add(rowName);
                columnNames.Add(columnName);
                Tuple<object, object> key = Tuple.Create(rowName, columnName);
                values.Add(key, vValue);
            }

            DataList keyColumn = DataList.Create(rowNamesColumn.Name, rowNamesColumn.StorageType);
            foreach (object rowName in rowNames)
            {
                keyColumn.AddItem(rowName);
            }

            // The new columns have the same type as the original, unless we need them to be nullable.
            Type valueType = valuesColumn.StorageType;
            if ((valueType.GetTypeInfo().IsValueType) && (values.Count < rowNames.Count * columnNames.Count) && (Nullable.GetUnderlyingType(valueType) == null))
            {
                valueType = typeof(Nullable<>).MakeGenericType(valueType);
            }

            // if values.Count < rowValues.Count * columnValues.Count, make nullable

            List<DataList> outputColumns = new List<DataList>();
            outputColumns.Add(keyColumn);

            foreach (object columnName in columnNames)
            {
                DataList outputColumn = DataList.Create(columnName.ToString(), valueType);
                for (int r = 0; r < keyColumn.Count; r++)
                {
                    object rowValue = keyColumn.GetItem(r);
                    Tuple<object, object> key = Tuple.Create(rowValue, columnName);
                    object value;
                    if (values.TryGetValue(key, out value))
                    {
                        outputColumn.AddItem(value);
                    }
                    else
                    {
                        outputColumn.AddItem(null);
                    }
                }

                outputColumns.Add(outputColumn);
            }


            return (new DataFrame(outputColumns));
        }
        */
        
        /// <summary>
        /// Add a computed column.
        /// </summary>
        /// <typeparam name="T">The type of the computed value.</typeparam>
        /// <param name="columnName">The name of the computed column.</param>
        /// <param name="function">The function that computes the column value.</param>
        public void AddComputedColumn<T>(string columnName, Func<FrameRow, T> function)
        {
            if (columnName == null) throw new ArgumentNullException(nameof(columnName));
            if (function == null) throw new ArgumentNullException(nameof(function));
            NamedList column = new ComputedDataList<T>(this, columnName, function);
            int columnIndex = columns.Count;
            columns.Add(column);
            columnMap[columnName] = columnIndex;
        }

        // CSV output method

        /// <summary>
        /// Write the data in the view to a comma-seperated-value file.
        /// </summary>
        /// <param name="writer">A writer to accept the data.</param>
        public void ToCsv (TextWriter writer) {
            if (writer == null) throw new ArgumentNullException(nameof(writer));

            string[] columnNames = new string[this.columns.Count];
            for (int c = 0; c < columns.Count; c++) {
                columnNames[c] = columns[c].Name;
            }
            writer.WriteLine(CsvHelper.WriteCells(columnNames));
            foreach (FrameRow row in this.Rows) {
                writer.WriteLine(CsvHelper.WriteCells(row));
            }
        }

        // JSON output method

        /// <summary>
        /// Writes the data in the view as a sequence of dictionaries.
        /// </summary>
        /// <returns>An enumerable sequence of dictionaries, one for each row.</returns>
        /// <remarks>
        /// <para>Each row produces one dictionary, whose keys are the column names and whose values
        /// are the corresponding column values in that row.</para>
        /// <para>This method can be used to produce a JSON seriealized form of the data.</para>
        /// </remarks>
        public IEnumerable<Dictionary<string, object>> ToDictionaries () {
            foreach (FrameRow row in this.Rows) {
                yield return (row.ToDictionary());
            }
        }

    }

}
