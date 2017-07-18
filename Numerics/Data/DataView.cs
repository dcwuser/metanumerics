using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;

namespace Meta.Numerics.Data
{
    /// <summary>
    /// A readable array of data.
    /// </summary>
    public class DataView {

        internal DataView() {

        }

        internal DataView(List<DataList> columns, List<int> map) {
            this.columns = columns;
            this.columnMap = new Dictionary<string, int>();
            for (int columnIndex = 0; columnIndex < columns.Count; columnIndex++) {
                this.columnMap[columns[columnIndex].Name] = columnIndex;
            }
            this.map = map;
        }


        internal List<DataList> columns;
        internal Dictionary<string, int> columnMap;
        internal List<int> map;

        internal int GetColumnIndex (string name) {
            int columnIndex;
            if (this.columnMap.TryGetValue(name, out columnIndex)) {
                return (columnIndex);
            } else {
                return (-1);
            }
        }

        /// <summary>
        /// Gets the rows of the table.
        /// </summary>
        public DataRowCollection Rows {
            get {
                return (new DataRowCollection(this));
            }
        }

        /// <summary>
        /// Gets the columns of the table.
        /// </summary>
        public DataColumnCollection Columns {
            get {
                return (new DataColumnCollection(this));
            }
        }

        /// <summary>
        /// Get a given column, cast to a given data type.
        /// </summary>
        /// <typeparam name="T">The type into which to cast the column.</typeparam>
        /// <param name="columnName">The name of the column.</param>
        /// <returns>The requested column, cast to the given data type.</returns>
        public DataColumn<T> Column<T>(string columnName) {
            if (columnName == null) throw new ArgumentNullException(nameof(columnName));

            int c = GetColumnIndex(columnName);
            return (new DataColumn<T>(this, c));
        }

        /// <summary>
        /// Constructs a new view, containing only the given columns.
        /// </summary>
        /// <param name="columnNames">The names of the columns to include.</param>
        /// <returns>A view containing only the given columns.</returns>
        public DataView Select(params string[] columnNames) {
            return (Select((ICollection<string>)columnNames));
        }

        /// <summary>
        /// Constructs a new view, containg only the given columns.
        /// </summary>
        /// <param name="columnNames">The names of the columns to include.</param>
        /// <returns>A new view containing only the given columns.</returns>
        public DataView Select (ICollection<string> columnNames) {
            if (columnNames == null) throw new ArgumentNullException(nameof(columnNames));

            List<DataList> newColumns = new List<DataList>();
            foreach (string columnName in columnNames) {
                int columnIndex = GetColumnIndex(columnName);
                newColumns.Add(columns[columnIndex]);
            }
            return (new DataView(newColumns, map));
        }

        /// <summary>
        /// Sort the rows by the values in the given column.
        /// </summary>
        /// <param name="columnName">The name of the column to sort by.</param>
        /// <returns>A new view, with rows sorted in ascending order by the values in the given column.</returns>
        public DataView OrderBy(string columnName) {
            return (OrderBy(columnName, SortOrder.Ascending));
        }

        /// <summary>
        /// Sort the rows by the values in the given column in the given direction.
        /// </summary>
        /// <param name="columnName">The name of the column to sort by.</param>
        /// <param name="direction">The direction in which the result should be ordered.</param>
        /// <returns>A new view, with rows sorted in the given direction by the values in the given column.</returns>
        public DataView OrderBy(string columnName, SortOrder direction) {
            if (columnName == null) throw new ArgumentNullException(nameof(columnName));

            int columnIndex = GetColumnIndex(columnName);
            DataList column = columns[columnIndex];
            List<int> newMap = new List<int>(map);
            NullTollerantComparer comparer = new NullTollerantComparer(direction);
            newMap.Sort((i, j) => comparer.Compare((IComparable)column.GetItem(i), (IComparable)column.GetItem(j)));
            return (new DataView(this.columns, newMap));
        }


        /// <summary>
        /// Sorts all rows by the given function.
        /// </summary>
        /// <param name="comparer">A row comparison function.</param>
        /// <returns>A new view, with rows sorted by the given function.</returns>
        public DataView OrderBy(Comparison<DataRow> comparer) {
            if (comparer == null) throw new ArgumentNullException(nameof(comparer));

            List<int> newMap = new List<int>(map);
            newMap.Sort((i, j) => comparer(new DataRow(this, i), new DataRow(this, j)));
            return (new DataView(this.columns, newMap));
        }

        /// <summary>
        /// Selects rows matching the given criteria.
        /// </summary>
        /// <param name="selector">A function which accepts or rejects rows.</param>
        /// <returns>A new view, containing only those rows which were accepted by the <paramref name="selector"/> function.</returns>
        public DataView Where (Func<DataRow, bool> selector) {
            if (selector == null) throw new ArgumentNullException(nameof(selector));

            List<int> newMap = new List<int>();
            for (int i = 0; i < map.Count; i++) {
                DataRow row = new DataRow(this, i);
                bool include = selector(row);
                if (include) newMap.Add(map[i]);
            }
            return (new DataView(this.columns, newMap));
        }

        /// <summary>
        /// Selects rows in which the value of the given column matches the given criteria.
        /// </summary>
        /// <typeparam name="T">The type of the row data.</typeparam>
        /// <param name="columnName">The name of the column.</param>
        /// <param name="selector">A function which accepts or rejects values.</param>
        /// <returns>A new view, containing only those rows in which the <paramref name="columnName"/> value was accepted by
        /// the <paramref name="selector"/>.</returns>
        public DataView Where<T> (string columnName, Func<T, bool> selector) {
            if (columnName == null) throw new ArgumentNullException(nameof(columnName));
            if (selector == null) throw new ArgumentNullException(nameof(selector));

            int columnIndex = GetColumnIndex(columnName);
            DataList column = columns[columnIndex];
            //IReadOnlyDataList<T> column = (IReadOnlyDataList<T>) columns[columnIndex];

            List<int> newMap = new List<int>();
            for (int i = 0; i < map.Count; i++) {
                T value = (T) column.GetItem(map[i]); // column[map[i]];
                bool include = selector(value);
                if (include) newMap.Add(map[i]);
            }

            return (new DataView(this.columns, newMap));
        }

        /// <summary>
        /// Removes any rows that have null values in the named columns.
        /// </summary>
        /// <param name="columnNames">The columns to check for null values.</param>
        /// <returns>A new view, without the rows that have null values in the named columns.</returns>
        /// <remarks>
        /// <para>If no column names are given, rows with nulls in any column are removed.</para>
        /// </remarks>
        public DataView WhereNotNull(params string[] columnNames)
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

            return (new DataView(columns, newMap));
        }

        /// <summary>
        /// Groups the data by the values in the given column, and computes the given aggregate quantity for each group.
        /// </summary>
        /// <typeparam name="T">The type of the aggregate output.</typeparam>
        /// <param name="groupByColumnName">The name of the column to group by.</param>
        /// <param name="aggregator">A function that computes the aggregate quantity.</param>
        /// <param name="aggregateColumnName">The name of the column for the aggregate output.</param>
        /// <returns>A new data frame containing the requested aggregate values for each group.</returns>
        public DataFrame GroupBy<T>(string groupByColumnName, Func<DataView, T> aggregator, string aggregateColumnName)
        {
            if (groupByColumnName == null) throw new ArgumentNullException(nameof(groupByColumnName));
            if (aggregator == null) throw new ArgumentNullException(nameof(aggregator));
            if (aggregateColumnName == null) throw new ArgumentNullException(nameof(aggregateColumnName));

            int groupByColumnIndex = GetColumnIndex(groupByColumnName);
            DataList groupByColumn = columns[groupByColumnIndex];

            // Create a lookup table that maps group column values to the indexes of rows with that value.
            // We have to be a little tricky to deal with null values, because null is not allowed to be
            // a dictionary key. We use a special internal null signifier object to get around this problem.
            Dictionary<object, List<int>> groups = new Dictionary<object, List<int>>();
            for (int r = 0; r < map.Count; r++)
            {
                int index = map[r];
                object value = groupByColumn.GetItem(index);
                if (value == null) value = NullSignifier.Value;
                List<int> members;
                if (!groups.TryGetValue(value, out members))
                {
                    members = new List<int>();
                    groups.Add(value, members);
                }
                members.Add(index);
            }

            // Form destination columns based on group aggregates.
            DataList groupsColumn = DataList.Create(groupByColumnName, groupByColumn.StorageType);
            DataList<T> aggregateColumn = new DataList<T>(aggregateColumnName);
            foreach (KeyValuePair<object, List<int>> group in groups)
            {
                DataView values = new DataView(this.columns, group.Value);
                T aggregateValue = aggregator(values);
                aggregateColumn.AddItem(aggregateValue);

                object groupKey = group.Key;
                if (groupKey == NullSignifier.Value) groupKey = null;
                groupsColumn.AddItem(groupKey);

            }

            DataFrame result = new DataFrame(groupsColumn, aggregateColumn);
            return (result);
        }

        public DataFrame GroupBy<T, U>(string groupByColumnName, Func<DataView,Tuple<T,U>> aggregator, string firstAggregateColumnName, string secondAggregateColumnName) {
            throw new NotImplementedException();
        }


        /// <summary>
        /// Joins this data view to another data view on the given column.
        /// </summary>
        /// <param name="other">The second data view.</param>
        /// <param name="columnName">The column name, which must be the same in both views.</param>
        /// <returns>A new data frame whose columns are combined columns of the original views, inner joined on the non-null values of the given shared column.</returns>
        public DataFrame Join(DataView other, string columnName) {
            int thisColumnIndex = this.GetColumnIndex(columnName);
            DataList thisColumn = this.columns[thisColumnIndex];
            Type thisType = thisColumn.StorageType;
            int otherColumnIndex = other.GetColumnIndex(columnName);
            DataList otherColumn = other.columns[otherColumnIndex];
            Type otherType = otherColumn.StorageType;

            // Form a lookup from the other table
            Dictionary<object, int> hash = new Dictionary<object, int>();
            for (int otherRowIndex = 0; otherRowIndex < other.Rows.Count; otherRowIndex++) {
                object key = otherColumn.GetItem(other.map[otherRowIndex]);
                if (key == null) continue;
                hash[key] = otherRowIndex;
            }

            // Construct the joined columns
            List<DataList> joinedColumns = new List<DataList>();
            for (int i = 0; i < this.columns.Count; i++) {
                DataList joinedColumn = DataList.Create(this.columns[i].Name, this.columns[i].StorageType);
                joinedColumns.Add(joinedColumn);
            }
            for (int j = 0; j < other.columns.Count; j++) {
                //if (j == otherColumnIndex) continue;
                DataList joinedColumn = DataList.Create(other.columns[j].Name, other.columns[j].StorageType);
                joinedColumns.Add(joinedColumn);
            }

            // Populate the joined columns
            for (int thisRowIndex = 0; thisRowIndex < this.map.Count; thisRowIndex++) {
                object thisValue = thisColumn.GetItem(this.map[thisRowIndex]);
                if (thisValue == null) continue;
                int otherRowIndex;
                if (hash.TryGetValue(thisValue, out otherRowIndex)) {
                    for (int i = 0; i < this.columns.Count; i++) {
                        joinedColumns[i].AddItem(this.columns[i].GetItem(this.map[i]));
                    }
                    for (int j = 0; j < other.columns.Count; j++) {
                        //if (j == otherColumnIndex) continue;
                        joinedColumns[this.columns.Count + j].AddItem(other.columns[j].GetItem(other.map[otherRowIndex]));
                    }
                }
            }

            DataFrame result = new DataFrame(joinedColumns);
            return (result);

        }

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

        public void AddComputedColumn<T>(string columnName, Func<DataRow, T> function)
        {
            if (columnName == null) throw new ArgumentNullException(nameof(columnName));
            if (function == null) throw new ArgumentNullException(nameof(function));
            DataList column = new ComputedColumn<T>(this, columnName, function);
            int columnIndex = columns.Count;
            columns.Add(column);
            columnMap[columnName] = columnIndex;
        }

        // CSV output method

        /// <summary>
        /// Write the data in the view to a comma-seperated-value file.
        /// </summary>
        /// <param name="writer">A writer to accept the data.</param>
        public void WriteCsv (TextWriter writer) {
            if (writer == null) throw new ArgumentNullException(nameof(writer));

            string[] columnNames = new string[this.columns.Count];
            for (int c = 0; c < columns.Count; c++) {
                columnNames[c] = columns[c].Name;
            }
            writer.WriteLine(CsvHelper.WriteCells(columnNames));
            foreach (DataRow row in this.Rows) {
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
            foreach (DataRow row in this.Rows) {
                IReadOnlyDictionary<string, object> rowAsDictionary = (IReadOnlyDictionary<string, object>) row;
                Dictionary<string, object> rowDictionary = new Dictionary<string, object>();
                foreach (KeyValuePair<string, object> entry in rowAsDictionary) {
                    rowDictionary.Add(entry.Key, entry.Value);
                }
                yield return (rowDictionary);
            }
        }

    }

    internal sealed class NullSignifier {

        private NullSignifier() { }

        public static readonly NullSignifier Value = new NullSignifier();

        public override int GetHashCode () {
            return (-12948583);
        }

    }

    internal class NullTollerantComparer {

        public NullTollerantComparer(SortOrder direction) {
            if (direction == SortOrder.Ascending) {
                ordered = +1;
            } else {
                ordered = -1;
            }
        }

        private readonly int ordered;

        public int Compare(IComparable a, IComparable b) {
            if (a == null) {
                if (b == null) {
                    return (0);
                }
                else {
                    return (-ordered);
                }
            }
            else {
                if (b == null) {
                    return (+ordered);
                }
                else {
                    return (a.CompareTo(b));
                }
            }
        }

    }

    /// <summary>
    /// Specifies the order of a sort.
    /// </summary>
    public enum SortOrder {

        /// <summary>
        /// Ascending sort order.
        /// </summary>
        Ascending,

        /// <summary>
        /// Descending sort order.
        /// </summary>
        Descending
    }

}
