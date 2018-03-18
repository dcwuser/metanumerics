using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

namespace Meta.Numerics.Data
{

    // Supported storage types: int, double, bool, datetime, string
    // Supported data types: continuous, count, ordinal, categorical
    // Supported addenda: iid, series, circular

    /// <summary>
    /// A modify-able array of data.
    /// </summary>
    public sealed partial class FrameTable : FrameView
    {
        private FrameTable () : base(new List<NamedList>(), new List<int>())
        {

        }

        /// <summary>
        /// Initializes a new data frame with the columns specified by the given headers.
        /// </summary>
        /// <param name="columnHeaders"></param>
        public FrameTable(params ColumnDefinition[] columnHeaders) : this()
        {
            if (columnHeaders == null) throw new ArgumentNullException(nameof(columnHeaders));
            foreach(ColumnDefinition header in columnHeaders)
            {
                AddColumn(header.CreateList(this));
            }
        }

        internal FrameTable(params NamedList[] columns) : this((IEnumerable<NamedList>) columns)
        {

        }

        internal FrameTable(IEnumerable<NamedList> columns) : this()
        {
            if (columns == null) throw new ArgumentNullException(nameof(columns));
            foreach (NamedList column in columns)
            {
                AddColumn(column);
            }
        }


        /// <summary>
        /// Gets or sets the value in a given cell.
        /// </summary>
        /// <param name="rowIndex">The row index of the cell.</param>
        /// <param name="columnIndex">The column index of the cell.</param>
        public new object this [int rowIndex, int columnIndex] {
            get {
                return (base[rowIndex, columnIndex]);
            }
            set {
                columns[columnIndex].SetItem(map[rowIndex], value);
            }
        }

        /*
        /// <summary>
        /// Joins this data view to another data view on the given column.
        /// </summary>
        /// <param name="other"></param>
        /// <param name="columnName"></param>
        /// <returns></returns>
        public DataFrame Join (DataFrame other, string columnName)
        {
            int thisColumnIndex = this.GetColumnIndex(columnName);
            Type thisType = this.columns[thisColumnIndex].StorageType;
            int otherColumnIndex = other.GetColumnIndex(columnName);
            Type otherType = other.columns[otherColumnIndex].StorageType;

            // Form a lookup from the other table
            Dictionary<object, int> hash = new Dictionary<object, int>();
            for (int otherRowIndex = 0; otherRowIndex < other.Rows.Count; otherRowIndex++)
            {
                hash[other.columns[otherColumnIndex].GetItem(other.map[otherRowIndex])] = otherRowIndex;
            }

            // Construct the joined columns
            List<DataList> joinedColumns = new List<DataList>();
            for (int i = 0; i < this.columns.Count; i++)
            {
                DataList joinedColumn = DataList.Create(this.columns[i].Name, this.columns[i].StorageType);
                joinedColumns.Add(joinedColumn);
            }
            for (int j = 0; j < other.columns.Count; j++)
            {
                DataList joinedColumn = DataList.Create(other.columns[j].Name, other.columns[j].StorageType);
                joinedColumns.Add(joinedColumn);
            }

            // Populate the joined columns
            for (int thisRowIndex = 0; thisRowIndex < this.map.Count; thisRowIndex++)
            {
                object thisValue = this.columns[thisColumnIndex].GetItem(this.map[thisRowIndex]);
                int otherRowIndex;
                if (hash.TryGetValue(thisValue, out otherRowIndex))
                {
                    for (int i = 0; i < this.columns.Count; i++)
                    {
                        joinedColumns[i].AddItem(this.columns[i].GetItem(this.map[i]));
                    }
                    for (int j = 0; j < other.columns.Count; j++)
                    {
                        joinedColumns[this.columns.Count + j].AddItem(other.columns[j].GetItem(other.map[otherRowIndex]));
                    }
                }
            }

            DataFrame result = new DataFrame(joinedColumns);
            return (result);

        }
        */

        // IReadableDataList
        // IDataList : IReadableDataList
        // IReadableDataList<T> : IReadableDataList
        // IDataList<T> : IReadableDataList<T>, IDataList

        // DataView
        //    DataTable
        //    VirtualDataTable

        // DataColumn
        //    DataList
        //    VirtualDataColumn
        //    ComputedDataColumn

        /// <summary>
        /// Adds the given data list as a column to the data frame.
        /// </summary>
        /// <param name="column">The data list to add.</param>
        /// <returns>The index of the new column.</returns>
        internal int AddColumn(NamedList column)
        {
            if (column == null) throw new ArgumentNullException(nameof(column));
            if (this.columns.Count == 0) {
                for (int i = 0; i < column.Count; i++) {
                    this.map.Add(i);
                }
            } else {
                if (!column.IsComputed && column.Count != map.Count) throw new DimensionMismatchException();
            }
            int columnCount = this.columns.Count;
            this.columns.Add(column);
            this.columnMap[column.Name] = columnCount;
            return (columnCount);
        }

        /// <summary>
        /// Adds the given column to the data frame.
        /// </summary>
        /// <param name="column">The column to add.</param>
        public void AddColumn(ColumnDefinition column) {
            if (column == null) throw new ArgumentNullException(nameof(column));
            AddColumn(column.CreateList(this));
        }

        /// <summary>
        /// Removes the column with the given index from the data frame.
        /// </summary>
        /// <param name="columnName">The name of the column to remove.</param>
        public void RemoveColumn (string columnName)
        {
            if (columnName == null) throw new ArgumentNullException(nameof(columnName));
            int columnIndex = GetColumnIndex(columnName);
            this.columns.RemoveAt(columnIndex);
            // Remove the name of the removed column from the index,
            // and fix the index values for the higher columns, which will have changed 
            bool removed = this.columnMap.Remove(columnName);
            Debug.Assert(removed);
            for (int index = columnIndex; index < columns.Count; index++)
            {
                columnMap[columns[index].Name] = index;
            }
        }
       
        /// <summary>
        /// Add a new row of data to the data frame.
        /// </summary>
        /// <param name="values">The values to add to each column.</param>
        public void AddRow(params object[] values)
        {
            AddRow((IReadOnlyList<object>) values);
        }
        
        /// <summary>
        /// Adds a new row of data to the frame.
        /// </summary>
        /// <typeparam name="T">The type of the data collection.</typeparam>
        /// <param name="values">The values to add to each column.</param>
        public void AddRow<T>(IReadOnlyList<T> values) {
            if (values == null) throw new ArgumentNullException(nameof(values));
            if (values.Count != columns.Count) throw new DimensionMismatchException();
            int r = -1;
            for (int i = 0; i < values.Count; i++) {
                int previous_r = r;
                r = columns[i].AddItem(values[i]);
                if (previous_r > 0) Debug.Assert(r == previous_r);
            }
            map.Add(r);
        }

        /// <summary>
        /// Adds a new row of data to the data frame.
        /// </summary>
        /// <param name="values"></param>
        public void AddRow(IReadOnlyDictionary<string, object> values)
        {
            if (values == null) throw new ArgumentNullException(nameof(values));
            int rowCount = map.Count;
            for (int columnIndex = 0; columnIndex < columns.Count; columnIndex++)
            {
                NamedList column = columns[columnIndex];
                if (column.IsComputed) continue;
                object value = values[column.Name];
                int rowIndex = column.AddItem(value);
                if (rowIndex != rowCount) throw new InvalidOperationException();
            }
            map.Add(rowCount);
        }

        // RemoveRow and InsertRowAt are difficult because we need to not only re-order all the columns,
        // but also fix the indexes pointed to by all the map entries.
      
        /// <summary>
        /// Removes all data from the table.
        /// </summary>
        /// <remarks>
        /// <para>This method removes all rows, but does not remove any columns.</para>
        /// </remarks>
        public void Clear () {
            foreach(NamedList column in columns) {
                column.Clear();
            }
            map.Clear();
        }
        
        /// <summary>
        /// Creates a new data frame by reading values from a file of comma-seperated values.
        /// </summary>
        /// <param name="reader"></param>
        /// <returns>A new data frame with data from the file.</returns>
        /// <remarks>The column names are taken from the first line of the file.</remarks>
        public static FrameTable FromCsv (TextReader reader) {
            if (reader == null) throw new ArgumentNullException(nameof(reader));

            NamedList<string>[] textColumns;
            DataAdaptor[] headers;
            ReadCsvAsStrings(reader, out textColumns, out headers);

            NamedList[] columns = new NamedList[headers.Length];
            for (int columnIndex = 0; columnIndex < columns.Length; columnIndex++) {
                DataAdaptor header = headers[columnIndex];
                if (header.TypeCandidates.Count == 0) {
                    columns[columnIndex] = textColumns[columnIndex];
                } else {
                    TypeParser adaptor = header.TypeCandidates.First.Value;
                    NamedList column = adaptor.CreateStorage(textColumns[columnIndex].Name, header.IsNullable);
                    foreach (string textValue in textColumns[columnIndex]) {
                        if (textValue == null) {
                            column.AddItem(null);
                        } else {
                            object value = adaptor.Parse(textValue);
                            column.AddItem(value);
                        }
                    }
                    columns[columnIndex] = column;
                }
            }

            FrameTable frame = new FrameTable(columns);
            return (frame);
        }

        private static void ReadCsvAsStrings (TextReader reader, out NamedList<string>[] columns, out DataAdaptor[] headers) {
            Debug.Assert(reader != null);

            // Get the column names from the first line.

            string firstline = reader.ReadLine();
            if (firstline == null) {
                columns = null;
                headers = null;
                return;
            }

            List<string> names = CsvHelper.ReadCells(firstline);
            int count = names.Count;

            // Put the columns into lists of strings, and as we do so, maintain the collection of
            // types it can be parsed into, and whether any entries are null.

            columns = new NamedList<string>[names.Count];
            headers = new DataAdaptor[names.Count];
            for (int columnIndex = 0; columnIndex < columns.Length; columnIndex++) {
                columns[columnIndex] = new NamedList<string>(names[columnIndex]);
                headers[columnIndex] = new DataAdaptor();
            }

            while (true) {
                string line = reader.ReadLine();
                if (line == null) break;
                List<string> cells = CsvHelper.ReadCells(line);
                if (cells.Count != count) throw new FormatException();
                for (int columnIndex = 0; columnIndex < count; columnIndex++) {
                    string cell = cells[columnIndex];
                    if (String.IsNullOrEmpty(cell)) {
                        headers[columnIndex].IsNullable = true;
                        columns[columnIndex].Add(null);
                    } else {
                        DataAdaptor header = headers[columnIndex];
                        header.TryParse(cell);
                        columns[columnIndex].Add(cell);
                    }
                }
            }

        }

    }

}

