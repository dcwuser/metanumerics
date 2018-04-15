using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

namespace Meta.Numerics.Data
{

    // Consider adding flags to indicate data types and properties, and having statistics
    // functions check for these flags.
    // Data types could include: continuous, count, ordinal, categorical
    // Additional flags could include: iid, series, circular

    /// <summary>
    /// A modify-able array of data.
    /// </summary>
    /// <remarks>
    /// <para>This is the central class for storing data in our data frame system.</para>
    /// <para>Use the <see cref="FromCsv"/> method to create a frame table from a comma-separated values
    /// file or the <see cref="FromDictionaries"/> method to create a frame table from JSON or another
    /// collection-of-dictionaries representation. Or create one programmatically by using the
    /// <see cref="FrameTable()"/> constructor to instantiate an empty data frame and
    /// <see cref="AddColumn{T}(string)"/> and <see cref="AddRow(object[])"/> to add columns and rows.</para>
    /// <para>Using the methods inherited from the <see cref="FrameView"/> class, you can filter, re-order,
    /// manipulate, and analyze data without incurring the space or time costs of copying the stored data.</para>
    /// </remarks>
    public sealed partial class FrameTable : FrameView {

        /// <summary>
        /// Initializes a new, empty frame table.
        /// </summary>
        public FrameTable () : base(new List<NamedList>(), new List<int>()) {

        }

        internal FrameTable(params NamedList[] columns) : this((IEnumerable<NamedList>) columns) {

        }

        internal FrameTable(IEnumerable<NamedList> columns) : this() {
            if (columns == null) throw new ArgumentNullException(nameof(columns));
            foreach (NamedList column in columns) {
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

        internal void AddColumn(NamedList column)
        {
            if (column == null) throw new ArgumentNullException(nameof(column));
            if (this.columns.Count == 0) {
                // This is the first column; create a row map.
                for (int i = 0; i < column.Count; i++) {
                    this.map.Add(i);
                }
            } else {
                // This is not the first column; if it isn't computed, it's length must match the existing columns.
                if (!column.IsComputed && column.Count != map.Count) throw new DimensionMismatchException();
            }
            this.columnMap.Add(column.Name, this.columns.Count);
            this.columns.Add(column);
        }

        /// <summary>
        /// Adds a new column with the given name and type.
        /// </summary>
        /// <typeparam name="T">The type of the column.</typeparam>
        /// <param name="name">The name of the column.</param>
        public void AddColumn<T>(string name) {
            if (name == null) throw new ArgumentNullException(nameof(name));
            AddColumn(new NamedList<T>(name));
        }

        /// <summary>
        /// Adds a new column with the given name and stored values.
        /// </summary>
        /// <typeparam name="T">The type of the column.</typeparam>
        /// <param name="name">The name of the column.</param>
        /// <param name="storage">The stored values.</param>
        /// <remarks>
        /// <para>Use this overload when you already have a collection of column values
        /// and want to use it as a column in an existing <see cref="FrameTable"/>.
        /// Note that the values are not copied; the collection you provide is used directly
        /// as the column storage. Therefore any changes you make to the collection after
        /// adding the column will affect the frame table. To avoid problems, do not
        /// alter the collection after submitting it to this method.</para>
        /// </remarks>
        public void AddColumn<T>(string name, List<T> storage) {
            if (name == null) throw new ArgumentNullException(nameof(name));
            if (storage == null) throw new ArgumentNullException(nameof(storage));
            AddColumn(new NamedList<T>(name, storage));
        }

        /// <summary>
        /// Adds the new columns with the given names.
        /// </summary>
        /// <typeparam name="T">The type of the columns.</typeparam>
        /// <param name="names">The names of the columns.</param>
        /// <remarks>
        /// <para>This method is useful for adding multiple columns of the same type.</para>
        /// </remarks>
        public void AddColumns<T>(params string[] names) {
            if (names == null) throw new ArgumentNullException(nameof(names));
            foreach (string name in names) {
                AddColumn(new NamedList<T>(name));
            }
        }

        /// <summary>
        /// Removes the column with the given index from the data frame.
        /// </summary>
        /// <param name="columnName">The name of the column to remove.</param>
        public void RemoveColumn (string columnName) {
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
        public void AddRow(params object[] values) {
            AddRow((IReadOnlyList<object>) values);
        }
        
        /// <summary>
        /// Adds a new row of data to the frame.
        /// </summary>
        /// <typeparam name="T">The type of the data collection.</typeparam>
        /// <param name="values">The values to add to each column.</param>
        /// <remarks>
        /// <para>Use this overload when the data types are uniform across all columns.</para>
        /// </remarks>
        /// <exception cref="DimensionMismatchException">The number of items in <paramref name="values"/>
        /// not equal the number of columns in the table.</exception>
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
        /// <param name="values">A dictionary that maps the existing column names to the cell values for the new row.</param>
        public void AddRow(IReadOnlyDictionary<string, object> values) {
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
        /// Creates a new frame table from a file of comma-separated values.
        /// </summary>
        /// <param name="reader">A reader positioned at the beginning of the file.</param>
        /// <returns>A new data frame with data from the file.</returns>
        /// <remarks>
        /// <para>The column names are taken from the first line of the file.</para>
        /// <para>The storage type of each column is inferred from the types of objects
        /// encountered are the frame table is constructed.</para>
        /// </remarks>
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

        /// <summary>
        /// Constructs a new frame table from a sequence of dictionaries.
        /// </summary>
        /// <param name="dictionaries">A enumerable set of dictionaries, one for each row, whose
        /// keys are the column names and whose values are the cell value of the column for that row.</param>
        /// <returns>A new frame table with data from the dictionaries.</returns>
        /// <remarks><para>The storage type of each column is inferred from the types of objects
        /// encountered are the frame table is constructed.</para></remarks>
        public static FrameTable FromDictionaries (IEnumerable<IReadOnlyDictionary<string, object>> dictionaries) {
            if (dictionaries == null) throw new ArgumentNullException(nameof(dictionaries));

            List<NamedList> columns = DictionaryHelper.ReadDictionaries(dictionaries);

            // Collect the results into a frame table.
            FrameTable result = new FrameTable(columns);
            return (result);
        }

    }

}

