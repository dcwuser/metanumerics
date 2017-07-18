using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Reflection;

namespace Meta.Numerics.Data
{

    // Supported storage types: int, double, bool, datetime, string
    // Supported data types: continuous, count, ordinal, categorical
    // Supported addenda: iid, series, circular

    /*
    public abstract class DataView<R, C>
    {
        public abstract DataView<R, C> Select(params C[] columns);

        public abstract DataView<R, C> Where<T>(C column, Func<T, bool> selector);

        public abstract DataView<R, C> Where(Func<DataRow, bool> selector);

        public abstract DataView<R, C> OrderBy(C column);

        public abstract DataView<R, C> OrderBy<T>(C column, Comparer<T> comparer);

        public abstract DataView<R, C> OrderBy(Comparer<DataRow> comparer);

        public abstract DataView<int, C> GroupBy<T>(C groupColumn, Func<DataView<R, C>, T> aggregator, C aggregateColumn);

        public abstract DataView<R, C> DiscardNulls(params C[] columns);

        public abstract object CrossTabs(C rowColumn, C columnColumn);

        public abstract DataList<T> Column<T>(C column);

        public abstract DataRow Row(R row);

    }
    */

    /*
    public class DataViewImplementation<R, C> : DataView<R, C>
    {
        private List<IDataList> columns;
        private Dictionary<C, int> columnIndex;
        private Dictionary<R, int> rowIndex;

        private int GetRowIndex (R row)
        {
            int index = -1;
            rowIndex.TryGetValue(row, out index);
            return (index);
        }

        public override DataFrameRow Row(R row)
        {
            int index = GetRowIndex(row);
            return (null);
            //return (new DataFrameRow(this, index));
        }

    }
    */

    /// <summary>
    /// A modifyable array of data.
    /// </summary>
    public sealed partial class DataFrame : DataView
    {
        private DataFrame ()
        {
            this.columns = new List<DataList>();
            this.columnMap = new Dictionary<string, int>();
            this.map = new List<int>();
        }

        /// <summary>
        /// Initializes a new data frame with the columns specifed by the given headers.
        /// </summary>
        /// <param name="columnHeaders"></param>
        public DataFrame(params DataHeader[] columnHeaders) : this()
        {
            if (columnHeaders == null) throw new ArgumentNullException(nameof(columnHeaders));
            foreach(DataHeader header in columnHeaders)
            {
                AddColumn(header.CreateList());
            }
        }

        /// <summary>
        /// Initializes a new data frame with the given data lists as columns.
        /// </summary>
        /// <param name="columns">The data lists.</param>
        public DataFrame(params DataList[] columns) : this((IList<DataList>) columns)
        {

        }

        /// <summary>
        /// Initializes a new data frame with the given sequence of data lists as columns.
        /// </summary>
        /// <param name="columns">An enumerable sequence of non-null data lists.</param>
        public DataFrame(IEnumerable<DataList> columns) : this()
        {
            if (columns == null) throw new ArgumentNullException(nameof(columns));
            foreach (DataList column in columns)
            {
                AddColumn(column);
            }
            int count = this.columns[0].Count;
            for (int i = 0; i < count; i++)
            {
                this.map.Add(i);
            }
        }

        internal DataFrame(List<DataList> columns, List<int> map) :base(columns, map)
        {

        }

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
        public int AddColumn(DataList column)
        {
            if (column == null) throw new ArgumentNullException(nameof(column));
            int columnCount = columns.Count;
            if ((columnCount > 0) && (columns[0].Count != column.Count)) throw new InvalidOperationException();
            this.columns.Add(column);
            this.columnMap[column.Name] = columnCount;
            return (columnCount);
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
            this.columnMap.Remove(columnName);
            for (int index = columnIndex; index < columns.Count; index++)
            {
                columnMap[columns[index].Name] = index;
            }
        }
       
        /// <summary>
        /// Add a new row of data to the data frame.
        /// </summary>
        /// <param name="values"></param>
        public void AddRow(params object[] values)
        {
            if (values == null) throw new ArgumentNullException(nameof(values));
            if (values.Length != columns.Count) throw new InvalidOperationException();
            int r = 0;
            for (int i = 0; i < values.Length; i++)
            {
                r = columns[i].AddItem(values[i]);
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
                DataList column = columns[columnIndex];
                // handle computed columns
                object value = values[column.Name];
                int rowIndex = column.AddItem(value);
                if (rowIndex != rowCount) throw new InvalidOperationException();
            }
            map.Add(rowCount);
        }


        // parse to integer -> int
        // parse to double -> double
        // parse to datetime -> datetime
        // parse to boolean -> boolean
        // parse to guid -> guid
        // leave as string        

    }

}
