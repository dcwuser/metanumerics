using System;
using System.Collections;
using System.Collections.Generic;

namespace Meta.Numerics.Data
{

    /// <summary>
    /// A readable, indexed set of data points. 
    /// </summary>
    /// <typeparam name="T">The type of the data points.</typeparam>
    public class DataColumn<T> : IReadOnlyDataList<T>
    {
        internal DataColumn(DataView frame, int c)
        {
            this.frame = frame;
            this.c = c;
        }

        private DataView frame;
        private int c;

        /// <summary>
        /// Gets the name of the column.
        /// </summary>
        public string Name
        {
            get
            {
                return (frame.columns[c].Name);
            }
            /*
            set
            {
                frame.columns[c].Name = value;
            }
            */
        }

        /// <summary>
        /// Gets the type of data stored in the column.
        /// </summary>
        public Type StorageType
        {
            get
            {
                return (frame.columns[c].StorageType);
            }
        }

        /// <summary>
        /// Gets the number of values in the column.
        /// </summary>
        public int Count
        {
            get
            {
                return (frame.map.Count);
            }
        }

        /// <summary>
        /// Gets the value at the given index.
        /// </summary>
        /// <param name="r">The (zero-based) index.</param>
        /// <returns>The value at index <paramref name="r"/>.</returns>
        public T this[int r]
        {
            get
            {
                return ((T)frame.columns[c].GetItem(frame.map[r]));
            }
        }

        /*
        object IDataList.this[int r]
        {
            get
            {
                return (this[r]);
            }
        }
        */

        IEnumerator<T> IEnumerable<T>.GetEnumerator()
        {
            return ((IEnumerable<T>)frame.columns[c]).GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return (((IEnumerable<T>)this).GetEnumerator());
        }

        /*
        int IDataList.Add(object value)
        {
            throw new NotImplementedException();
        }
        */

    }
}
