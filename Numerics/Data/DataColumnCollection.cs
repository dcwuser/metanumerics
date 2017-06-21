using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Data
{

    /// <summary>
    /// A collection of data columns.
    /// </summary>
    public sealed class DataColumnCollection : IReadOnlyCollection<DataColumn<object>>, IReadOnlyList<DataColumn<object>>
    {
        internal DataColumnCollection(DataView frame)
        {
            Debug.Assert(frame != null);
            this.frame = frame;
        }

        private readonly DataView frame;

        /// <summary>
        /// Gets the number of columns.
        /// </summary>
        public int Count
        {
            get
            {
                return(frame.columns.Count);
            }
        }

        /// <summary>
        /// Gets the column with the given index.
        /// </summary>
        /// <param name="index">The index of the column.</param>
        /// <returns>The column with index <paramref name="index"/>.</returns>
        public DataColumn<object> this [int index]
        {
            get
            {
                if ((index < 0) || (index >= frame.columns.Count)) throw new ArgumentOutOfRangeException(nameof(index));
                return (new DataColumn<object>(frame, index));
            }
        }

        IEnumerator<DataColumn<object>> IEnumerable<DataColumn<object>>.GetEnumerator()
        {
            throw new NotImplementedException();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            throw new NotImplementedException();
        }
    }
}
