using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Data
{

    /// <summary>
    /// A collection of data columns.
    /// </summary>
    public sealed class DataColumnCollection : IReadOnlyCollection<DataColumn<object>>, IReadOnlyList<DataColumn<object>>, IReadOnlyDictionary<string,DataColumn<object>>
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
        public int Count {
            get {
                return(frame.columns.Count);
            }
        }

        /// <summary>
        /// Gets the column with the given index.
        /// </summary>
        /// <param name="index">The index of the column.</param>
        /// <returns>The column with index <paramref name="index"/>.</returns>
        public DataColumn<object> this [int index] {
            get {
                if ((index < 0) || (index >= frame.columns.Count)) throw new ArgumentOutOfRangeException(nameof(index));
                return (new DataColumn<object>(frame, index));
            }
        }

        /// <summary>
        /// Gets the column with the given name.
        /// </summary>
        /// <param name="name">The name of the column.</param>
        /// <returns>The column with the given name.</returns>
        public DataColumn<object> this [string name] {
            get {
                int index = frame.GetColumnIndex(name);
                return (new DataColumn<object>(frame, index));
            }
        }

        IEnumerator<DataColumn<object>> IEnumerable<DataColumn<object>>.GetEnumerator()
        {
            for (int i = 0; i < frame.columns.Count; i++) {
                yield return (new DataColumn<object>(frame, i));
            }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return (((IEnumerable<DataColumn<object>>) this).GetEnumerator());
        }

        bool IReadOnlyDictionary<string, DataColumn<object>>.ContainsKey (string key) {
            return (frame.GetColumnIndex(key) >= 0);
        }

        bool IReadOnlyDictionary<string, DataColumn<object>>.TryGetValue (string key, out DataColumn<object> value) {
            throw new NotImplementedException();
        }

        IEnumerator<KeyValuePair<string, DataColumn<object>>> IEnumerable<KeyValuePair<string, DataColumn<object>>>.GetEnumerator () {
            foreach(DataColumn<object> column in ((IEnumerable<DataColumn<object>>) this)) {
                yield return (new KeyValuePair<string, DataColumn<object>>(column.Name, column));
            }
        }

        IEnumerable<string> IReadOnlyDictionary<string, DataColumn<object>>.Keys {
            get {
                return (frame.columnMap.Keys);
            }
        }

        IEnumerable<DataColumn<object>> IReadOnlyDictionary<string, DataColumn<object>>.Values {
            get {
                throw new NotImplementedException();
            }
        }

    }
}
