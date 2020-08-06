using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace Meta.Numerics.Data
{

    /// <summary>
    /// Represents one row of a data frame.
    /// </summary>
    public sealed class FrameRow : IReadOnlyDictionary<string, object>, IReadOnlyList<object> {

        internal FrameRow (FrameView frame, int r) {
            Debug.Assert(frame != null);
            Debug.Assert(r >= 0);
            this.frame = frame;
            this.r = r;
        }

        private readonly FrameView frame;
        private readonly int r;

        /// <summary>
        /// Gets the number of columns in the row.
        /// </summary>
        public int Count {
            get {
                return (frame.columns.Count);
            }
        }

        IEnumerable<string> IReadOnlyDictionary<string, object>.Keys {
            get {
                foreach (NamedList column in frame.columns) {
                    yield return (column.Name);
                }
            }
        }

        IEnumerable<object> IReadOnlyDictionary<string, object>.Values {
            get {
                foreach (NamedList column in frame.columns) {
                    yield return (column.GetItem(frame.map[r]));
                }
            }
        }

        /// <summary>
        /// Gets the value of the column with the given index.
        /// </summary>
        /// <param name="c">The index of the column.</param>
        /// <returns>The value of the <paramref name="c"/>th column.</returns>
        /// <exception cref="IndexOutOfRangeException"><paramref name="c"/> is less than zero or greater than or equal to the number of columns <see cref="Count"/>.</exception>
        public object this[int c] {
            get {
                if ((c < 0) || (c >= Count)) throw new IndexOutOfRangeException();
                return (frame.columns[c].GetItem(frame.map[r]));
            }
        }

        /// <summary>
        /// Gets the value of the column with the given name.
        /// </summary>
        /// <param name="name">The name of the column.</param>
        /// <returns>The value.</returns>
        /// <exception cref="IndexOutOfRangeException">There is no column in the row with the given <paramref name="name"/>.</exception>
        public object this[string name] {
            get {
                int c = frame.GetColumnIndex(name);
                return (this[c]);
            }
        }


        /// <summary>
        /// Gets the index of the column with the given name.
        /// </summary>
        /// <param name="name">The name of the column.</param>
        /// <returns>The index of the column, of -1 if no such column exists.</returns>
        public int GetIndexOf (string name) {
            return (frame.GetColumnIndex(name));
        }

        /// <summary>
        /// Returns a dictionary containing the row data.
        /// </summary>
        /// <returns>A dictionary whose keys are the column names
        /// and whose values are the corresponding cell values for the row.
        /// This dictionary is an independent copy, so any subsequent changes to the
        /// frame data will not change it.</returns>
        public Dictionary<string, object> ToDictionary () {
            Dictionary<string, object> result = new Dictionary<string, object>();
            foreach (KeyValuePair<string, object> entry in (IReadOnlyDictionary<string, object>) this) {
                result.Add(entry.Key, entry.Value);
            }
            return (result);
        }

        bool IReadOnlyDictionary<string, object>.ContainsKey (string key) {
            int c = frame.GetColumnIndex(key);
            return (c >= 0);
        }

        bool IReadOnlyDictionary<string, object>.TryGetValue (string key, out object value) {
            int c = frame.GetColumnIndex(key);
            if (c < 0) {
                value = null;
                return (false);
            } else {
                value = this[c];
                return (true);
            }
        }

        IEnumerator<KeyValuePair<string, object>> IEnumerable<KeyValuePair<string, object>>.GetEnumerator () {
            foreach (NamedList column in frame.columns) {
                yield return (new KeyValuePair<string, object>(column.Name, column.GetItem(frame.map[r])));
            }
        }

        IEnumerator<object> IEnumerable<object>.GetEnumerator () {
            foreach (NamedList column in frame.columns) {
                yield return (column.GetItem(frame.map[r]));
            }
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (((IEnumerable<KeyValuePair<string, object>>) this).GetEnumerator());
        }

    }

}
