using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace Meta.Numerics.Data
{
    [DebuggerDisplay("{DebuggerDisplay()}")]
    public sealed class DataRow : IReadOnlyDictionary<string, object>, IReadOnlyList<object> {

        internal DataRow (DataView frame, int r) {
            Debug.Assert(frame != null);
            Debug.Assert(r >= 0);
            this.frame = frame;
            this.r = r;
        }

        private readonly DataView frame;
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
                foreach (DataList column in frame.columns) {
                    yield return (column.Name);
                }
            }
        }

        IEnumerable<object> IReadOnlyDictionary<string, object>.Values {
            get {
                foreach (DataList column in frame.columns) {
                    yield return (column.GetItem(frame.map[r]));
                }
            }
        }

        /// <summary>
        /// Gets the value of the column with the given index.
        /// </summary>
        /// <param name="c">The index of the column.</param>
        /// <returns>The value of the <paramref name="c"/>th column.</returns>
        public object this[int c] {
            get {
                if ((c < 0) || (c >= Count)) throw new ArgumentOutOfRangeException(nameof(c));
                return (frame.columns[c].GetItem(frame.map[r]));
            }
        }

        /// <summary>
        /// Gets the value of the column with the given name.
        /// </summary>
        /// <param name="name">The name of the column.</param>
        /// <returns>The value.</returns>
        public object this[string name] {
            get {
                int c = frame.GetColumnIndex(name);
                return (this[c]);
            }
        }

        private string DebuggerDisplay () {
            StringBuilder builder = new StringBuilder();
            builder.Append("{ ");
            for (int c = 0; c < Math.Min(4, frame.columns.Count - 1); c++) {
                builder.AppendFormat("{0} ", this[c]);
            }
            if (frame.columns.Count > 4) builder.Append("... ");
            if (frame.columns.Count > 0) builder.Append(this[frame.columns.Count - 1]);
            builder.Append(" }");
            return (builder.ToString());
        }

        bool IReadOnlyDictionary<string, object>.ContainsKey (string key) {
            throw new NotImplementedException();
        }

        bool IReadOnlyDictionary<string, object>.TryGetValue (string key, out object value) {
            throw new NotImplementedException();
        }

        IEnumerator<KeyValuePair<string, object>> IEnumerable<KeyValuePair<string, object>>.GetEnumerator () {
            foreach (DataList column in frame.columns) {
                yield return (new KeyValuePair<string, object>(column.Name, column.GetItem(frame.map[r])));
            }
        }

        IEnumerator<object> IEnumerable<object>.GetEnumerator () {
            foreach (DataList column in frame.columns) {
                yield return (column.GetItem(frame.map[r]));
            }
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (((IEnumerable<KeyValuePair<string, object>>) this).GetEnumerator());
        }

    }

}
