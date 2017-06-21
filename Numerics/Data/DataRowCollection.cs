using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Data
{

    /// <summary>
    /// A collection of rows from a table.
    /// </summary>
    public sealed class DataRowCollection : IReadOnlyCollection<DataRow>, IReadOnlyList<DataRow> {
        internal DataRowCollection (DataView frame) {
            Debug.Assert(frame != null);
            this.frame = frame;
        }

        private readonly DataView frame;

        /// <summary>
        /// Gets the number of rows.
        /// </summary>
        public int Count {
            get {
                return (frame.map.Count);
            }
        }

        /// <summary>
        /// Gets the row with the given index.
        /// </summary>
        /// <param name="index">The (zero-based) index.</param>
        /// <returns></returns>
        public DataRow this[int index] {
            get {
                return (new DataRow(frame, index));
            }
        }

        public IEnumerator<DataRow> GetEnumerator () {
            for (int i = 0; i < frame.map.Count; i++) {
                yield return (new DataRow(frame, i));
            }
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (((IEnumerable<DataRow>) this).GetEnumerator());
        }

    }

}
