using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Data
{

    /// <summary>
    /// Represents a collection of columns of a data frame.
    /// </summary>
    public sealed class FrameColumnCollection : IReadOnlyCollection<FrameColumn>, IReadOnlyList<FrameColumn>, IEnumerable<FrameColumn> /*, IReadOnlyDictionary<string,FrameColumn> */
    {
        internal FrameColumnCollection(FrameView frame)
        {
            Debug.Assert(frame != null);
            this.frame = frame;
        }

        private readonly FrameView frame;

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
        public FrameColumn this [int index] {
            get {
                if ((index < 0) || (index >= frame.columns.Count)) throw new ArgumentOutOfRangeException(nameof(index));
                return (new FrameColumn(frame, index));
            }
        }

        /// <summary>
        /// Gets the column with the given name.
        /// </summary>
        /// <param name="name">The name of the column.</param>
        /// <returns>The column with the given name.</returns>
        public FrameColumn this [string name] {
            get {
                int index = frame.GetColumnIndex(name);
                return (new FrameColumn(frame, index));
            }
        }

        /// <summary>
        /// Gets the index of the column with the given name.
        /// </summary>
        /// <param name="name">The name of the column.</param>
        /// <returns>The index of the column, or -1 if no such column exists.</returns>
        public int GetIndexOf (string name) {
            return (frame.GetColumnIndex(name));
        }

        /// <summary>
        /// Gets an enumerator that iterates through the columns.
        /// </summary>
        /// <returns>The requested enumerator.</returns>
        IEnumerator<FrameColumn> IEnumerable<FrameColumn>.GetEnumerator () {
            for (int i = 0; i < frame.columns.Count; i++) {
                yield return (new FrameColumn(frame, i));
            }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return (((IEnumerable<FrameColumn>) this).GetEnumerator());
        }

    }
}
