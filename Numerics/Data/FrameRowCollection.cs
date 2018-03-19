using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Data
{

    /// <summary>
    /// A collection of rows from a frame.
    /// </summary>
    public sealed class FrameRowCollection : IReadOnlyCollection<FrameRow>, IReadOnlyList<FrameRow>, IEnumerable<FrameRow> {

        internal FrameRowCollection (FrameView frame) {
            Debug.Assert(frame != null);
            this.frame = frame;
        }

        private readonly FrameView frame;

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
        public FrameRow this[int index] {
            get {
                return (new FrameRow(frame, index));
            }
        }

        /// <summary>
        /// Returns an enumerator that iterates through the rows.
        /// </summary>
        /// <returns>The requested enumerator.</returns>
        IEnumerator<FrameRow> IEnumerable<FrameRow>.GetEnumerator () {
            for (int i = 0; i < frame.map.Count; i++) {
                yield return (new FrameRow(frame, i));
            }
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (((IEnumerable<FrameRow>) this).GetEnumerator());
        }

    }

}
