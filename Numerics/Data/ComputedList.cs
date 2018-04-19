using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Data
{

    internal class ComputedList<T> : NamedList, IEnumerable, IEnumerable<T>, IReadOnlyList<T>
    {
        internal ComputedList (FrameView frame, string name, Func<FrameRow,T> function) : base(name) {
            Debug.Assert(frame != null);
            Debug.Assert(function != null);
            this.frame = frame;
            this.function = function;
        }

        private readonly FrameView frame;
        private readonly Func<FrameRow, T> function;

        public T this[int index]
        {
            get
            {
                FrameRow row = new FrameRow(this.frame, index);
                return (function(row));
            }
        }

        public override int Count
        {
            get
            {
                return (frame.map.Count);
            }
        }

        public override Type StorageType
        {
            get
            {
                return (typeof(T));
            }
        }

        internal override bool IsComputed {
            get {
                return (true);
            }
        }

        public IEnumerator<T> GetEnumerator()
        {
            for (int r = 0; r < frame.map.Count; r++)
            {
                yield return function(new FrameRow(frame, r));
            }
        }
        
        IEnumerator IEnumerable.GetEnumerator() {
            return (((IEnumerable<T>) this).GetEnumerator());
        }

        internal override object GetItem (int index)
        {
            return (this[index]);
        }

        internal override void SetItem(int index, object value)
        {
            throw new InvalidOperationException();
        }

        internal override int AddItem(object value)
        {
            throw new InvalidOperationException();
        }

        public override void Clear () {
            // no-op
        }

    }
}
