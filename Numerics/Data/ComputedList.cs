using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Data
{

    internal class ComputedList<T> : NamedList, IEnumerable, IReadOnlyList<T>
    {
        internal ComputedList (FrameView frame, string name, Func<FrameRow,T> function) : base(name) {
            this.frame = frame;
            this.name = name;
            this.function = function;
        }

        private readonly FrameView frame;
        private readonly string name;
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
