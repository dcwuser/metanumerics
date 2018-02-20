using System;
using System.Diagnostics;

namespace Meta.Numerics.Data {

    public abstract class AggregateColumn {

        internal AggregateColumn (string name) {
            if (name == null) throw new ArgumentNullException(nameof(name));
            Name = name;
        }

        public string Name { get; private set; }

        internal abstract Type StorageType { get; }

        internal abstract object ApplyAggregator (DataView view);

    }

    public sealed class AggregateColumn<T> : AggregateColumn {

        public AggregateColumn(string name, Func<DataView, T> aggregator) : base(name) {
            if (aggregator == null) throw new ArgumentNullException(nameof(aggregator));
            Aggregator = aggregator;
        }

        public Func<DataView, T> Aggregator { get; private set; }

        internal override Type StorageType {
            get {
                return (typeof(T));
            }
        }

        internal override object ApplyAggregator (DataView view) {
            Debug.Assert(view != null);
            Debug.Assert(Aggregator != null);
            return (Aggregator(view));
        }

    }
}
