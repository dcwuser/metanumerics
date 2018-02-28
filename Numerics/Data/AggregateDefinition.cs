using System;
using System.Diagnostics;

namespace Meta.Numerics.Data {

    /// <summary>
    /// Represents a column defined by an aggregate function.
    /// </summary>
    public abstract class AggregateDefinition {

        internal AggregateDefinition (string name) {
            if (name == null) throw new ArgumentNullException(nameof(name));
            Name = name;
        }

        /// <summary>
        /// Gets the name of the aggregate column.
        /// </summary>
        public string Name { get; private set; }

        internal abstract Type StorageType { get; }

        internal abstract object ApplyAggregator (DataView view);

    }

    /// <summary>
    /// Represents a column defined by an aggregate function.
    /// </summary>
    /// <typeparam name="T">The type of the data.</typeparam>
    public sealed class AggregateDefinition<T> : AggregateDefinition {

        /// <summary>
        /// Creates a new instance of an aggregate column.
        /// </summary>
        /// <param name="name">The name of the aggregate column.</param>
        /// <param name="aggregator">The function that computes the aggregate value.</param>
        public AggregateDefinition(string name, Func<DataView, T> aggregator) : base(name) {
            if (aggregator == null) throw new ArgumentNullException(nameof(aggregator));
            Aggregator = aggregator;
        }

        /// <summary>
        /// Gets the function that computes the aggregate value.
        /// </summary>
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
