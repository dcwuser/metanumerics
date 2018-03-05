using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Data {

    /// <summary>
    /// Defines a column of a data frame.
    /// </summary>
    /// <remarks><para>This is an abstract marker class. To define a concrete column, use <see cref="ColumnDefinition{T}"/>.</para></remarks>
    public abstract class ColumnDefinition {

        internal ColumnDefinition(string name) {
            Debug.Assert(name != null);
            this.name = name;
        }

        private readonly string name;

        /// <summary>
        /// Gets the name of the column.
        /// </summary>
        public string Name {
            get {
                return (name);
            }
        }

        /// <summary>
        /// Gets the type of data stored in the column.
        /// </summary>
        public abstract Type StorageType { get; }

        internal abstract NamedList CreateList(FrameView view);

    }

    /// <summary>
    /// Defines a column of a data frame that stores the given type of data.
    /// </summary>
    /// <typeparam name="T">The type of data stored in the column.</typeparam>
    /// <remarks>
    /// <para>Pass instances of this type to the constructor <see cref="FrameTable.FrameTable(ColumnDefinition[])"/> to instantiate a new <see cref="FrameTable"/>.</para>
    /// <para>If you want a column to support null value types, use <see cref="Nullable{T}"/> types for T.</para>
    /// </remarks>
    public sealed class ColumnDefinition<T> : ColumnDefinition {

        /// <summary>
        /// Initializes a new definition of a column with the given name.
        /// </summary>
        /// <param name="name">The name of the column.</param>
        public ColumnDefinition (string name) : this(name, null) {
        }

        /// <summary>
        /// Initializes a new definition of a column with the given name and storage.
        /// </summary>
        /// <param name="name">The name of the column.</param>
        /// <param name="storage">The storage to be used.</param>
        public ColumnDefinition (string name, List<T> storage) : base(name) {
            this.storage = storage;
        }


        private readonly List<T> storage;

        /// <inheritdoc/>
        public override Type StorageType {
            get {
                return (typeof(T));
            }
        }

        internal override NamedList CreateList (FrameView view) {
            if (storage == null) {
                return (new NamedList<T>(this.Name));
            } else {
                return (new NamedList<T>(this.Name, storage));
            }
        }

    }

    /// <summary>
    /// Represents a column of a data frame that is computed from other column values.
    /// </summary>
    /// <typeparam name="T">The type of data computed.</typeparam>
    public sealed class ComputedColumnDefinition<T> : ColumnDefinition {

        private readonly Func<FrameRow, T> function;

        /// <summary>
        /// Initializes a new definition of a computed column.
        /// </summary>
        /// <param name="name"></param>
        /// <param name="function"></param>
        public ComputedColumnDefinition (string name, Func<FrameRow, T> function) : base(name) {
            if (function == null) throw new ArgumentNullException(nameof(function));
            this.function = function;
        }

        /// <inheritdoc/>
        public override Type StorageType {
            get {
                return (typeof(T));
            }
        }

        internal override NamedList CreateList (FrameView view) {
            Debug.Assert(view != null);
            return (new ComputedDataList<T>(view, this.Name, function));
        }

    }
}
