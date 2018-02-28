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
    public class ColumnDefinition {

        /// <summary>
        /// Initializes a new data header with the given name and data type.
        /// </summary>
        /// <param name="name">The name of the column.</param>
        /// <param name="type">The type of data stored in the column.</param>
        public ColumnDefinition (string name, Type type) {
            if (name == null) throw new ArgumentNullException(nameof(name));
            if (type == null) throw new ArgumentNullException(nameof(type));
            this.name = name;
            this.type = type;
        }

        private readonly string name;

        private readonly Type type;

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
        public Type StorageType {
            get {
                return (type);
            }
        }

        internal virtual DataList CreateList (DataView view) {
            return (DataList.Create(name, type));
        }
    }

    /// <summary>
    /// Represents a column of a data frame that stores the given type of data.
    /// </summary>
    /// <typeparam name="T">The type of data stored in the column.</typeparam>
    public sealed class ColumnDefinition<T> : ColumnDefinition {

        /// <summary>
        /// Initializes a new data header for a column with the given name.
        /// </summary>
        /// <param name="name">The name of the column.</param>
        public ColumnDefinition (string name) : base(name, typeof(T)) {
        }

        internal override DataList CreateList (DataView view) {
            return (new DataList<T>(this.Name));
        }

    }

    /// <summary>
    /// Represents a column of a data frame that is computed from other column values.
    /// </summary>
    /// <typeparam name="T">The type of data computed.</typeparam>
    public sealed class ComputedColumnDefinition<T> : ColumnDefinition {

        private readonly Func<DataRow, T> function;

        public ComputedColumnDefinition (string name, Func<DataRow, T> function) : base(name, typeof(T)) {
            if (function == null) throw new ArgumentNullException(nameof(function));
            this.function = function;
        }

        internal override DataList CreateList (DataView view) {
            Debug.Assert(view != null);
            return (new ComputedColumn<T>(view, this.Name, function));
        }
    }
}
