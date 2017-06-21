using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Data
{

    /// <summary>
    /// The header of a data column.
    /// </summary>
    /// <remarks>
    /// <para>The header contains basic data about a column, such as its name and the type
    /// of data it stores.</para>
    /// </remarks>
    public class DataHeader {
        /// <summary>
        /// Initializes a new data header with the given name and data type.
        /// </summary>
        /// <param name="name">The name of the column.</param>
        /// <param name="type">The type of data stored in the column.</param>
        public DataHeader(string name, Type type) {
            if (name == null) throw new ArgumentNullException(nameof(name));
            if (type == null) throw new ArgumentNullException(nameof(type));
            this.name = name;
            this.type = type;
        }

        private string name;

        private Type type;

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

        internal virtual DataList CreateList() {
            return (DataList.Create(name, type));
        }
    }

    /// <summary>
    /// The header of a data column of a specifc type.
    /// </summary>
    /// <typeparam name="T">The type of data stored in the column.</typeparam>
    public sealed class DataHeader<T> : DataHeader {

        /// <summary>
        /// Initializes a new data header for a column with the given name.
        /// </summary>
        /// <param name="name">The name of the column.</param>
        public DataHeader(string name) : base(name, typeof(T)) {
        }

        internal override DataList CreateList () {
            return (new DataList<T>(this.Name));
        }

    }
}
