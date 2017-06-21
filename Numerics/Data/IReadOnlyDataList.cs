using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Data
{
    /// <summary>
    /// Represents the interface of a readable data list. 
    /// </summary>
    /// <typeparam name="T">The type of the data stored in the list.</typeparam>
    /// <remarks>
    /// <para>This interface is a simple extension of the <see cref="IReadOnlyList{T}"/> interface, which carries
    /// along a name for the data list.</para>
    /// </remarks>
    public interface IReadOnlyDataList<out T> : IReadOnlyList<T> {

        /// <summary>
        /// Gets the name of the data list.
        /// </summary>
        string Name { get; }

        /// <summary>
        /// Gets the type of data stored in the list.
        /// </summary>
        Type StorageType { get; }

    }

}
