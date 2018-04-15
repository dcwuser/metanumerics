using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Reflection;

namespace Meta.Numerics.Data
{

    /// <summary>
    /// Represents one column of a data frame.
    /// </summary>
    /// <remarks>
    /// <para>To obtained strongly typed values of the data in the column,
    /// use the <see cref="As{T}" autoUpgrade="true"/> method.</para>
    /// </remarks>
    public sealed class FrameColumn : IEnumerable {

        internal FrameColumn(FrameView frame, int c) {
            this.frame = frame;
            this.column = frame.columns[c];
            this.map = frame.map;
        }

        private readonly FrameView frame;
        private readonly NamedList column;
        private readonly List<int> map;

        /// <summary>
        /// Gets the name of the column.
        /// </summary>
        public string Name {
            get {
                return (column.Name);
            }
        }

        /// <summary>
        /// Gets the type of data stored in the column.
        /// </summary>
        public Type StorageType {
            get {
                return (column.StorageType);
            }
        }

        /// <summary>
        /// Gets the number of values in the column.
        /// </summary>
        public int Count {
            get {
                return (map.Count);
            }
        }

        /// <summary>
        /// Gets the value at the given index.
        /// </summary>
        /// <param name="r">The (zero-based) index.</param>
        /// <returns>The value at index <paramref name="r"/>.</returns>
        public object this[int r] {
            get {
                return (column.GetItem(map[r]));
            }
        }

        /// <summary>
        /// Exposes the column as a list of uniform type.
        /// </summary>
        /// <typeparam name="T">The type of the list.</typeparam>
        /// <returns>An object that allows the column to be accessed as a typed list.</returns>
        /// <remarks>
        /// <para>This operation is both time and memory efficient. It does not convert every element
        /// into a new list, but instead simply instantiates a converter that implements the list
        /// operations directly against the already existing storage.</para>
        /// <para>The type specified by <typeparamref name="T"/> need not be the same as the
        /// underlying stored type of the column. As long as the exposed values can be cast or
        /// converted to the specified type, the typed list returned will work as required.</para>
        /// </remarks>
        public IReadOnlyList<T> As<T> () {
            // If the requested column is of the requested type, expose it directly.
            IReadOnlyList<T> typedColumn = column as IReadOnlyList<T>;
            if (typedColumn != null) {
                return (new TypedFrameColumn<T>(column.Name, typedColumn, map));
            } else {
                // If the requested column is cast-able to the requested type, cast it.
                // Also, allow casts from null-able where the underlying type is cast-able.
                TypeInfo targetTypeInfo = typeof(T).GetTypeInfo();
                if (targetTypeInfo.IsAssignableFrom(column.StorageType.GetTypeInfo()) ||
                    targetTypeInfo.IsAssignableFrom(Nullable.GetUnderlyingType(column.StorageType)?.GetTypeInfo())) {
                    return (new CastFrameColumn<T>(column, map));
                } else {
                    // If nothing else works, use convert.
                    return (new ConvertedFrameColumn<T>(column, map));
                }
            }
        }

        /// <summary>
        /// Returns a transformation of the column as a typed list.
        /// </summary>
        /// <typeparam name="TIn">The input type of the transformation.</typeparam>
        /// <typeparam name="TOut">The output type of the transformation.</typeparam>
        /// <param name="transformation">The transformation.</param>
        /// <returns>An object that allows a transformation of the column to be accessed as a typed list.</returns>
        public IReadOnlyList<TOut> As<TIn, TOut>(Func<TIn, TOut> transformation) {
            IReadOnlyList<TIn> inputColumn = this.As<TIn>();
            return (new TransformedList<TIn, TOut>(this.Name, inputColumn, transformation));
        }

        IEnumerator IEnumerable.GetEnumerator() {
            return (((IEnumerable) column).GetEnumerator());
        }

    }

    internal class TypedFrameColumn<T> : IReadOnlyList<T>, INamed {

        public TypedFrameColumn(string name, IReadOnlyList<T> column, List<int> map) {
            Debug.Assert(name != null);
            Debug.Assert(column != null);
            Debug.Assert(map != null);
            this.name = name;
            this.column = column;
            this.map = map;
        }

        private readonly string name;
        private readonly IReadOnlyList<T> column;
        private readonly List<int> map;

        public string Name {
            get {
                return (name);
            }
        }

        public T this[int index] {
            get {
                return (column[map[index]]);
            }
        }

        public int Count {
            get {
                return (map.Count);
            }
        }

        public IEnumerator<T> GetEnumerator () {
            foreach (int index in map) {
                yield return (column[index]);
            }
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (((IEnumerable<T>) this).GetEnumerator());
        }
    }

    internal class CastFrameColumn<T> : IReadOnlyList<T>, INamed {

        public CastFrameColumn (NamedList column, List<int> map) {
            Debug.Assert(column != null);
            Debug.Assert(map != null);
            this.column = column;
            this.map = map;
        }

        private readonly NamedList column;
        private readonly List<int> map;

        public string Name {
            get {
                return (column.Name);
            }
        }

        public T this[int index] {
            get {
                object value = column.GetItem(map[index]);
                return ((T) value);
            }
        }

        public int Count {
            get {
                return (map.Count);
            }
        }

        public IEnumerator<T> GetEnumerator () {
            foreach (int index in map) {
                object value = column.GetItem(index);
                yield return ((T) value);
            }
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (((IEnumerable<T>) this).GetEnumerator());
        }

    }

    internal class ConvertedFrameColumn<T> : IReadOnlyList<T>, INamed {

        public ConvertedFrameColumn(NamedList column, List<int> map) {
            Debug.Assert(column != null);
            Debug.Assert(map != null);
            this.column = column;
            this.map = map;
        }

        private readonly NamedList column;
        private readonly List<int> map;

        public string Name {
            get {
                return (column.Name);
            }
        }

        public T this[int index] {
            get {
                object value = column.GetItem(map[index]);
                return ((T) Convert.ChangeType(value, typeof(T)));
            }
        }

        public int Count {
            get {
                return (map.Count);
            }
        }

        public IEnumerator<T> GetEnumerator () {
            foreach (int index in map) {
                object value = column.GetItem(index);
                yield return ((T) Convert.ChangeType(value, typeof(T)));
            }
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (((IEnumerable<T>) this).GetEnumerator());
        }

    }

    internal class TransformedList<TIn, TOut> : IReadOnlyList<TOut>, INamed {

        public TransformedList(string name, IReadOnlyList<TIn> input, Func<TIn, TOut> transformation) {
            Debug.Assert(name != null);
            Debug.Assert(input != null);
            Debug.Assert(transformation != null);
            this.name = name;
            this.input = input;
            this.transformation = transformation;
        }

        private readonly string name;
        private readonly IReadOnlyList<TIn> input;
        private readonly Func<TIn, TOut> transformation;

        public string Name {
            get {
                return (name);
            }
        }

        public int Count {
            get {
                return (input.Count);
            }
        }

        public TOut this[int index] {
            get {
                return (transformation(input[index]));
            }
        }

        public IEnumerator<TOut> GetEnumerator () {
            foreach (TIn entry in input) {
                yield return (transformation(entry));
            }
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (((IEnumerable<TOut>) this).GetEnumerator());
        }
    }

    internal interface INamed {

        string Name { get; }

    }

}
