using System;
using System.Collections;
using System.Collections.Generic;

namespace Meta.Numerics.Data
{

    /// <summary>
    /// One column of a data frame.
    /// </summary>
    public sealed class FrameColumn : IEnumerable
    {
        internal FrameColumn(FrameView frame, int c)
        {
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
        public string Name
        {
            get
            {
                return (column.Name);
            }
        }

        /// <summary>
        /// Gets the type of data stored in the column.
        /// </summary>
        public Type StorageType
        {
            get
            {
                return (column.StorageType);
            }
        }

        /// <summary>
        /// Gets the number of values in the column.
        /// </summary>
        public int Count
        {
            get
            {
                return (map.Count);
            }
        }

        /// <summary>
        /// Gets the value at the given index.
        /// </summary>
        /// <param name="r">The (zero-based) index.</param>
        /// <returns>The value at index <paramref name="r"/>.</returns>
        public object this[int r]
        {
            get
            {
                return (column.GetItem(map[r]));
            }
        }

        /// <summary>
        /// Returns the column as a typed list.
        /// </summary>
        /// <typeparam name="T">The type of the list.</typeparam>
        /// <returns>An object that allows the column to be accessed as a typed list.</returns>
        /// <remarks>
        /// <para>This operation is both time and memory efficient. It does not convert every element
        /// into a new list, but instead simply instantiates a converter that implements the list
        /// operations directly against the already existing storage.</para>
        /// </remarks>
        public IReadOnlyList<T> As<T> () {
            NamedList<T> typedColumn = column as NamedList<T>;
            if (typedColumn != null) {
                return (new TypedFrameColumn<T>(typedColumn, map));
            } else {
                return (new ConvertedFrameColumn<T>(column, map));
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
            return (new TransformedList<TIn, TOut>(inputColumn, transformation));
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return (((IEnumerable) column).GetEnumerator());
        }

    }

    internal class TypedFrameColumn<T> : IReadOnlyList<T> {

        public TypedFrameColumn(NamedList<T> column, List<int> map) {
            this.column = column;
            this.map = map;
        }

        private NamedList<T> column;
        private List<int> map;

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

    internal class ConvertedFrameColumn<T> : IReadOnlyList<T> {

        public ConvertedFrameColumn(NamedList column, List<int> map) {
            this.column = column;
            this.map = map;
        }

        private NamedList column;
        private List<int> map;

        public T this[int index] {
            get {
                return ((T) column.GetItem(map[index]));
            }
        }

        public int Count {
            get {
                return (map.Count);
            }
        }

        public IEnumerator<T> GetEnumerator () {
            foreach (int index in map) {
                yield return ((T) column.GetItem(index));
            }
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (((IEnumerable<T>) this).GetEnumerator());
        }

    }

    internal class TransformedList<TIn, TOut> : IReadOnlyList<TOut> {

        public TransformedList(IReadOnlyList<TIn> input, Func<TIn, TOut> transformation) {
            this.input = input;
            this.transformation = transformation;
        }

        private readonly IReadOnlyList<TIn> input;

        private readonly Func<TIn, TOut> transformation;

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

}
