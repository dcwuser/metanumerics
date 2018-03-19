using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Reflection;

namespace Meta.Numerics.Data
{

    /// <summary>
    /// Represents a data list of any stored data type.
    /// </summary>
    internal abstract class NamedList /*: IReadOnlyDataList<object> */
    {

        internal NamedList (string name) {
            if (name == null) throw new ArgumentNullException(nameof(name));
            this.name = name;
        }

        private string name;

        /// <summary>
        /// Gets the name of the list.
        /// </summary>
        public virtual string Name {
            get {
                return (name);
            }
        }

        /// <summary>
        /// Gets the type of values stored in the list.
        /// </summary>
        public abstract Type StorageType { get; }

        /// <summary>
        /// Gets a value indicating whether the entries can have null values.
        /// </summary>
        public virtual bool IsNullable {
            get {
                Type type = StorageType;
                if (type.GetTypeInfo().IsValueType) {
                    return (Nullable.GetUnderlyingType(type) != null);
                } else {
                    return (true);
                }
            }
        }

        internal abstract bool IsComputed { get; }

        /// <summary>
        /// Gets the number of entries in the list.
        /// </summary>
        public abstract int Count { get; }

        internal abstract object GetItem (int index);

        internal abstract void SetItem (int index, object value);

        public abstract void Clear ();

        internal abstract int AddItem (object value);

        /*
        object IReadOnlyList<object>.this[int index]
        {
            get
            {
                return (GetItem(index));
            }
        }

        IEnumerator<object> IEnumerable<object>.GetEnumerator()
        {
            for (int index = 0; index < this.Count; index++)
            {
                yield return (GetItem(index));
            }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return (((IEnumerable<object>)this).GetEnumerator());
        }
        */
        internal static NamedList Create (string name, Type type) {
            Debug.Assert(name != null);
            Debug.Assert(type != null);
            Type listType = genericListType.MakeGenericType(type);
            object list = Activator.CreateInstance(listType, name);
            return ((NamedList) list);
        }

        private static readonly Type genericListType = typeof(NamedList<>);

    }

    /// <summary>
    /// Represents a data list of a particular stored type.
    /// </summary>
    /// <typeparam name="T">The type of value in the list.</typeparam>
    internal class NamedList<T> : NamedList, /*IReadOnlyDataList<T>, IList<T>, ICollection<T>,*/ IEnumerable<T>, IEnumerable {

        /// <summary>
        /// Initializes a new data list with the given name.
        /// </summary>
        /// <param name="name">The name of the data list.</param>
        public NamedList (string name) : base(name) {
            storage = new List<T>();
        }

        /// <summary>
        /// Initializes a new data list with the given storage.
        /// </summary>
        /// <param name="name">The name of the data list.</param>
        /// <param name="storage">The stored list.</param>
        public NamedList (string name, List<T> storage) : base(name) {
            this.storage = storage;
        }


        private readonly List<T> storage;

        /// <summary>
        /// Gets the number of values in the list.
        /// </summary>
        public override int Count {
            get {
                return (storage.Count);
            }
        }

        /// <summary>
        /// Gets the type of stored data.
        /// </summary>
        public override Type StorageType {
            get {
                return (typeof(T));
            }
        }

        internal override bool IsComputed {
            get {
                return (false);
            }
        }

        internal override object GetItem (int index) {
            return (storage[index]);
        }

        internal override void SetItem (int index, object value) {
            storage[index] = (T) value;
        }

        public override void Clear () {
            storage.Clear();
        }

        internal override int AddItem (object value) {
            int index = storage.Count;
            storage.Add((T) value);
            return (index);
        }

        /// <summary>
        /// Gets the value at the given index.
        /// </summary>
        /// <param name="index">The (zero-based) index.</param>
        /// <returns>The value at the index.</returns>
        public virtual T this[int index] {
            get {
                return (storage[index]);
            }
            set {
                storage[index] = value;
            }
        }

        /// <summary>
        /// Adds the given value to the list.
        /// </summary>
        /// <param name="value">The value to add.</param>
        public void Add (T value) {
            //int count = storage.Count;
            storage.Add(value);
            //return (count);
        }

        IEnumerator<T> IEnumerable<T>.GetEnumerator () {
            return (storage.GetEnumerator());
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (storage.GetEnumerator());
        }

        /// <summary>
        /// Determines the first index at which the given value occurs.
        /// </summary>
        /// <param name="value">The value to search for.</param>
        /// <returns>The index at which <paramref name="value"/> first occurs, or -1 if the value
        /// is not in the data list.</returns>
        public int IndexOf (T value) {
            return (storage.IndexOf(value));
        }

        /// <summary>
        /// Inserts the value at the given index.
        /// </summary>
        /// <param name="index">The index at which to insert <paramref name="value"/>.</param>
        /// <param name="value">The value to insert.</param>
        public void InsertAt (int index, T value) {
            storage.Insert(index, value);
        }


        /// <summary>
        /// Removes the value at the given index.
        /// </summary>
        /// <param name="index">The index of the value to remove.</param>
        /// <remarks>
        /// <para>This is a O(n) operation, where n is the number of entries after the given index.</para>
        /// </remarks>
        public void RemoveAt (int index) {
            storage.RemoveAt(index);
        }

    }

}
