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
    public abstract class DataList /*: IReadOnlyDataList<object> */
    {

        internal DataList (string name) {
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

        /// <summary>
        /// Gets the number of entries in the list.
        /// </summary>
        public abstract int Count { get; }

        internal abstract object GetItem (int index);

        internal abstract void SetItem (int index, object value);

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
        internal static DataList Create (string name, Type type) {
            Debug.Assert(name != null);
            Debug.Assert(type != null);
            Type listType = genericListType.MakeGenericType(type);
            object list = Activator.CreateInstance(listType, name);
            return ((DataList) list);
        }

        private static readonly Type genericListType = typeof(DataList<>);

    }

    /// <summary>
    /// Represents a data list of a particular stored type.
    /// </summary>
    /// <typeparam name="T">The type of value in the list.</typeparam>
    public class DataList<T> : DataList, IReadOnlyDataList<T>, IList<T> {


        /// <summary>
        /// Initializes a new data list with the given name.
        /// </summary>
        /// <param name="name">The name of the data list.</param>
        public DataList (string name) : base(name) {
            storage = new List<T>();
        }

        /// <summary>
        /// Initializes a new data list with the given storage.
        /// </summary>
        /// <param name="name">The name of the data list.</param>
        /// <param name="storage">The stored list.</param>
        public DataList (string name, List<T> storage) : base(name) {
            this.storage = storage;
        }

        /// <summary>
        /// Initializes a new data list with the given values.
        /// </summary>
        /// <param name="name">The name of the data list.</param>
        /// <param name="values">The values to be stored.</param>
        public DataList (string name, IEnumerable<T> values) : base(name) {
            this.storage = values.ToList();
        }

        /// <summary>
        /// Initalizes a new data list with the given values.
        /// </summary>
        /// <param name="name">The name of the data list.</param>
        /// <param name="values">The values to be stored.</param>
        public DataList (string name, params T[] values) : this(name, (IEnumerable<T>) values) {

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

        internal override object GetItem (int index) {
            return (storage[index]);
        }

        internal override void SetItem (int index, object value) {
            storage[index] = (T) value;
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

        // IList (mutable) operations

        /// <summary>
        /// Removes all values from the data list.
        /// </summary>
        public void Clear () {
            storage.Clear();
        }

        /// <summary>
        /// Determines whether the data list contains the given value.
        /// </summary>
        /// <param name="value">The value to search for.</param>
        /// <returns>True if the data list contains the given value, otherwise false.</returns>
        public bool Contains (T value) {
            return (storage.Contains(value));
        }

        void ICollection<T>.CopyTo (T[] array, int start) {
            storage.CopyTo(array, start);
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

        void IList<T>.Insert (int index, T value) {
            InsertAt(index, value);
        }

        bool ICollection<T>.IsReadOnly {
            get {
                return (false);
            }
        }

        /// <summary>
        /// Removes the first occurance of the given value.
        /// </summary>
        /// <param name="value">The value to remove.</param>
        /// <returns><see langword="true"/> if the value was found and removed, otherwise <see langword="false"/>.</returns>
        public bool Remove (T value) {
            return (storage.Remove(value));
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
