using System;
using System.Collections;
using System.Collections.Generic;

namespace Meta.Numerics.Statistics {

    // Since this is internal and we don't use all IDictionary methods, there is no need for
    // us to implement them all.

    internal class NullableDictionary<K, V> : /* IDictionary<K, V>, */ IEnumerable<KeyValuePair<K, V>> {

        public NullableDictionary () {
            dictionary = new Dictionary<K, V>();
        }

        private readonly Dictionary<K, V> dictionary;

        private bool nullKey;

        private V nullValue;

        public IReadOnlyCollection<K> Keys {
            get {
                return (new KeyCollection(this));
            }
        }

        public int Count {
            get {
                int count = dictionary.Count;
                if (nullKey) count++;
                return (count);
            }
        }

        public V this[K key] {
            get {
                if (key == null) {
                    if (nullKey) {
                        return (nullValue);
                    } else {
                        throw new KeyNotFoundException();
                    }
                } else {
                    return (dictionary[key]);
                }
            }
            set {
                if (key == null) {
                    nullKey = true;
                    nullValue = value;
                } else {
                    dictionary[key] = value;
                }
            }
        }

        public void Add (K key, V value) {
            if (key == null) {
                if (nullKey) {
                    throw new ArgumentException();
                } else {
                    nullKey = true;
                    nullValue = value;
                }
            } else {
                dictionary.Add(key, value);
            }
        }

        public bool ContainsKey (K key) {
            if (key == null) {
                return (nullKey);
            } else {
                return (dictionary.ContainsKey(key));
            }
        }

        public bool TryGetValue (K key, out V value) {
            if (key == null) {
                if (nullKey) {
                    value = nullValue;
                    return (true);
                } else {
                    value = default(V);
                    return (false);
                }
            } else {
                return (dictionary.TryGetValue(key, out value));
            }
        }

        public IEnumerator<KeyValuePair<K, V>> GetEnumerator () {
            foreach (KeyValuePair<K, V> item in dictionary) {
                yield return item;
            }
            if (nullKey) yield return new KeyValuePair<K, V>(default(K), nullValue);
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (((IEnumerable<KeyValuePair<K, V>>) this).GetEnumerator());
        }

        internal class KeyCollection : IReadOnlyCollection<K> {

            internal KeyCollection (NullableDictionary<K, V> parent) {
                this.parent = parent;
            }

            private NullableDictionary<K, V> parent;

            public int Count {
                get {
                    return (parent.Count);
                }
            }

            public IEnumerator<K> GetEnumerator () {
                foreach (KeyValuePair<K, V> item in parent) {
                    yield return item.Key;
                }
            }

            IEnumerator IEnumerable.GetEnumerator () {
                return (((IEnumerable<K>) this).GetEnumerator());
            }

        }

    }
}
