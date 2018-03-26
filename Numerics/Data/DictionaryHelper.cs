using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Reflection;

namespace Meta.Numerics.Data {

    internal static class DictionaryHelper {

        public static List<NamedList> ReadDictionaries (IEnumerable<IReadOnlyDictionary<string, object>> dictionaries) {

            Debug.Assert(dictionaries != null);

            // Iterate through the dictionaries, creating header objects that contain the un-cast values
            // and some information about them.
            List<DictionaryColumn> headers = null;
            foreach (IReadOnlyDictionary<string, object> dictionary in dictionaries) {

                // From the first row, create the headers list based on key names.
                if (headers == null) {
                    headers = new List<DictionaryColumn>(dictionary.Count);
                    foreach (string key in dictionary.Keys) {
                        DictionaryColumn header = new DictionaryColumn() {
                            Name = key, IsNullable = false, Type = null, Data = new List<object>()
                        };
                        headers.Add(header);
                    }
                }

                if (dictionary.Count != headers.Count) throw new InvalidOperationException();

                // For all rows, check for null, record the type if we haven't found it yet, and store the value.
                for (int i = 0; i < headers.Count; i++) {
                    DictionaryColumn header = headers[i];
                    object value = dictionary[header.Name];
                    if (value == null) {
                        header.IsNullable = true;
                    } else {
                        if (header.Type == null) header.Type = value.GetType();
                    }
                    header.Data.Add(value);
                }
            }

            // Arrange the columns into named lists of the appropriate type
            List<NamedList> columns = new List<NamedList>(headers.Count);
            foreach (DictionaryColumn header in headers) {
                NamedList column;
                if (header.Type == null) {
                    // If no non-null value was ever found, we can't infer a type, so just make an object-column.
                    column = new NamedList<object>(header.Name, header.Data);
                } else {
                    // Based on null-ability and observed type, create the appropriate storage.
                    Type type = header.Type;
                    if (header.IsNullable && type.GetTypeInfo().IsValueType) type = typeof(Nullable<>).MakeGenericType(type);
                    column = NamedList.Create(header.Name, type);
                    // Copy the objects into the storage, which will cast them to the storage type.
                    foreach (object value in header.Data) column.AddItem(value);
                }
                columns.Add(column);
            }
            Debug.Assert(columns.Count == headers.Count);

            return (columns);

        }

        // This class is used internally as we parse dictionaries into columns.
        private class DictionaryColumn {

            public string Name;

            public bool IsNullable;

            public Type Type;

            public List<object> Data;

        }

    }

}
