using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Data
{
    public partial class FrameTable
    {

        /// <summary>
        /// Constructs a new frame table from a sequence of dictionaries.
        /// </summary>
        /// <param name="dictionaries">A enumerable set of dictionaries, one for each row, whose
        /// keys are the column names and whose values are the cell value of the column for that row.</param>
        /// <returns>A new frame table with data from the dictionaries.</returns>
        /// <remarks><para>The storage type of each column is inferred from the types of objects
        /// encountered are the frame table is constructed.</para></remarks>
        public static FrameTable FromDictionaries (IEnumerable<IReadOnlyDictionary<string, object>> dictionaries)
        {
            if (dictionaries == null) throw new ArgumentNullException(nameof(dictionaries));

            FrameTable frame = null;
            foreach(IReadOnlyDictionary<string, object> dictionary in dictionaries)
            {
                if (dictionary == null) throw new ArgumentNullException(nameof(dictionary));

                // Create columns based on first types encountered
                // Change to a Dictionary<string,DataColumn> to avoid Frame overhead?
                if (frame == null)
                {
                    List<NamedList> columns = new List<NamedList>();
                    foreach(KeyValuePair<string, object> entry in dictionary)
                    {
                        NamedList column = NamedList.Create(entry.Key, entry.Value.GetType());
                        columns.Add(column);
                    }
                    frame = new FrameTable(columns);
                }

                // Add in each row
                // If a new key/value pair appears, this will fail with KeyNotFound.
                // If an expected key/value pair is missing, it should fail too.
                foreach(KeyValuePair<string, object> entry in dictionary)
                {
                    int columnIndex = frame.columnMap[entry.Key];
                    NamedList column = frame.columns[columnIndex];
                    // If we encounter a null in a column which was created as non-nullable, re-create the column as nullable.
                    if ((entry.Value == null) && !column.IsNullable) {
                        Type nullableType = typeof(Nullable<>).MakeGenericType(column.StorageType);
                        NamedList nullableColumn = NamedList.Create(entry.Key, nullableType);
                        for (int i = 0; i < column.Count; i++)
                        {
                            nullableColumn.AddItem(column.GetItem(i));
                        }
                        frame.columns[columnIndex] = nullableColumn;
                        column = nullableColumn;
                    }
                    column.AddItem(entry.Value);
                }
                frame.map.Add(frame.map.Count);
            }
            return (frame);
        }

    }
}
