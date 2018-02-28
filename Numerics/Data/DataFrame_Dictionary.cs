using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Data
{
    public partial class DataFrame
    {

        /// <summary>
        /// Constructs a new data frame from a sequence of dictionaries.
        /// </summary>
        /// <param name="dictionaries"></param>
        /// <returns></returns>
        public static DataFrame FromDictionaries (IEnumerable<IReadOnlyDictionary<string, object>> dictionaries)
        {
            if (dictionaries == null) throw new ArgumentNullException(nameof(dictionaries));

            DataFrame frame = null;
            foreach(IReadOnlyDictionary<string, object> dictionary in dictionaries)
            {
                if (frame == null)
                {
                    List<ColumnDefinition> headers = new List<ColumnDefinition>();
                    foreach(KeyValuePair<string, object> entry in dictionary)
                    {
                        ColumnDefinition header = new ColumnDefinition(entry.Key, entry.Value.GetType());
                        headers.Add(header);
                        frame = new DataFrame(headers.ToArray());
                    }
                }

                foreach(KeyValuePair<string, object> entry in dictionary)
                {
                    int columnIndex = frame.columnMap[entry.Key];
                    DataList column = frame.columns[columnIndex];
                    if ((entry.Value == null) && !column.IsNullable) {
                        // fix column
                        Type nullableType = typeof(Nullable<>).MakeGenericType(column.StorageType);
                        DataList nullableColumn = DataList.Create(entry.Key, nullableType);
                        for (int i = 0; i < column.Count; i++)
                        {
                            nullableColumn.AddItem(column.GetItem(i));
                        }
                        frame.columns[columnIndex] = nullableColumn;
                        column = nullableColumn;
                    }
                    int rowIndex = column.AddItem(entry.Value);
                }
                frame.map.Add(frame.map.Count);
            }
            return (frame);
        }

    }
}
