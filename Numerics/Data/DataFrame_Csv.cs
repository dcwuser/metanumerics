using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;

namespace Meta.Numerics.Data {

    internal static class CsvHelper {

        private static char delimiter = ',';
        private static char escape = '"';
        private static char[] delimiterOrEscape = new char[] { delimiter, escape };

        public static List<string> ReadCells (string text) {
            Debug.Assert(text != null);

            List<string> cells = new List<string>();
            int startIndex = 0;
            bool escaped = false;
            for (int index = 0; index < text.Length; index++) {
                if (text[index] == delimiter) {
                    if (!escaped) {
                        string cell = ReadCell(text, startIndex, index);
                        cells.Add(cell);
                        startIndex = index + 1;
                    }
                } else if (text[index] == escape) {
                    escaped = !escaped;
                }
            }
            if (escaped) throw new FormatException();
            if (startIndex == text.Length) {
                cells.Add(String.Empty);
            } else if (startIndex < text.Length) {
                string value = ReadCell(text, startIndex, text.Length);
                cells.Add(value);
            }

            return (cells);
        }

        private static string ReadCell (string text, int startIndex, int index) {
            Debug.Assert(text != null);
            Debug.Assert((0 <= startIndex) && (startIndex < text.Length));
            Debug.Assert((startIndex <= index) && (index <= text.Length));

            string cell;
            int endIndex = index - 1;
            if (text[startIndex] == escape) {
                startIndex++;
                if (text[endIndex] != escape) throw new FormatException();
                endIndex--;
                cell = text.Substring(startIndex, endIndex - startIndex + 1);
                cell = cell.Replace($"{escape}{escape}", Char.ToString(escape));
            } else {
                cell = text.Substring(startIndex, endIndex - startIndex + 1);
                int forbiddenIndex = cell.IndexOfAny(delimiterOrEscape);
                if (forbiddenIndex >= 0) throw new FormatException();
            }
            return (cell);
        }

        public static string WriteCells (IEnumerable<object> cells) {
            Debug.Assert(cells != null);

            StringBuilder line = new StringBuilder();
            bool first = true;
            foreach (object cell in cells) {
                if (first) {
                    first = false;
                } else {
                    line.Append(',');
                }
                line.Append(WriteCell(cell));
            }
            return (line.ToString());
        }

        private static string WriteCell (object cell) {
            if (cell == null) {
                return (String.Empty);
            } else {
                string text = cell.ToString();
                int forbiddenIndex = text.IndexOfAny(delimiterOrEscape);
                if (forbiddenIndex >= 0) {
                    text = text.Replace(Char.ToString(escape), $"{escape}{escape}");
                    text = String.Format("\"{0}\"", text);
                }
                return (text);
            }
        }

    }

    public sealed partial class DataFrame
    {

        /// <summary>
        /// Creates a new data frame by reading values from a file of comma-seperated values.
        /// </summary>
        /// <param name="reader"></param>
        /// <returns>A new data frame with data from the file.</returns>
        /// <remarks>The column names are taken from the first line of the file.</remarks>
        public static DataFrame ReadCsv (TextReader reader) {
            if (reader == null) throw new ArgumentNullException(nameof(reader));

            DataList<string>[] textColumns;
            DataAdaptor[] headers;
            ReadCsvAsStrings(reader, out textColumns, out headers);

            DataList[] columns = new DataList[headers.Length];
            for (int columnIndex = 0; columnIndex < columns.Length; columnIndex++) {
                DataAdaptor header = headers[columnIndex];
                if (header.TypeCandidates.Count == 0) {
                    columns[columnIndex] = textColumns[columnIndex];
                } else {
                    TypeAdaptor adaptor = header.TypeCandidates.First.Value;
                    DataList column = adaptor.CreateStorage(textColumns[columnIndex].Name, header.IsNullable);
                    foreach (string textValue in textColumns[columnIndex]) {
                        if (textValue == null) {
                            column.AddItem(null);
                        } else {
                            object value = adaptor.Parse(textValue);
                            column.AddItem(value);
                        }
                    }
                    columns[columnIndex] = column;
                }
            }

            DataFrame frame = new DataFrame(columns);
            return (frame);
        }

        private static void ReadCsvAsStrings (TextReader reader, out DataList<string>[] columns, out DataAdaptor[] headers) {
            if (reader == null) throw new ArgumentNullException(nameof(reader));

            string firstline = reader.ReadLine();
            if (firstline == null) {
                columns = null;
                headers = null;
                return;
            }

            List<string> names = CsvHelper.ReadCells(firstline);
            int count = names.Count;

            columns = new DataList<string>[names.Count];
            headers = new DataAdaptor[names.Count];
            for (int columnIndex = 0; columnIndex < columns.Length; columnIndex++) {
                columns[columnIndex] = new DataList<string>(names[columnIndex]);
                headers[columnIndex] = new DataAdaptor();
            }

            while (true) {
                string line = reader.ReadLine();
                if (line == null) break;
                List<string> cells = CsvHelper.ReadCells(line);
                if (cells.Count != count) throw new FormatException();
                for (int columnIndex = 0; columnIndex < count; columnIndex++) {
                    string cell = cells[columnIndex];
                    if (String.IsNullOrEmpty(cell)) {
                        headers[columnIndex].IsNullable = true;
                        columns[columnIndex].Add(null);
                    } else {
                        DataAdaptor header = headers[columnIndex];
                        LinkedListNode<TypeAdaptor> adaptorNode = header.TypeCandidates.First;
                        while (adaptorNode != null) {
                            if (!adaptorNode.Value.IsParsable(cell)) header.TypeCandidates.Remove(adaptorNode);
                            adaptorNode = adaptorNode.Next;
                        }
                        columns[columnIndex].Add(cell);
                    }
                }
            }

        }

    }

    // Type adaptors parse text into objects of a particular type.

    internal abstract class TypeAdaptor
    {
        public abstract bool IsParsable(string text);

        public abstract object Parse(string text);

        public abstract DataList CreateStorage(string name, bool nullable);

    }

    internal class DoubleAdaptor : TypeAdaptor
    {
        public override bool IsParsable(string text)
        {
            double value;
            bool isParsable = Double.TryParse(text, out value);
            return (isParsable);
        }

        public override object Parse(string text)
        {
            return (Double.Parse(text));
        }

        public override DataList CreateStorage (string name, bool nullable) {
            if (nullable) {
                return (new DataList<double?>(name));
            } else {
                return (new DataList<double>(name));
            }
        }

    }

    internal class IntAdaptor : TypeAdaptor
    {
        public override bool IsParsable(string text)
        {
            int value;
            bool isParsable = Int32.TryParse(text, out value);
            return (isParsable);
        }

        public override object Parse(string text)
        {
            return (Int32.Parse(text));
        }

        public override DataList CreateStorage (string name, bool nullable) {
            if (nullable) {
                return (new DataList<int?>(name));
            } else {
                return (new DataList<int>(name));
            }
        }
    }

    internal class DateTimeAdaptor : TypeAdaptor
    {

        public override bool IsParsable(string text)
        {
            DateTime value;
            bool isParsable = DateTime.TryParse(text, out value);
            return (isParsable);
        }

        public override object Parse(string text)
        {
            return (DateTime.Parse(text));
        }

        public override DataList CreateStorage(string name, bool nullable)
        {
            if (nullable) {
                return (new DataList<DateTime?>(name));
            } else {
                return (new DataList<DateTime>(name));
            }
        }

    }

    internal class TimeSpanAdaptor : TypeAdaptor {

        public override bool IsParsable (string text) {
            Debug.Assert(text != null);
            // TimeSpan's parse will accept pure ints or doubles as timespans.
            // Since we don't want a coulumn of ints or doubles to become timespans,
            // we require a colon before we even try.
            if (!text.Contains(":")) return (false);
            TimeSpan value;
            bool isParsalbe = TimeSpan.TryParse(text, out value);
            return (isParsalbe);
        }

        public override object Parse (string text) {
            TimeSpan value = TimeSpan.Parse(text);
            return (value);
        }

        public override DataList CreateStorage (string name, bool nullable) {
            if (nullable) {
                return (new DataList<TimeSpan?>(name));
            } else {
                return (new DataList<TimeSpan>(name));
            }
        }

    }

    // This type is used to keep track of which types are compatible with the text
    // representations that have been read in for a column.

    internal class DataAdaptor
    {
        private static TypeAdaptor[] knownTypes = new TypeAdaptor[] {
            new DateTimeAdaptor(), /* new TimeSpanAdaptor(), */ new IntAdaptor(), new DoubleAdaptor()
        };

        public bool IsNullable = false;

        internal LinkedList<TypeAdaptor> TypeCandidates = new LinkedList<TypeAdaptor>(knownTypes);

    }

}
