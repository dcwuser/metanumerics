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

            string cell = text.Substring(startIndex, index - startIndex);
            cell = cell.Trim();
            if (cell.Length > 0) {
                if (cell[0] == escape) {
                    if (cell.Length < 2 || cell[cell.Length - 1] != escape) throw new FormatException();
                    cell = cell.Substring(1, cell.Length - 2);
                    cell = cell.Replace($"{escape}{escape}", Char.ToString(escape));
                } else {
                    if (cell.LastIndexOfAny(delimiterOrEscape) >= 0) throw new FormatException();
                }
            }
            /*
            int endIndex = index - 1;
            while (startIndex < endIndex && Char.IsWhiteSpace(text[startIndex])) startIndex++;
            while (startIndex < endIndex && Char.IsWhiteSpace(text[endIndex])) endIndex--;

            string cell;
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
            */
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

}
