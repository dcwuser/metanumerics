using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;

namespace Meta.Numerics.Data {

    // DataAdaptor and TypeAdaptors are used internally to auto-recognize data read as strings.
    // For each column of string data of unknown type, instantiate a new DataAdaptor. As data
    // is read, call IsParsable to determine whether the string is compatible with each type
    // candidate. If it is not, remove that TypeAdaptor from the candidates list. If all types
    // a removed, just leave the type as a string.

    internal class DataAdaptor {

        private static readonly TypeParser[] knownTypes = new TypeParser[] {
            new DateTimeParser(), new TimeSpanParser(), new IntParser(), new DoubleParser(), new BoolParser()
        };

        public bool IsNullable = false;

        internal LinkedList<TypeParser> TypeCandidates = new LinkedList<TypeParser>(knownTypes);

        public void TryParse (string text) {
            LinkedListNode<TypeParser> adaptorNode = TypeCandidates.First;
            while (adaptorNode != null) {
                // It's necessary to extract the next node first because calling Remove sets Next to null
                LinkedListNode<TypeParser> nextAdaptorNode = adaptorNode.Next;
                bool parsable = adaptorNode.Value.IsParsable(text);
                if (!parsable) TypeCandidates.Remove(adaptorNode);
                adaptorNode = nextAdaptorNode;
            }
        }

    }

    internal abstract class TypeParser {

        public abstract bool IsParsable(string text);

        public abstract object Parse(string text);

        public abstract NamedList CreateStorage(string name, bool nullable);

    }

    internal abstract class TypeParser<T> : TypeParser where T : struct {

        public override NamedList CreateStorage(string name, bool nullable) {
            if (nullable) {
                return (new NamedList<T?>(name));
            } else {
                return (new NamedList<T>(name));
            }
        }

    }

    internal class DoubleParser : TypeParser<double> {

        public override bool IsParsable(string text) {
            double value;
            bool isParsable = Double.TryParse(text, out value);
            return (isParsable);
        }

        public override object Parse(string text) {
            return (Double.Parse(text, NumberFormatInfo.InvariantInfo));
        }

    }

    internal class IntParser : TypeParser<int> {
        public override bool IsParsable(string text) {
            int value;
            bool isParsable = Int32.TryParse(text, NumberStyles.Integer, NumberFormatInfo.InvariantInfo, out value);
            return (isParsable);
        }

        public override object Parse(string text) {
            return (Int32.Parse(text, NumberFormatInfo.InvariantInfo));
        }

    }

    internal class DateTimeParser : TypeParser<DateTime> {

        public override bool IsParsable(string text) {
            DateTime value;
            bool isParsable = DateTime.TryParse(text, DateTimeFormatInfo.InvariantInfo, DateTimeStyles.None, out value);
            return (isParsable);
        }

        public override object Parse(string text) {
            return (DateTime.Parse(text, DateTimeFormatInfo.InvariantInfo));
        }

    }

    internal class TimeSpanParser : TypeParser<TimeSpan> {

        public override bool IsParsable(string text)
        {
            Debug.Assert(text != null);
            // TimeSpan's parse will accept pure ints or doubles as timespans.
            // Since we don't want a coulumn of ints or doubles to become timespans,
            // we require a colon before we even try.
            if (!text.Contains(":")) return (false);
            TimeSpan value;
            bool isParsalbe = TimeSpan.TryParse(text, out value);
            return (isParsalbe);
        }

        public override object Parse(string text) {
            TimeSpan value = TimeSpan.Parse(text);
            return (value);
        }

    }

    internal class BoolParser : TypeParser<bool> {

        public override bool IsParsable(string text) {
            bool value;
            bool isParsable = Boolean.TryParse(text, out value);
            return (isParsable);
        }

        public override object Parse(string text) {
            return (Boolean.Parse(text));
        }

    }
}
