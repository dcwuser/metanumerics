using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Data {

    // DataAdaptor and TypeAdaptors are used internally to auto-recognize data read as strings.
    // For each column of string data of unknown type, instantiate a new DataAdaptor. As data
    // is read, call IsParsable to determine whether the string is compatible with each type
    // candidate. If it is not, remove that TypeAdaptor from the candidates list. If all types
    // a removed, just leave the type as a string.

    internal class DataAdaptor {

        private static TypeAdaptor[] knownTypes = new TypeAdaptor[] {
            new DateTimeAdaptor(), new TimeSpanAdaptor(), new IntAdaptor(), new DoubleAdaptor(), new BooleanAdaptor()
        };

        public bool IsNullable = false;

        internal LinkedList<TypeAdaptor> TypeCandidates = new LinkedList<TypeAdaptor>(knownTypes);

    }

    internal abstract class TypeAdaptor {

        public abstract bool IsParsable(string text);

        public abstract object Parse(string text);

        public abstract DataList CreateStorage(string name, bool nullable);

    }

    internal abstract class TypeAdaptor<T> : TypeAdaptor where T : struct {

        public override DataList CreateStorage(string name, bool nullable) {
            if (nullable) {
                return (new DataList<T?>(name));
            } else {
                return (new DataList<T>(name));
            }
        }

    }

    internal class DoubleAdaptor : TypeAdaptor<double> {

        public override bool IsParsable(string text) {
            double value;
            bool isParsable = Double.TryParse(text, out value);
            return (isParsable);
        }

        public override object Parse(string text) {
            return (Double.Parse(text));
        }

    }

    internal class IntAdaptor : TypeAdaptor<int> {
        public override bool IsParsable(string text) {
            int value;
            bool isParsable = Int32.TryParse(text, out value);
            return (isParsable);
        }

        public override object Parse(string text) {
            return (Int32.Parse(text));
        }

    }

    internal class DateTimeAdaptor : TypeAdaptor<DateTime> {

        public override bool IsParsable(string text) {
            DateTime value;
            bool isParsable = DateTime.TryParse(text, out value);
            return (isParsable);
        }

        public override object Parse(string text) {
            return (DateTime.Parse(text));
        }

    }

    internal class TimeSpanAdaptor : TypeAdaptor<TimeSpan> {

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

    internal class BooleanAdaptor : TypeAdaptor<bool> {

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
