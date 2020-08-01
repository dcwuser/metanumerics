
namespace Meta.Numerics.Extended {

    // Used by our parsers to return the result of attempted parsing.

    internal enum ParseResult {
        Success,
        Null,
        Empty,
        Sign,
        Format,
        Overflow
    }
}
