using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace Meta.Numerics.Data {

    /// <summary>
    /// Contains types used to import and export, filter, and transform data sets.
    /// </summary>
    /// <remarks>
    /// <para>This namespace contains types for data wrangling, i.e. reading data and preparing it for analysis.</para>
    /// <para>Our data wrangling system is based on the data frames concept. Data frames are table-based
    /// data stores with the following characteristics:</para>
    /// <list type="bullet">
    /// <item>Columns are addressable by name.</item>
    /// <item>All data in a column is of the same type, but different columns in the same table can have different types.</item>
    /// <item>Columns can added, removed, computed, and otherwise manipulated.</item>
    /// <item>Rows can be added, removed, filtered, re-ordered, and otherwise manipulated.</item>
    /// <item>Tables can be imported and exported in multiple formats.</item>
    /// </list>
    /// <para>The central class for storing data is the <see cref="FrameTable"/> class.
    /// Methods such a <see cref="FrameTable.FromCsv(System.IO.TextReader)"/> and
    /// <see cref="FrameTable.FromDictionaries(IEnumerable{IReadOnlyDictionary{string, object}})"/>
    /// allow data to be imported into a table from CSV or JSON formatted text files.</para>
    /// <para>The central class for interrogating data is the <see cref="FrameView"/> class.
    /// Each <see cref="FrameTable"/> is itself a view, and most manipulations
    /// of a view produce a new view, which can itself be further manipulated.
    /// When a new <see cref="FrameView"/> is returned, the original
    /// data is not copied, but simply referenced, so these manipulations are typically fast
    /// and have a much smaller memory footprint than the data they expose.</para>
    /// <para>After using the types of this namespace to prepare your data, you can
    /// extract columns from your views as strongly-typed lists (which implement
    /// <see cref="IReadOnlyList{T}"/>, <see cref="IReadOnlyCollection{T}"/>, and
    /// <see cref="IEnumerable{T}"/>) and hand them to any analysis API.
    /// This includes our own statistical analysis APIs in the <see cref="Meta.Numerics.Data"/>
    /// namespace, as well as any APIs from other libraries that operate of lists of strongly typed
    /// collections.</para>
    /// </remarks>
    [CompilerGenerated]
    internal class NamespaceDoc {
    }

}
