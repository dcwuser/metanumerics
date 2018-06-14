using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Contains types for the statistical analysis of data.
    /// </summary>
    /// <remarks>
    /// <para>This namespace contains types for doing basic and advanced statistics. It contains
    /// APIs for finding moments and percentiles, measuring associations, fitting to models, and
    /// other statistical operations.
    /// </para>
    /// <para>The central class for operating on samples consisting of independent measurements of a single variable is the
    /// <see cref="Univariate"/> class. The central class for operating on samples of independent bivariate
    /// data (i.e. paired measurements) is the <see cref="Bivariate"/> class. The central class
    /// for operating on samples of independent measurements of multiple variables is the <see cref="Multivariate"/>
    /// class. The central class for operating on time series data is the <see cref="Series"/> class.</para>
    /// <para>All of these central classes consist of static methods that accept one or more columns of data. Each column
    /// can be of any type that implements the appropriate collection interface
    /// (e.g. <see cref="IReadOnlyList{T}"/>).
    /// Many of the methods are extension methods, so they effectively become instance methods on all such types.</para>
    /// <para>Some of the classes in this namespace are left over from earlier versions of Meta.Numerics
    /// which required users to store each kind of data in a particular storage class. Examples
    /// of these storage classes include <see cref="Sample"/>, <see cref="BivariateSample"/>,
    /// <see cref="MultivariateSample"/>, and <see cref="TimeSeries"/>. These storage classes
    /// each expose methods appropriate for the analysis of a particular type of data. The
    /// advantage of such a system is that it makes immediately clear to the user which
    /// methods are appropriate for which types of data. The disadvantage is that it requires
    /// users to transfer their data into our containers before it can be analyzed. You can still
    /// use these classes, if you prefer, but essentially all of their functionality is also exposed
    /// in the new central, static classes that can be applied to any appropriate collection.</para>
    /// </remarks>
    [CompilerGenerated]
    internal static class NamespaceDoc {
    }
}
