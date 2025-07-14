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
    /// <para>The data storage classes Sample, BivariateSample, MultivariateSample, and TimeSeries,
    /// which appeared in earlier versions of Meta.Numerics, have been removed.
    /// All the summary statistics, statistical tests, fits, and other functionality that they provided
    /// can now be accessed using the methods of the <see cref="Univariate"/>, <see cref="Bivariate"/>,
    /// <see cref="Multivariate"/> and <see cref="Series"/> classes. The new methods can be applied
    /// to data in arbitrary collection types, so specialized storage types are no longer necessary.
    /// We expect to remove the storage types in a future release.</para>
    /// </remarks>
    [CompilerGenerated]
    internal static class NamespaceDoc {
    }
}
