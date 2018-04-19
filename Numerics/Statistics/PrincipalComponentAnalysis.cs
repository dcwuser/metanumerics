using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a principal component analysis.
    /// </summary>
    /// <remarks>
    /// <para>Principal component analysis (PCA) is a form of factor analysis or dimension reduction.
    /// It attempts to identify a small number of factors which allow most of the variation in the
    /// data to be explained by giving the vales for the factor dimensions.</para>
    /// <para>Mathematically, PCA constructs an alternative set of orthonormal basis vectors for a multivariate data set. These
    /// basis vectors, called the principal components, are ordered by the total variance explained by each.</para>
    /// <para>Suppose, for example, you measure the value of different possessions possessions for a sample of people:
    /// home value, car value, furniture value, etc. You might expect that much of the variation in these numbers can
    /// be explained by one underlying factor, which you might call "richness". If this is true, then a PCA analysis will
    /// show that the most principal component explains a very large faction of the total variance, and the other less
    /// principal components will explain only small fractions of the total variance.</para>
    /// <para>Note that PCA is not invariant with respect to the re-scaling of individual variables.</para>
    /// <para>Note that PCA is an exploratory technique, not a hypothesis test.</para>
    /// <para>To carry out a principal component analysis, call the <see cref="Multivariate.PrincipalComponentAnalysis(IReadOnlyList{IReadOnlyList{double}})" autoUpgrade="true"/> method.</para>
    /// </remarks>
    /// <seealso href="https://en.wikipedia.org/wiki/Principal_component_analysis"/>
    public sealed class PrincipalComponentAnalysis {

        internal PrincipalComponentAnalysis (double[] utStore, double[] wStore, double[] vStore, int rows, int cols) {
            Debug.Assert(utStore != null);
            Debug.Assert(wStore != null);
            Debug.Assert(vStore != null);
            Debug.Assert(rows > 0);
            Debug.Assert(cols > 0);
            this.rows = rows;
            this.cols = cols;
            this.utStore = utStore;
            this.wStore = wStore;
            this.vStore = vStore;
            // keep track of cumulative sum of squares, which is proportional to the cumulative explained variance
            wSquaredSum = new double[wStore.Length];
            wSquaredSum[0] = MoreMath.Sqr(wStore[0]);
            for (int i = 1; i < wSquaredSum.Length; i++) {
                wSquaredSum[i] = wSquaredSum[i - 1] + MoreMath.Sqr(wStore[i]);
            }
        }

        internal readonly int rows, cols;
        internal readonly double[] utStore, wStore, vStore;
        internal readonly double[] wSquaredSum;

        /// <summary>
        /// Gets the number of components.
        /// </summary>
        /// <remarks>
        /// <para>The number of components is equal to the number of variables used in the analysis.</para>
        /// </remarks>
        public int Dimension {
            get {
                return(cols);
            }
        }

        /// <summary>
        /// Gets the number of data entries.
        /// </summary>
        public int Count {
            get {
                return(rows);
            }
        }

        /// <summary>
        /// Gets a collection of the principal components.
        /// </summary>
        public PrincipalComponentCollection Components {
            get {
                return (new PrincipalComponentCollection(this));
            }
        }

        /// <summary>
        /// Gets the minimum number of principal components that must be included to explain the given fraction of the total variance.
        /// </summary>
        /// <param name="P">The fraction of the variance to explain, which must lie between zero and one.</param>
        /// <returns>The required number of components.</returns>
        public int MinimumDimension (double P) {
            if ((P < 0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            // binary search would be faster, but linear search will do for now
            double wss = P * wSquaredSum[wSquaredSum.Length-1];
            for (int i = 0; i < wSquaredSum.Length - 1; i++) {
                if (wSquaredSum[i] >= wss) return (i + 1);
            }
            return (wSquaredSum.Length);
        }

        /// <summary>
        /// Represents the original data in terms of principal components.
        /// </summary>
        /// <returns>A multivariate sample whose columns are the weights of each principal component in each entry of the
        /// originally analyzed sample.</returns>
        public MultivariateSample TransformedSample () {
            double[] entry = new double[Dimension];
            MultivariateSample scores = new MultivariateSample(Dimension);
            for (int i = 0; i < rows; i++) {
                Array.Copy(utStore, rows * i, entry, 0, entry.Length);
                scores.Add(entry);
            }
            return (scores);
        }

    }



}
