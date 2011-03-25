using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a principal component analysis.
    /// </summary>
    /// <remarks>
    /// <para>Principal component analysis constructs an alternative set of orthonormal basis vectors for a multi-variate data set. These  represents multi-variate data set (<see cref="MultivariateSample"/>) 
    /// It is a form of factor analysis, ...</para>
    /// <para>Suppose, for example, that most of the data in a two-dimensional data set lies near the line y = -x. One way to explain
    /// this state of affairs is to invoke a single underlying factor. The factor increases the value of y and decreases the value
    /// of x. Other factors values of this
    /// factor increases the value of y and decreases the value of x, and this is the most important factor  there is one underlying factor, which mostly determines the values of x and y. Increasing
    /// values of this factor increase </para>
    /// </remarks>
    public class PrincipalComponentAnalysis {

        internal PrincipalComponentAnalysis (double[] utStore, double[] wStore, double[] vStore, int rows, int cols) {
            this.rows = rows;
            this.cols = cols;
            this.utStore = utStore;
            this.wStore = wStore;
            this.vStore = vStore;
            // keep track of cumulative sum of squares, which is proportional to the cumulative explained variance
            wSquaredSum = new double[wStore.Length];
            wSquaredSum[0] = MoreMath.Pow2(wStore[0]);
            for (int i = 1; i < wSquaredSum.Length; i++) {
                wSquaredSum[i] = wSquaredSum[i - 1] + MoreMath.Pow2(wStore[i]);
            }
        }

        internal int rows, cols;
        internal double[] utStore, wStore, vStore;
        internal double[] wSquaredSum;

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
        /// Gets the requested principal component.
        /// </summary>
        /// <param name="componentIndex">The (zero-based) index of the principal component.</param>
        /// <returns>The requested principal component.</returns>
        /// <remarks>
        /// <para>Principal components are ordered by strength. The most principal component, i.e. the component which explains
        /// the most variance, has index zero. The least principal component has the highest index.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="componentIndex"/> lies outside the range [0, <see cref="Dimension"/>-1].</exception>
        public PrincipalComponent Component (int componentIndex) {
            if ((componentIndex < 0) || (componentIndex >= Dimension)) throw new ArgumentOutOfRangeException("componentIndex");
            return(new PrincipalComponent(componentIndex, this));
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

    /// <summary>
    /// Represents a component of a principal component analysis.
    /// </summary>
    /// <seealso cref="PrincipalComponentAnalysis"/>
    public class PrincipalComponent {

        internal PrincipalComponent (int index, PrincipalComponentAnalysis analysis) {
            this.index = index;
            this.analysis = analysis;
        }

        private int index;
        private PrincipalComponentAnalysis analysis;

        /// <summary>
        /// Gets the index of the principal component.
        /// </summary>
        public int Index { 
            get {
                return(index);
            }
        }

        /// <summary>
        /// Gets the principal component analysis to which the component belongs.
        /// </summary>
        public PrincipalComponentAnalysis Analysis {
            get {
                return(analysis);
            }
        }

        /// <summary>
        /// Gets the weight of the component.
        /// </summary>
        public double Weight {
            get {
                return (analysis.wStore[index]);
            }
        }

        /// <summary>
        /// Gets the fraction of the total variance accounted for by the principal component.
        /// </summary>
        public double VarianceFraction {
            get {
                return (MoreMath.Pow2(analysis.wStore[index]) / analysis.wSquaredSum[analysis.wSquaredSum.Length - 1]);
            }
        }

        /// <summary>
        /// Gets the fraction of the total variance accounted for by the principal component and all strong principal components.
        /// </summary>
        public double CumulativeVarianceFraction {
            get {
                return (analysis.wSquaredSum[index] / analysis.wSquaredSum[analysis.wSquaredSum.Length - 1]);
            }
        }

        /// <summary>
        /// Gets the normalized component vector.
        /// </summary>
        public RowVector NormalizedVector {
            get {
                double[] pc = new double[analysis.cols];
                Blas1.dCopy(analysis.vStore, analysis.cols * index, 1, pc, 0, 1, analysis.cols);
                return (new RowVector(pc, pc.Length));
            }
        }

        /// <summary>
        /// Gets the scaled component vector.
        /// </summary>
        public RowVector ScaledVector {
            get {
                return (analysis.wStore[index] * NormalizedVector);
            }
        }

    }

}
