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
    /// <para>Principal component analysis represents multi-variate data set (<see cref="MultivariateSample"/>)...
    /// It is a form of factor analysis, ...</para>
    /// </remarks>
    public class PrincipalComponentAnalysis {

        internal PrincipalComponentAnalysis (double[] uStore, double[] wStore, double[] vStore, int rows, int cols) {
            this.rows = rows;
            this.cols = cols;
            this.uStore = uStore;
            this.wStore = wStore;
            this.vStore = vStore;
            wSquaredSum = new double[wStore.Length];
            wSquaredSum[0] = MoreMath.Pow2(wStore[0]);
            for (int i = 1; i < wSquaredSum.Length; i++) {
                wSquaredSum[i] = wSquaredSum[i - 1] + MoreMath.Pow2(wStore[i]);
            }
        }

        internal int rows, cols;
        internal double[] uStore, wStore, vStore;
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
            throw new NotImplementedException();
        }

        /// <summary>
        /// Gets the requested principal component.
        /// </summary>
        /// <param name="componentIndex">The index of the requested principal component.</param>
        /// <returns>The requested principal component.</returns>
        /// <remarks>
        /// <para>Principal components are ordered by strength. The strongest component, i.e. the component which explains
        /// the most variance, has index zero. The weakest component has the highest index.</para>
        /// </remarks>
        public PrincipalComponent Component (int componentIndex) {
            if ((componentIndex < 0) || (componentIndex >= Dimension)) throw new ArgumentOutOfRangeException("componentIndex");
            return(new PrincipalComponent(componentIndex, this));
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
        public RowVector Vector {
            get {
                double[] pc = new double[analysis.cols];
                Blas1.dCopy(analysis.vStore, analysis.cols * index, 1, pc, 0, 1, analysis.cols);
                return (new RowVector(pc, pc.Length));
            }
        }

    }

}
