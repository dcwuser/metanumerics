using System;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a component of a principal component analysis.
    /// </summary>
    /// <seealso cref="PrincipalComponentAnalysis"/>
    public sealed class PrincipalComponent {

        internal PrincipalComponent (int index, PrincipalComponentAnalysis analysis) {
            this.index = index;
            this.analysis = analysis;
        }

        private readonly int index;
        private readonly PrincipalComponentAnalysis analysis;

        /// <summary>
        /// Gets the index of the principal component.
        /// </summary>
        public int Index {
            get {
                return (index);
            }
        }

        /// <summary>
        /// Gets the principal component analysis to which the component belongs.
        /// </summary>
        public PrincipalComponentAnalysis Analysis {
            get {
                return (analysis);
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
                return (MoreMath.Sqr(analysis.wStore[index]) / analysis.wSquaredSum[analysis.wSquaredSum.Length - 1]);
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
        /// <value>A unit-length vector that defines the direction of the principal component.</value>
        public RowVector NormalizedVector {
            get {
                return (new RowVector(analysis.vStore, analysis.cols * index, 1, analysis.cols, true));
            }
            //get {
            //double[] pc = new double[analysis.cols];
            //Blas1.dCopy(analysis.vStore, analysis.cols * index, 1, pc, 0, 1, analysis.cols);
            //return (new RowVector(pc, pc.Length));
            //}
        }

        /// <summary>
        /// Gets the scaled component vector.
        /// </summary>
        /// <returns>The <see cref="NormalizedVector"/> multiplied by its associated weight.</returns>
        public RowVector ScaledVector () {
            //get {
            return (analysis.wStore[index] * NormalizedVector);
            //}
        }

    }

}
