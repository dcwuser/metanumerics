using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;


namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a collection of principal components.
    /// </summary>
    public sealed class PrincipalComponentCollection : IReadOnlyList<PrincipalComponent>, IReadOnlyCollection<PrincipalComponent>, IEnumerable<PrincipalComponent> {

        internal PrincipalComponentCollection(PrincipalComponentAnalysis analysis) {
            Debug.Assert(analysis != null);
            this.analysis = analysis;
        }

        private readonly PrincipalComponentAnalysis analysis;

        /// <summary>
        /// Gets the principal component with the given index.
        /// </summary>
        /// <param name="index">The (zero-based) index.</param>
        /// <returns>The principal component with the given index.</returns>
        /// <remarks>
        /// <para>Principal components are ordered from most to least significant.
        /// The most principal component, i.e. the component which explains the most variance, has index zero.
        /// The least principal component has the highest index.</para>
        /// </remarks>
        public PrincipalComponent this[int index] {
            get {
                return (new PrincipalComponent(index, analysis));
            }
        }

        /// <summary>
        /// Gets the count of principal components.
        /// </summary>
        public int Count {
            get {
                return (analysis.cols);
            }
        }

        IEnumerator<PrincipalComponent> IEnumerable<PrincipalComponent>.GetEnumerator () {
            for (int i = 0; i < this.Count; i++) {
                yield return (this[i]);
            }
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (((IEnumerable<PrincipalComponent>) this).GetEnumerator());
        }
    }
}
