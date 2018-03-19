using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {
    /// <summary>
    /// Represents the result of a k-means clustering analysis.
    /// </summary>
    /// <seealso cref="MultivariateSample.MeansClustering(int)"/>
    public sealed class MeansClusteringResult {

        internal MeansClusteringResult (double[][] centroids) {
            Debug.Assert(centroids != null);
            Debug.Assert(centroids.Length > 0);
            this.centroids = centroids;
        }

        private readonly double[][] centroids;

        /// <summary>
        /// Gets the dimension of the space.
        /// </summary>
        public int Dimension {
            get {
                return (centroids[0].Length);
            }
        }

        /// <summary>
        /// Gets the number of clusters.
        /// </summary>
        public int Count {
            get {
                return (centroids.Length);
            }
        }

        /// <summary>
        /// Assign a vector to a cluster.
        /// </summary>
        /// <param name="values">The vector to classify.</param>
        /// <returns>The index of the cluster to which the vector belongs.</returns>
        public int Classify (IReadOnlyList<double> values) {

            if (values == null) throw new ArgumentNullException(nameof(values));
            if (values.Count != this.Dimension) throw new DimensionMismatchException();

            double minDistance = Double.MaxValue;
            int minK = -1;
            for (int k = 0; k < centroids.Length; k++) {
                double distance = 0.0;
                for (int j = 0; j < values.Count; j++) {
                    distance += MoreMath.Sqr(values[j] - centroids[k][j]);
                }
                if (distance < minDistance) {
                    minDistance = distance;
                    minK = k;
                }
            }

            return (minK);
        }

        /// <summary>
        /// Get the centroid of the given cluster.
        /// </summary>
        /// <param name="k">The index of the cluster, which must lie between zero and the cluster count.</param>
        /// <returns>The centroid of the given index.</returns>
        public ColumnVector Centroid (int k) {
            if ((k < 0) || (k >= centroids.Length)) throw new ArgumentOutOfRangeException(nameof(k));
            return (new ColumnVector(centroids[k], 0, 1, centroids[k].Length, true));
        }

    }

}
