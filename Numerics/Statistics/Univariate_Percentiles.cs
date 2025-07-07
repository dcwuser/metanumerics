﻿using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Statistics {

    public static partial class Univariate {

        /// <summary>
        /// Finds the minimum value.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>The minimum value in the sample.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        public static double Minimum (this IReadOnlyCollection<double> sample) {
            if (sample is null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 1) return Double.NaN;
            double minimum = Double.MaxValue;
            foreach (double value in sample) {
                if (value < minimum) minimum = value;
            }
            return minimum;
        }

        /// <summary>
        /// Finds the maximum value.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>The maximum value in the sample.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        public static double Maximum (this IReadOnlyCollection<double> sample) {
            if (sample is null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 1) return Double.NaN;
            double maximum = Double.MinValue;
            foreach (double value in sample) {
                if (value > maximum) maximum = value;
            }
            return maximum;
        }

        /// <summary>
        /// Finds the median.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>The median of the sample.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Median"/> 
        public static double Median (this IReadOnlyList<double> sample) {
            if (sample is null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 1) return Double.NaN;

            // Change this to use quickselect with deterministic pivots
            int[] order = GetSortOrder(sample);

            int m = sample.Count / 2;
            if (sample.Count % 2 == 0) {
                // Even n; median is average of middle two entries
                return 0.5 * (sample[order[m - 1]] + sample[order[m]]);
            } else {
                // Odd n; median is middle entry
                return sample[order[m]];
            }

        }

        /// <summary>
        /// Finds the interquartile range.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>The interquartile range of the sample.</returns>
        /// <remarks>
        /// <para>The interquartile range is the interval between the 25th and the 75th percentile.
        /// Statistics textbooks define this range as the width of this interval, but since we
        /// have to compute the endpoints to get the width, we return the full <see cref="Interval"/>.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <seealso cref="InverseLeftProbability(IReadOnlyList{Double},Double)"/>
        /// <seealso href="https://en.wikipedia.org/wiki/Interquartile_range"/>
        public static Interval InterquartileRange (this IReadOnlyList<double> sample) {
            if (sample is null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 2) throw new InsufficientDataException();
            int[] order = GetSortOrder(sample);
            return Interval.FromEndpoints(InverseLeftProbability(sample, 0.25, order), InverseLeftProbability(sample, 0.75, order));
        }

        /// <summary>
        /// Finds the tri-mean.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>Tukey's tri-mean for the sample.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Trimean"/>
        public static double Trimean (this IReadOnlyList<double> sample) {
            if (sample is null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 1) return Double.NaN;
            int[] order = GetSortOrder(sample);
            double x1 = InverseLeftProbability(sample, 0.25, order);
            double x2 = InverseLeftProbability(sample, 0.50, order);
            double x3 = InverseLeftProbability(sample, 0.75, order);
            return 0.5 * x2 + 0.25 * (x1 + x3);
        }

        /// <summary>
        /// Finds the sample value corresponding to a given percentile.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <param name="P">The percentile, which must lie between zero and one.</param>
        /// <returns>The corresponding value.</returns>
        /// <remarks>
        /// <para>For percentile inputs which do not lie exactly on a sample value,
        /// this method interpolates linearly between the two nearest sample values.
        /// Therefore the value returned is not guaranteed to be an acutal value in the sample.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="P"/> lies outside [0,1].</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Quantile_function"/>
        public static double InverseLeftProbability (this IReadOnlyList<double> sample, double P) {
            if (sample is null) throw new ArgumentNullException(nameof(sample));
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            if (sample.Count < 2) throw new InsufficientDataException();
            int[] order = GetSortOrder(sample);
            return InverseLeftProbability(sample, P, order);
        }

        private static double InverseLeftProbability(IReadOnlyList<double> sample, double P, IReadOnlyList<int> order) {
            Debug.Assert(sample != null);
            Debug.Assert((0.0 <= P) && (P <= 1.0));
            Debug.Assert(order != null);
            Debug.Assert(order.Count == sample.Count);
            double n = P * (sample.Count - 1);
            int n1 = (int) Math.Floor(n);
            int n2 = (int) Math.Ceiling(n);
            double w1 = n2 - n;
            double w2 = 1.0 - w1;
            return w1 * sample[order[n1]] + w2 * sample[order[n2]];
        }

        /// <summary>
        /// Gets the fraction of values equal to or less than the given value.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <param name="value">The reference value.</param>
        /// <returns>The empirical distribution function.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <see href="https://en.wikipedia.org/wiki/Empirical_distribution_function"/>
        public static double LeftProbability (this IReadOnlyList<double> sample, double value) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            int[] order = GetSortOrder(sample);
            Debug.Assert(order.Length == sample.Count);
            for (int i = 0; i < order.Length; i++) {
                if (sample[order[i]] > value) {
                    double P = ((double) i) / order.Length;
                    return (P);
                }
            }
            return (1.0);
        }

        private static int[] GetSortOrder (IReadOnlyList<double> sample) {
            Debug.Assert(sample != null);
            int[] order = new int[sample.Count];
            for (int i = 0; i < sample.Count; i++) {
                order[i] = i;
            }
            Array.Sort(order, (int i, int j) => sample[i].CompareTo(sample[j]));
            return (order);
        }

        internal static int[] GetRanks (IReadOnlyList<double> sample) {
            Debug.Assert(sample != null);
            int[] order = GetSortOrder(sample);
            int[] ranks = new int[order.Length];
            for (int i = 0; i < order.Length; i++) {
                ranks[order[i]] = i;
            }
            return (ranks);
            // Is there a way we can do this without two array allocations?
        }

    }

}
