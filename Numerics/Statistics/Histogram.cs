using System;
using System.Collections.Generic;

using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// A histogram.
    /// </summary>
    public class Histogram {

        /// <summary>
        /// Creates a histogram with the given number of bins.
        /// </summary>
        /// <param name="bins">The number of bins.</param>
        public Histogram (int bins) {
            if (bins < 1) throw new ArgumentOutOfRangeException("bins");
            counts = new int[bins];
            borders = new double[bins+1];
            for (int i = 0; i < borders.Length; i++) {
                borders[i] = i;
            }
        }

        /// <summary>
        /// Create a histogram with the given bin borders.
        /// </summary>
        /// <param name="borders">A list of bin borders.</param>
        public Histogram (IList<double> borders) {
            if (borders == null) throw new ArgumentNullException("borders");
            if (borders.Count < 2) throw new InvalidOperationException();
            this.counts = new int[borders.Count - 1];
            this.borders = new double[borders.Count];
            borders.CopyTo(this.borders, 0);
            if (!CheckBorders()) throw new InvalidOperationException();
        }

        internal int[] counts;
        internal double[] borders;
        private int extraCounts;
        internal int total;

        private bool CheckBorders () {
            for (int i = 1; i < borders.Length; i++) {
                if (borders[i-1] >= borders[i]) return(false);
            }
            return(true);
        }

        // FindIndex is a simple binary search; we re-implement it here because Array.BinarySearch requires exact matching

        private int FindIndex (double value) {
            int min = 0;
            int max = counts.Length;
            while (max > min) {
                int mid = (min + max) / 2;
                if (value < borders[mid]) {
                    max = mid;
                } else if (value >= borders[mid + 1]) {
                    min = mid + 1;
                } else {
                    return (mid);
                }
            }
            return (-1);
        }

        /// <summary>
        /// Adds a new value to the histogram.
        /// </summary>
        /// <param name="value">The value to add.</param>
        public void Add (double value) {
            int index = FindIndex(value);
            if (index < 0) {
                extraCounts++;
                total++;
            } else {
                counts[index]++;
                total++;
            }
        }

        public void Add (int k) {
            if ((k < 0) || (k >= counts.Length)) {
                extraCounts++;
            } else {
                counts[k]++;
            }
            total++;
        }

        /// <summary>
        /// Sets all bin counts to zero.
        /// </summary>
        public void Clear () {
            counts = new int[counts.Length];
            extraCounts = 0;
            total = 0;
        }

        /// <summary>
        /// Gets the specified bin.
        /// </summary>
        /// <param name="k">The index of the bin.</param>
        /// <returns>The specified bin.</returns>
        public HistogramBin this[int k] {
            get {
                if ((k < 0) || (k >= counts.Length)) throw new ArgumentOutOfRangeException("k");
                return (new HistogramBin(this, k));
            }
        }

        /// <summary>
        /// Gets the total number of counts in the histogram.
        /// </summary>
        public int TotalCounts {
            get {
                return (total);
            }
        }

        /// <summary>
        /// Gets the number of counts out of the range of the histogram.
        /// </summary>
        public int OutOfRangeCounts {
            get {
                return (extraCounts);
            }
        }

        public int Count {
            get {
                return(counts.Length);
            }
        }

        /// <summary>
        /// The range of values in the histogram.
        /// </summary>
        public Interval Range {
            get {
                return (Interval.FromEndpoints(borders[0], borders[borders.Length - 1]));
            }
        }

        public TestResult ChiSquaredTest (DiscreteDistribution distribution) {

            double chi2 = 0.0;
            int dof = 0;

            //int lastObservedCounts = extraCounts;
            //double lastExpectedCounts = (distribution.LeftExclusiveProbability(0) + distribution.RightExclusiveProbability(counts.Length - 1)) * total;
            int lastObservedCounts = 0;
            double lastExpectedCounts = 0.0;

            int observedCounts = 0;
            double expectedCounts = 0.0;

            for (int i = 0; i < counts.Length; i++) {
                observedCounts += counts[i];
                expectedCounts += distribution.ProbabilityMass(i) * total;
                if (expectedCounts > 4.0) {

                    if (lastExpectedCounts > 0.0) {
                        chi2 += MoreMath.Sqr(lastObservedCounts - lastExpectedCounts) / lastExpectedCounts;
                        dof++;
                    }

                    lastObservedCounts = observedCounts;
                    lastExpectedCounts = expectedCounts;

                    observedCounts = 0;
                    expectedCounts = 0.0;
                }
            }

            lastObservedCounts += observedCounts;
            lastExpectedCounts += expectedCounts;
            if (lastExpectedCounts > 0.0) {
                chi2 += MoreMath.Sqr(lastObservedCounts - lastExpectedCounts) / lastExpectedCounts;
            }

            Distribution nullDistribution = new ChiSquaredDistribution(dof);
            return (new TestResult(chi2, nullDistribution));

        }

        /*
       public TestResult ChiSquaredTest (DiscreteDistribution distribution) {

           int N = TotalCounts;

           double chi2 = 0.0;
           for (int i = 0; i < counts.Length; i++) {
               HistogramBin bin = this[i];
               double p = 0.0;
               foreach (int k in bin.Range.GetContainedIntegers()) {
                   p += distribution.ProbabilityMass(k);
               }
               double n = p * N;
               if (n > 0.0) {
                   chi2 += MoreMath.Pow2(bin.Counts - n) / n;
               } else {
                   if (bin.Counts != 0) chi2 += Double.PositiveInfinity;
               }
           }

           return (new TestResult(chi2, new ChiSquaredDistribution(counts.Length - 1)));

        }
        */

    }
    
    /// <summary>
    /// One bin of a histogram.
    /// </summary>
    public struct HistogramBin {

        internal HistogramBin (Histogram histogram, int index) {
            this.histogram = histogram;
            this.index = index;
        }

        private Histogram histogram;
        private int index;

        /// <summary>
        /// Gets the number of counts in the bin.
        /// </summary>
        public int Counts {
            get {
                return (histogram.counts[index]);
            }
        }

        /// <summary>
        /// Gets the range of values stored in the bin.
        /// </summary>
        public Interval Range {
            get {
                return (Interval.FromEndpoints(histogram.borders[index], histogram.borders[index + 1]));
            }
        }

        /// <summary>
        /// Increments (increases by one) the bin count.
        /// </summary>
        public void Increment () {
            histogram.counts[index]++;
            histogram.total++;
        }

        /// <summary>
        /// De-increments (reduces by one) the bin count.
        /// </summary>
        /// <exception cref="InvalidOperationException">The bin count was already zero.</exception>
        public void Deincrement () {
            if (histogram.counts[index] <= 0) {
                throw new InvalidOperationException();
            } else {
                histogram.counts[index]--;
                histogram.total--;
            }
        }

    }

}
