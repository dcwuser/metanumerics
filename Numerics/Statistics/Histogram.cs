using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    // The internal HistorgramStorage stores all actual historgram data.
    // The Histogram, HistogramBinCollection, and HistogramBin objects associated with a histogram all just hold references to the
    // storage object and surface methods from it appropriate to their roles.

    // To accomodate BelowRangeBin and AboveRangeBin, the internal index used to indicate a bin is one off from its external index.
    // Thus bins 0..(Count - 1) have indexes 1..(Count); index 0 represents the BelowRangeBin and index (Count + 1) represents the
    // AboveRangeBin.

    internal class HistogramStorage {

        public HistogramStorage (double[] borders) {
            this.counts = new int[borders.Length + 1];
            this.borders = borders;
            this.total = 0;
        }

        private int[] counts;

        private double[] borders;

        private int total;

        public int Count {
            get {
                return (borders.Length - 1);
            }
        }

        public int Total {
            get {
                return (total);
            }
        }

        public int GetCounts (int index) {
            return (counts[index]);
        }

        public void SetCounts (int index, int value) {
            counts[index] = value;
        }

        public Interval GetBorders (int index) {
            if (index == 0) {
                return (Interval.FromEndpoints(Double.NegativeInfinity, borders[0]));
            } else if (index < borders.Length) {
                return (Interval.FromEndpoints(borders[index - 1], borders[index]));
            } else {
                return (Interval.FromEndpoints(borders[borders.Length - 1], Double.PositiveInfinity));
            }
        }

        public Interval Range {
            get {
                return (Interval.FromEndpoints(borders[0], borders[borders.Length - 1]));
            }
        }

        public void Increment (int index) {
            counts[index]++;
            total++;
        }

        public void Decriment (int index) {
            if (counts[index] > 0) {
                counts[index]--;
                total--; 
            } else {
                throw new InvalidOperationException();
            }
        }

        public void Clear () {
            counts = new int[counts.Length];
            total = 0;
        }


        public int FindIndex (double value) {
            int min = -1;
            int max = borders.Length;
            while (max > min + 1) {
                int mid = (min + max) / 2;
                if (value < borders[mid]) {
                    max = mid;
                } else {
                    min = mid;
                }
            }
            return (max);
        }

    }

    /// <summary>
    /// Represents a histogram.
    /// </summary>
    /// <remarks>
    /// <para>A histogram stores the number of occurances of a value (called the count) within a set of ranges (called bins). It is the
    /// natural system to store discrete univariate data. The natural system to store continuous univariate data is the <see cref="Sample"/>
    /// class, but often histograms are also used for continuous data, for example because there is too much data to record each value
    /// individually.</para>
    /// </remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/Histogram"/>
    public sealed class Histogram {

        /// <summary>
        /// Creates a histogram with the given number of bins.
        /// </summary>
        /// <param name="binCount">The number of bins.</param>
        /// <remarks>
        /// <para>By default, the bin borders are the integers, i.e. values in [0, 1) fall into bin 0, values in [1, 2) fall
        /// into bin 1, etc.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="binCount"/> is less than 1.</exception>
        public Histogram (int binCount) {
            if (binCount < 1) throw new ArgumentOutOfRangeException("binCount");
            double[] borders = new double[binCount + 1];
            for (int i = 0; i < borders.Length; i++) {
                borders[i] = i;
            }
            storage = new HistogramStorage(borders);
            bins = new HistogramBinsCollection(storage);
        }

        /// <summary>
        /// Create a histogram with the given bin borders.
        /// </summary>
        /// <param name="binBorders">A list of bin borders.</param>
        /// <exception cref="ArgumentNullException"><paramref name="binBorders"/> is null.</exception>
        /// <exception cref="ArgumentException"><paramref name="binBorders"/> has less than two values, or the values are not ordered.</exception>
        public Histogram (IList<double> binBorders) {
            if (binBorders == null) throw new ArgumentNullException("binBorders");
            if (binBorders.Count < 2) throw new ArgumentException("The bin borders list must be at least two elements long.", "binBorders");
            if (!ValidateBorderMononicity(binBorders)) throw new ArgumentException("The values in the bin borders list must increase strictly monotonically.", "binBorders");
            double[] borders = new double[binBorders.Count];
            binBorders.CopyTo(borders, 0);
            storage = new HistogramStorage(borders);
            bins = new HistogramBinsCollection(storage);
        }

        private readonly HistogramStorage storage;

        // We construct the HistogramBinsCollection object once when we construct the Histogram object
        // instead of creating a new one every time the Bins property is read.
        private readonly HistogramBinsCollection bins;

        private static bool ValidateBorderMononicity (IList<double> borders) {
            Debug.Assert(borders != null);
            for (int i = 1; i < borders.Count; i++) {
                if (borders[i-1] >= borders[i]) return(false);
            } 
            return(true);
        }

        // FindIndex is a simple binary search; we re-implement it here because Array.BinarySearch requires exact matching
        // A bin runs from [a,b), i.e. given value a, the count should fall into that bin; given value b, the count should
        // fall into the next bin.

        /// <summary>
        /// Adds a new value to the histogram.
        /// </summary>
        /// <param name="value">The value to add.</param>
        public void Add (double value) {
            storage.Increment(storage.FindIndex(value));
        }

        /// <summary>
        /// Sets all bin counts to zero.
        /// </summary>
        public void Clear () {
            storage.Clear();
        }

        /// <summary>
        /// Gets a collection of histogram bins.
        /// </summary>
        public HistogramBinsCollection Bins {
            get {
                return (bins);
            }
        }

        /// <summary>
        /// Gets the total number of counts in the histogram.
        /// </summary>
        public int TotalCounts {
            get {
                return (storage.Total);
            }
        }

        /// <summary>
        /// Gets a bin for all values below the histogram range.
        /// </summary>
        public HistogramBin BelowRangeBin {
            get {
                return (new HistogramBin(storage, 0));
            }
        }

        /// <summary>
        /// Gets a bin for all values above the histogram range.
        /// </summary>
        public HistogramBin AboveRangeBin {
            get {
                return (new HistogramBin(storage, storage.Count + 1));
            }
        }


        /// <summary>
        /// The range of values in the histogram.
        /// </summary>
        public Interval Range {
            get {
                return (storage.Range);
            }
        }

        /// <summary>
        /// Performs a  &#x3C7;<sup>2</sup> test comparing the histogram to the given distribution.
        /// </summary>
        /// <param name="distribution">The distribution against which to test the histogram.</param>
        /// <returns>The test result.</returns>
        public TestResult ChiSquaredTest (DiscreteDistribution distribution) {

            double chi2 = 0.0;
            int dof = 0;

            //int lastObservedCounts = extraCounts;
            //double lastExpectedCounts = (distribution.LeftExclusiveProbability(0) + distribution.RightExclusiveProbability(counts.Length - 1)) * total;
            int lastObservedCounts = 0;
            double lastExpectedCounts = 0.0;

            int observedCounts = 0;
            double expectedCounts = 0.0;

            for (int i = 0; i < storage.Count; i++) {
                observedCounts += storage.GetCounts(i + 1);
                expectedCounts += distribution.ProbabilityMass(i) * storage.Total;
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
            return (new TestResult("ChiSquared", chi2, TestType.RightTailed, nullDistribution));

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
    /// Represents one bin in a histogram.
    /// </summary>
    /// <seealso cref="Histogram"/>
    public struct HistogramBin {

        internal HistogramBin (HistogramStorage storage, int index) {
            this.storage = storage;
            this.index = index;
        }

        internal readonly HistogramStorage storage;
        private readonly int index;

        /// <summary>
        /// Gets or sets the number of counts in the bin.
        /// </summary>
        public int Counts {
            get {
                return (storage.GetCounts(index));
            }
            set {
                if (value < 0) throw new ArgumentOutOfRangeException("value");
                storage.SetCounts(index, value);
            }
        }

        /// <summary>
        /// Gets the range of values stored in the bin.
        /// </summary>
        public Interval Range {
            get {
                return (storage.GetBorders(index));
            }
        }

        /// <summary>
        /// Increments (increases by one) the bin count.
        /// </summary>
        public void Increment () {
            storage.Increment(index);
        }

        /// <summary>
        /// De-increments (reduces by one) the bin count.
        /// </summary>
        /// <exception cref="InvalidOperationException">The bin count was already zero.</exception>
        public void Deincrement () {
            storage.Decriment(index);
        }

    }

    /// <summary>
    /// Represents a collection of histogram bins.
    /// </summary>
    /// <seealso cref="Histogram"/>
    public sealed class HistogramBinsCollection : ICollection<HistogramBin>, IEnumerable<HistogramBin> {

        internal HistogramBinsCollection (HistogramStorage storage) {
            this.storage = storage;
        }

        private readonly HistogramStorage storage;

        /// <summary>
        /// Gets the given histogram bin.
        /// </summary>
        /// <param name="index">The bin number.</param>
        /// <returns>The specified bin.</returns>
        public HistogramBin this[int index] {
            get {
                if ((index < 0) || (index >= storage.Count)) throw new ArgumentOutOfRangeException("index");
                return (new HistogramBin(storage, index + 1));
            }
        }

        void ICollection<HistogramBin>.Add (HistogramBin item) {
            throw new InvalidOperationException();
        }

        void ICollection<HistogramBin>.Clear () {
            throw new InvalidOperationException();
        }

        bool ICollection<HistogramBin>.Contains (HistogramBin item) {
            return (item.storage == this.storage);
        }

        void ICollection<HistogramBin>.CopyTo (HistogramBin[] array, int arrayIndex) {
            for (int index = 0; index < storage.Count; index++) {
                array[arrayIndex + index] = new HistogramBin(storage, index + 1);
            }
        }

        /// <summary>
        /// Gets the number of histogram bins in the collection.
        /// </summary>
        public int Count {
            get {
                return(storage.Count);
            }
        }

        bool ICollection<HistogramBin>.IsReadOnly {
            get {
                return(true);
            }
        }

        bool ICollection<HistogramBin>.Remove (HistogramBin item) {
            throw new InvalidOperationException();
        }

        IEnumerator<HistogramBin> IEnumerable<HistogramBin>.GetEnumerator () {
            for (int index = 0; index < storage.Count; index++) {
                yield return (new HistogramBin(storage, index + 1));
            }
        }

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator () {
            return (((IEnumerable<HistogramBin>) this).GetEnumerator());
        }
    }

}
