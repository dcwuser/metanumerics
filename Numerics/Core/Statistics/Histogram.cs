using System;
using System.Collections.Generic;


namespace Meta.Numerics.Statistics {

    public class Histogram {

        public Histogram (int bins) {
            if (bins < 1) throw new ArgumentOutOfRangeException("bins");
            counts = new int[bins];
            borders = new double[bins+1];
            for (int i = 0; i < borders.Length; i++) {
                borders[i] = i;
            }
        }

        public Histogram (double[] borders) {
            if (borders == null) throw new ArgumentNullException("borders");
            if (borders.Length < 2) throw new InvalidOperationException();
            this.counts = new int[borders.Length - 1];
            this.borders = borders;
            if (!CheckBorders()) throw new InvalidOperationException();
        }

        internal int[] counts;
        internal double[] borders;

        private bool CheckBorders () {
            for (int i = 1; i < borders.Length; i++) {
                if (borders[i-1] >= borders[i]) return(false);
            }
            return(true);
        }

        // FindIndex is a simple binary search; we re-implement it here because Array.BinarySearch does exact matching

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

        public void Add (double value) {
            int index = FindIndex(value);
            if (index < 0) {
                Console.WriteLine(value);
                throw new ArgumentOutOfRangeException("value");
            } else {
                counts[index]++;
            }
        }

        public int this[int k] {
            get {
                if ((k < 0) || (k >= counts.Length)) throw new ArgumentOutOfRangeException("k");
                return (counts[k]);
            }
            set {
                if (k < 0) throw new ArgumentOutOfRangeException("k");
                if (value < 0) throw new InvalidOperationException();
                counts[k] = value;
            }
        }

        public int Counts {
            get {
                int total = 0;
                for (int i = 0; i < counts.Length; i++) {
                    total += counts[i];
                }
                return (total);
            }
        }

        public Interval Range {
            get {
                return (Interval.FromEndpoints(borders[0], borders[borders.Length - 1]));
            }
        }


    }

    public struct HistogramBin {

        internal HistogramBin (Histogram histogram, int index) {
            this.histogram = histogram;
            this.index = index;
        }

        private Histogram histogram;
        private int index;

        public int Counts {
            get {
                return (histogram.counts[index]);
            }
        }

        public Interval Range {
            get {
                return (Interval.FromEndpoints(histogram.borders[index], histogram.borders[index + 1]));
            }
        }

        public void Increment () {
            histogram.counts[index]++;
        }

    }

}
