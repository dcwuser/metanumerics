using System;
using System.Collections.Generic;


namespace Meta.Numerics.Statistics {

    public class Histogram {

        private int[] counts;

        public int this[int k] {
            get {
                if ((k < 0) || (k >= counts.Length)) throw new ArgumentOutOfRangeException("k");
                return (counts[k]);
            }
            set {
                if ((k < 0) || (k >= counts.Length)) throw new ArgumentOutOfRangeException("k");
                if (value < 0) throw new InvalidOperationException();
                counts[k] = value;
            }
        }


    }

}
