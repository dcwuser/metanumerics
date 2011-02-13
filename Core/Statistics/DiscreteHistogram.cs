using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    public class DiscreteHistogram {

        public DiscreteHistogram (DiscreteInterval interval) {
            this.interval = interval;
            this.counts = new int[interval.Width + 1];
        }

        private DiscreteInterval interval;
        private int[] counts;
        private int extra;
        private int total;

        public void Add (int k) {
            if ((k < interval.LeftEndpoint) || (k > interval.RightEndpoint)) {
                extra++;
            } else {
                counts[interval.LeftEndpoint - k]++;
            }
            total++;
        }



    }

}
