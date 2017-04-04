using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics {

    /// <summary>
    /// Represents an interval on the integers.
    /// </summary>
    public struct DiscreteInterval {

        public DiscreteInterval (int min, int max) {
            this.min = min;
            this.max = max;
        }

        private readonly int min;

        private readonly int max;

        /// <summary>
        /// The minimum value on the interval.
        /// </summary>
        public int LeftEndpoint {
            get {
                return (min);
            }
        }

        /// <summary>
        /// The maximum value on the interval.
        /// </summary>
        public int RightEndpoint {
            get {
                return (max);
            }
        }

    }
}
