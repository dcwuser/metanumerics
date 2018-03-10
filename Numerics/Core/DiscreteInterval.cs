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

        internal DiscreteInterval (int min, int max) {
            this.min = min;
            this.max = max;
        }

        private readonly int min;

        private readonly int max;

        /// <summary>
        /// Instantiates a new discrete interval from the given endpoints.
        /// </summary>
        /// <param name="a">One endpoint.</param>
        /// <param name="b">The other endpoint.</param>
        /// <returns>A discrete interval between the given endpoints.</returns>
        public static DiscreteInterval FromEndpoints (int a, int b) {
            return (new DiscreteInterval(Math.Min(a, b), Math.Max(a, b)));
        }

        /// <summary>
        /// Gets the minimum value on the interval.
        /// </summary>
        public int LeftEndpoint {
            get {
                return (min);
            }
        }

        /// <summary>
        /// Gets the maximum value on the interval.
        /// </summary>
        public int RightEndpoint {
            get {
                return (max);
            }
        }

        /// <summary>
        /// Gets the width of the interval.
        /// </summary>
        internal int Width {
            get {
                // For [0, Int32.MaxValue], w will overflow. Since Int32.MaxValue is our integer "infinity",
                // just reset it to that value.
                int w = max - min + 1;
                if (w < 0) w = Int32.MaxValue;
                return (w);
            }
        }

    }
}
