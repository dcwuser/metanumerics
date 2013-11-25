using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Meta.Numerics {

    /// <summary>
    /// Represents a two-dimensional point.
    /// </summary>
    public struct XY {

        /// <summary>
        /// Initializes a new point with the given coordinates.
        /// </summary>
        /// <param name="x">The X-corrdinate.</param>
        /// <param name="y">The Y-coordinate.</param>
        public XY (double x, double y) {
            this.x = x;
            this.y = y;
        }

        private double x, y;

        /// <summary>
        /// Gets the X-coordinate of the point.
        /// </summary>
        public double X {
            get {
                return (x);
            }
        }

        /// <summary>
        /// Gets the Y-coordinate of the point.
        /// </summary>
        public double Y {
            get {
                return (y);
            }
        }

    }

}
