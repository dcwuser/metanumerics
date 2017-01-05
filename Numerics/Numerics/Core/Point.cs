using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Meta.Numerics {

    /// <summary>
    /// Represents a two-dimensional point.
    /// </summary>
    public struct XY : IEquatable<XY> {

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

        /// <summary>
        /// Converts a point into a two-tuple.
        /// </summary>
        /// <param name="point">The point to convert.</param>
        /// <returns>The equivilent two-tuple, with <see cref="X"/> as the first item and
        /// <see cref="Y"/> as the second item.</returns>
        public static implicit operator Tuple<double, double> (XY point) {
            return (new Tuple<double, double>(point.x, point.y));
        }

        /// <summary>
        /// Converts a two-tuple into a point.
        /// </summary>
        /// <param name="point">The two-tuple to convert.</param>
        /// <returns>The equivilent point, with <see cref="X"/> equalto first item and
        /// <see cref="Y"/> equal to the second item.</returns>
        public static implicit operator XY (Tuple<double, double> point) {
            return (new XY(point.Item1, point.Item2));
        }

        private static bool Equals (XY a, XY b) {
            return ((a.x == b.x) && (a.y == b.y));
        }

        /// <summary>
        /// Determines whether two points are equal.
        /// </summary>
        /// <param name="a">The first point.</param>
        /// <param name="b">The second point.</param>
        /// <returns><see langword="true"/> if <paramref name="a"/> and <paramref name="b"/> refer to the same
        /// point, otherwise <see langword="false"/>.</returns>
        public static bool operator == (XY a, XY b) {
            return (Equals(a, b));
        }

        /// <summary>
        /// Determines whether two points are unequal.
        /// </summary>
        /// <param name="a">The first point.</param>
        /// <param name="b">The second point.</param>
        /// <returns><see langword="true"/> if <paramref name="a"/> and <paramref name="b"/> refer to different
        /// points, otherwise <see langword="false"/>.</returns>
        public static bool operator != (XY a, XY b) {
            return (!Equals(a, b));
        }

        /// <summary>
        /// Determines whether the given point is the same.
        /// </summary>
        /// <param name="point">The point to compare.</param>
        /// <returns><see langword="true"/> if <paramref name="point"/> refers to the same point, otherwise <see langword="false"/>.</returns>
        public bool Equals (XY point) {
            return (Equals(this, point));
        }

        /// <inheritdoc/>
        public override bool Equals (object obj) {
            if (obj is XY) {
                XY other = (XY) obj;
                return (Equals(this, other));
            } else {
                return (false);
            }
        }

        /// <inheritdoc/>
        public override int GetHashCode () {
            return (x.GetHashCode() + 17 * y.GetHashCode());
        }

        /// <inheritdoc/>
        public override string ToString () {
            return ($"({x},{y})");
        }

    }

}
