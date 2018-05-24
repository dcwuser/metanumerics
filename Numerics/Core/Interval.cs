using System;
using System.Collections.Generic;
using System.Globalization;

namespace Meta.Numerics {

    /// <summary>
    /// Represents an interval on the real number line.
    /// </summary>
    /// <remarks>
    /// <para>Use the static methods <see cref="FromEndpoints"/>, <see cref="FromMidpointAndWidth"/>,
    /// and <see cref="FromEndpointAndWidth"/> to instantiate intervals.</para>
    /// </remarks>
    public struct Interval {

        private double a, b, w;

        private Interval (double a, double b, double w) {
            this.a = a;
            this.b = b;
            this.w = w;
        }

        /// <summary>
        /// Gets the left (lower) endpoint of the interval.
        /// </summary>
        public double LeftEndpoint {
            get {
                return (a);
            }
        }

        /// <summary>
        /// Gets the right (upper) endpoint of the interval.
        /// </summary>
        public double RightEndpoint {
            get {
                return (b);
            }
        }

        /// <summary>
        /// Determines whether the argument lies in the open interval.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>True if <paramref name="x"/> lies in (a,b), otherwise False.</returns>
        public bool OpenContains (double x) {
            return ((x > a) && (x < b));
        }
        /// <summary>
        /// Determines whether the argument lies in the closed interval.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>True if <paramref name="x"/> lies in [a,b], otherwise False.</returns>
        public bool ClosedContains (double x) {
            return ((x >= a) && (x <= b));
        }

        /// <summary>
        /// Gets the width of the interval.
        /// </summary>
        public double Width {
            get {
                return (w);
            }
        }

        /// <summary>
        /// Gets the mid-point of the interval.
        /// </summary>
        public double Midpoint {
            get {
                return ((a + b) / 2.0);
            }
        }

        /// <summary>
        /// Creates a new interval, given its endpoints.
        /// </summary>
        /// <param name="a">The left (lower) endpoint of the interval.</param>
        /// <param name="b">The right (upper) endpoint of the interval.</param>
        /// <returns>The specified interval.</returns>
        /// <remarks>If width of the interval is very much smaller than its endpoint values, accuracy will be better maintained by constructing the interval using one endpoint and its width.</remarks>
        public static Interval FromEndpoints (double a, double b) {
            if (b >= a) {
                return (new Interval(a, b, b - a));
            } else {
                return (new Interval(b, a, a - b));
            }
        }

        /// <summary>
        /// Creates a new interval, given its lower endpoint and width.
        /// </summary>
        /// <param name="endpoint">The left (lower) endpoint of the interval.</param>
        /// <param name="width">The width of the interval.</param>
        /// <returns>The specified interval.</returns>
        public static Interval FromEndpointAndWidth (double endpoint, double width) {
            if (width < 0.0) {
                return (new Interval(endpoint + width, endpoint, -width));
            } else {
                return (new Interval(endpoint, endpoint+width, width));
            }
        }


        /// <summary>
        /// Creates a new interval, given its midpoint and width.
        /// </summary>
        /// <param name="midpoint">The midpoint of the interval.</param>
        /// <param name="width">The width of the interval.</param>
        /// <returns>The specified interval.</returns>
        public static Interval FromMidpointAndWidth (double midpoint, double width) {
            if (width < 0.0) {
                return (FromMidpointAndWidth(midpoint, -width));
            } else {
                return (new Interval(midpoint - width / 2.0, midpoint + width / 2.0, width));
            }
        }

        // equality

        /// <summary>
        /// Determines whether two intervals are equal.
        /// </summary>
        /// <param name="i1">The first interval.</param>
        /// <param name="i2">The second interval.</param>
        /// <returns>True if <paramref name="i1"/> and <paramref name="i2"/> are equal, otherwise false.</returns>
        public static bool operator == (Interval i1, Interval i2) {
            // struct, so no null testing
            return ((i1.a == i2.a) && (i1.w == i2.w));
        }

        /// <summary>
        /// Determines whether two intervals are not equal.
        /// </summary>
        /// <param name="i1">The first interval.</param>
        /// <param name="i2">The second interval.</param>
        /// <returns>True if <paramref name="i1"/> and <paramref name="i2"/> are not equal, otherwise false.</returns>
        public static bool operator != (Interval i1, Interval i2) {
            return (!(i1 == i2));
        }

        /// <summary>
        /// Determines whether a given object is an equal interval.
        /// </summary>
        /// <param name="obj">An object.</param>
        /// <returns>True if <paramref name="obj"/> is an equal <see cref="Interval"/>, otherwise false.</returns>
        public override bool Equals (object obj) {
            if (obj is Interval) {
                Interval i = (Interval) obj;
                return (i == this);
            } else {
                return (false);
            }
        }

        /// <inheritdoc />
        public override int GetHashCode () {
            return (a.GetHashCode() ^ w.GetHashCode());
        }

        // text representation
        /// <summary>
        /// Produces a string representation of the interval.
        /// </summary>
        /// <returns>A string representation of the interval.</returns>
        public override string ToString () {
            return (String.Format(CultureInfo.CurrentCulture, "[{0},{1}]", a, b));
        }

#if SHO
        /// <summary>
        /// Produces a representation of the interval for the Python interactive console.
        /// </summary>
        /// <returns>A string representation of the inverval.</returns>
        [CLSCompliant(false)]
        public string __repr__ () {
            return(ToString());
        }
#endif

#if FUTURE
        // get the integers within a given range; this is used for fitting histograms to discrete distributions

        internal IEnumerable<int> GetContainedIntegers () {

            // lower limit is inclusive; for example: 3.0 -> 3, 3.1 -> 4
            int min = (int) Math.Ceiling(a);
            // upper limit is exclusive: for example: 3.9 -> 3, 4.0 -> 3
            int max = (int) Math.Floor(b);
            if (!(max < b)) max--;
            // iterate over integers in range
            for (int i = min; i <= max; i++) {
                yield return(i);
            }

        }
#endif

    }

}
