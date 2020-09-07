using System;
using System.Diagnostics;
using System.Globalization;

namespace Meta.Numerics {

    /// <summary>
    /// Represents an interval on the real number line.
    /// </summary>
    /// <remarks>
    /// <para>Use the static methods <see cref="FromEndpoints"/>, <see cref="FromMidpointAndWidth"/>,
    /// and <see cref="FromEndpointAndWidth"/> to instantiate intervals.</para>
    /// </remarks>
    public readonly struct Interval : IEquatable<Interval> {

        private readonly double a, b;

        private Interval (double a, double b) {
            Debug.Assert(a <= b);
            this.a = a;
            this.b = b;
        }

        /// <summary>
        /// Gets the left (lower) endpoint of the interval.
        /// </summary>
        public double LeftEndpoint {
            get {
                return a;
            }
        }

        /// <summary>
        /// Gets the right (upper) endpoint of the interval.
        /// </summary>
        public double RightEndpoint {
            get {
                return b;
            }
        }

        /// <summary>
        /// Determines whether a given point is contained in the closed interval.
        /// </summary>
        /// <param name="x">The point.</param>
        /// <returns><see langword="true"/> if <paramref name="x"/> is contained in the closed interval [a,b],
        /// otherwise <see langword="false"/>.</returns>
        public bool Contains (double x) {
            return (x >= a) && (x <= b);
        }

        /// <summary>
        /// Determines whether a given point is contained in the closed or open interval.
        /// </summary>
        /// <param name="x">The point.</param>
        /// <param name="type">Specifier of whether the interval should be considered closed or open.</param>
        /// <returns><see langword="true"/> if <paramref name="x"/> is contained in the interval as specified,
        /// otherwise <see langword="false"/></returns>
        public bool Contains (double x, IntervalType type) {
            if (type == IntervalType.Closed) {
                return (x >= a) && (x <= b);
            } else {
                return (x > a) && (x < b);
            }
        }

        /// <summary>
        /// Determines whether a given point is contained in the interval, with the left and right endpoints
        /// separately specififed as closed or open.
        /// </summary>
        /// <param name="x">The point.</param>
        /// <param name="leftEndpointType">Specifier of whether the interval should be considered as including its left endpoint.</param>
        /// <param name="rightEndpointType">Specifier of whether the interval should be considered as including its right endpoint.</param>
        /// <returns><see langword="true"/> if <paramref name="x"/> is contained in the interval as specified,
        /// otherwise <see langword="false"/></returns>
        public bool Contains (double x, IntervalType leftEndpointType, IntervalType rightEndpointType) {
            bool leftSatisfied = leftEndpointType == IntervalType.Closed ? x >= a : x > a;
            bool rightSatisfied = rightEndpointType == IntervalType.Closed ? x <= b : x < b;
            return leftSatisfied && rightSatisfied;
        }

        /// <summary>
        /// Determines whether the argument lies in the open interval.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>True if <paramref name="x"/> lies in (a,b), otherwise False.</returns>
        public bool OpenContains (double x) {
            return (x > a) && (x < b);
        }
        /// <summary>
        /// Determines whether the argument lies in the closed interval.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>True if <paramref name="x"/> lies in [a,b], otherwise False.</returns>
        public bool ClosedContains (double x) {
            return (x >= a) && (x <= b);
        }

        /// <summary>
        /// Gets the width of the interval.
        /// </summary>
        public double Width {
            get {
                return b - a;
            }
        }

        /// <summary>
        /// Gets the mid-point of the interval.
        /// </summary>
        public double Midpoint {
            get {
                return 0.5 * (a + b);
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
                return new Interval(a, b);
            } else {
                return new Interval(b, a);
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
                return new Interval(endpoint + width, endpoint);
            } else {
                return new Interval(endpoint, endpoint+width);
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
                return FromMidpointAndWidth(midpoint, -width);
            } else {
                return new Interval(midpoint - 0.5 * width, midpoint + 0.5 * width);
            }
        }

        // equality

        private static bool Equals (Interval u, Interval v) {
            return (u.a == v.a) && (u.b == v.b);
        }

        /// <summary>
        /// Tests wither the current instance is equal to another interval.
        /// </summary>
        /// <param name="other">Another iteterval.</param>
        /// <returns><see langword="true"/> if the current instance equals <paramref name="other"/>, otherwise <see langword="false"/>.</returns>
        public bool Equals (Interval other) {
            return Equals(this, other);
        }

        /// <summary>
        /// Tests whether two intervals are equal.
        /// </summary>
        /// <param name="u">The first interval.</param>
        /// <param name="v">The second interval.</param>
        /// <returns><see langword="true"/> if <paramref name="u"/> and <paramref name="v"/> are equal, otherwise <see langword="false"/>.</returns>
        public static bool operator == (Interval u, Interval v) {
            return Equals(u, v);
        }

        /// <summary>
        /// Tests whether two intervals are not equal.
        /// </summary>
        /// <param name="u">The first interval.</param>
        /// <param name="v">The second interval.</param>
        /// <returns><see langword="false"/> if <paramref name="u"/> and <paramref name="v"/> are equal, otherwise <see langword="true"/>.</returns>
        public static bool operator != (Interval u, Interval v) {
            return !Equals(u, v);
        }

        /// <summary>
        /// Tests whether a given object is equal to the current interval.
        /// </summary>
        /// <param name="obj">An object.</param>
        /// <returns><see langword="true"/> if <paramref name="obj"/> is an equal <see cref="Interval"/>, otherwise <see langword="false"/>.</returns>
        public override bool Equals (object obj) {
            if (obj is Interval other) {
                return Equals(this, other);
            } else {
                return false;
            }
        }

        /// <inheritdoc />
        public override int GetHashCode () {
            return unchecked(a.GetHashCode() + 31 * b.GetHashCode());
        }

        /// <summary>
        /// Produces a string representation of the interval.
        /// </summary>
        /// <returns>A string representation of the interval.</returns>
        public override string ToString () {
            return ToString(CultureInfo.CurrentCulture);
        }

        private string ToString (IFormatProvider provider) {
            return String.Format(provider, "[{0},{1}]", a, b);
        }

        /// <summary>
        /// The entire number line.
        /// </summary>
        public static readonly Interval Infinite = new Interval(Double.NegativeInfinity, Double.PositiveInfinity);

        /// <summary>
        /// The number line from zero to positive infinity.
        /// </summary>
        public static readonly Interval Semiinfinite = new Interval(0.0, Double.PositiveInfinity);

        /// <summary>
        /// The unit interval.
        /// </summary>
        public static readonly Interval Unit = new Interval(0.0, 1.0);

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

    /// <summary>
    /// Indicates whether an interval should be considered closed or open.
    /// </summary>
    public enum IntervalType {

        /// <summary>
        /// Endpoints should be considered within the interval.
        /// </summary>
        Closed,

        /// <summary>
        /// Endpoints should be considered outside the interval.
        /// </summary>
        Open
    }

}
