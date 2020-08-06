using System;
using System.Diagnostics;
using System.Globalization;

namespace Meta.Numerics {

    /// <summary>
    /// Represents an interval on the integers.
    /// </summary>
    public struct DiscreteInterval : IEquatable<DiscreteInterval> {

        internal DiscreteInterval (int min, int max) {
            Debug.Assert(max >= min);
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

        // Equality

        private static bool Equals (DiscreteInterval u, DiscreteInterval v) {
            return (u.min == v.min) && (u.max == v.max);
        }


        /// <summary>
        /// Tests whether two discrete intervals are equal.
        /// </summary>
        /// <param name="a">The first discrete interval.</param>
        /// <param name="b">The second discrete interval.</param>
        /// <returns><see langword="true"/> if <paramref name="a"/> and <paramref name="b"/> are equal, otherwise <see langword="false"/>.</returns>
        public static bool operator == (DiscreteInterval a, DiscreteInterval b) {
            return Equals(a, b);
        }

        /// <summary>
        /// Tests whether two intervals are not equal.
        /// </summary>
        /// <param name="a">The first interval.</param>
        /// <param name="b">The second interval.</param>
        /// <returns><see langword="false"/> if <paramref name="a"/> and <paramref name="b"/> are equal, otherwise <see langword="true"/>.</returns>
        public static bool operator != (DiscreteInterval a, DiscreteInterval b) {
            return !Equals(a, b);
        }

        /// <summary>
        /// Tests wither the current instance is equal to another discrete interval.
        /// </summary>
        /// <param name="other">Another discrete interval.</param>
        /// <returns><see langword="true"/> if the current instance equals <paramref name="other"/>, otherwise <see langword="false"/>.</returns>
        public bool Equals (DiscreteInterval other) {
            return Equals(this, other);
        }

        /// <summary>
        /// Tests whether a given object is equal to the current interval.
        /// </summary>
        /// <param name="obj">An object.</param>
        /// <returns><see langword="true"/> if <paramref name="obj"/> is an equal <see cref="DiscreteInterval"/>, otherwise <see langword="false"/>.</returns>
        public override bool Equals (object obj) {
            if (obj is DiscreteInterval other) {
                return Equals(this, other);
            } else {
                return false;
            }
        }

        /// <inheritdoc />
        public override int GetHashCode () {
            unchecked {
                return (min.GetHashCode() + 31 * max.GetHashCode());
            }
        }

        /// <inheritdoc />
        public override string ToString () {
            return ToString(CultureInfo.CurrentCulture);
        }

        private string ToString (IFormatProvider format) {
            return String.Format(format, "{{{0}..{1}}}", min, max);
        }

    }
}
