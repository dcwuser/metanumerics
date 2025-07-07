using System;
using System.Diagnostics;
using System.Globalization;

namespace Meta.Numerics {

    /// <summary>
    /// Represents an interval on the integers.
    /// </summary>
    public readonly struct DiscreteInterval : IEquatable<DiscreteInterval> {

        internal DiscreteInterval (int min, int max) {
            Debug.Assert(max >= min);
            if (min == Int32.MinValue) throw new ArgumentOutOfRangeException(nameof(min));
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
        /// <remarks><para>In order that width calculations always fit within a <see cref="System.UInt32"/>, we limit discrete intervals
        /// to the symmetric range of signed integers -<see cref="Int32.MaxValue"/> to +<see cref="Int32.MaxValue"/>. That is,
        /// the value <see cref="Int32.MinValue"/> is not allowed. This should almost never have a practical effect, unless you
        /// attempt to instantiate the full interval. For such a purpose, use <see cref="DiscreteInterval.Infinite"/> instead.</para></remarks>
        public static DiscreteInterval FromEndpoints (int a, int b) {
            if (a > b) {
                return new DiscreteInterval(b, a);
            } else {
                return new DiscreteInterval(a, b);
            }
        }

        /// <summary>
        /// Gets the minimum value on the interval.
        /// </summary>
        public int LeftEndpoint {
            get {
                return min;
            }
        }

        /// <summary>
        /// Gets the maximum value on the interval.
        /// </summary>
        public int RightEndpoint {
            get {
                return max;
            }
        }

        /// <summary>
        /// Gets the width of the interval.
        /// </summary>
        /// <remarks>
        /// <para>The width of the interval is the difference between <see cref="LeftEndpoint"/> and <see cref="RightEndpoint"/>.
        /// This is one less than the number of integers in the interval. Thus an interval with equal left and right endpoints
        /// has width 0, not width 1. And and the interval {0 .. MaxValue} has width MaxValue, not MaxValue + 1.</para>
        /// </remarks>
        internal uint Width {
            get {
                return (uint) (max - min);
            }
        }

        /// <summary>
        /// The semi-infinite discrete interval representing all non-negative integers.
        /// </summary>
        public static readonly DiscreteInterval Semiinfinite = new DiscreteInterval(0, Int32.MaxValue);

        /// <summary>
        /// The infinite discrete interval representing all integers.
        /// </summary>
        public static readonly DiscreteInterval Infinite = new DiscreteInterval(-Int32.MaxValue, Int32.MaxValue);

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
            return unchecked(min.GetHashCode() + 31 * max.GetHashCode());
        }

        /// <inheritdoc />
        public override string ToString () {
            return ToString(CultureInfo.CurrentCulture);
        }

        private string ToString (IFormatProvider provider) {
            return String.Format(provider, "{{{0}..{1}}}", min, max);
        }

    }
}
