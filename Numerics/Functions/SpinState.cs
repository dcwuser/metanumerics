using System;
using System.Diagnostics;
using System.Globalization;
using System.Text;

namespace Meta.Numerics.Functions {

    /// <summary>
    /// Represents the state of a spinor.
    /// </summary>
    public struct SpinState : IEquatable<SpinState> {

        /// <summary>
        /// Instantiates a new SpinState with the given spin and magnetic quantum numbers.
        /// </summary>
        /// <param name="j">The spin number.</param>
        /// <param name="m">The magnetic number.</param>
        public SpinState (double j, double m) : this(new Spin(j), m) { }

        /// <summary>
        /// Instantiates a new SpinState with the given spin and magnetic quantum number.
        /// </summary>
        /// <param name="s">The spin.</param>
        /// <param name="m">The magnetic quantum number.</param>
        public SpinState (Spin s, double m) {

            spin = s;

            // 2M must be an integer
            double tm = 2.0 * m;
            double tmt = Math.Floor(tm);
            if (tmt != tm) throw new ArgumentOutOfRangeException(nameof(m));

            int tmti = (int)tmt;

            twoM = (int)tmt;
            // -J <= M <= J
            if (Math.Abs(tmti) > s.TwoJ) throw new ArgumentOutOfRangeException(nameof(m));

            // half-integer J requires half-integer M; integer J requires integer M
            if ((s.TwoJ % 2) != Math.Abs(twoM % 2)) throw new ArgumentOutOfRangeException(nameof(m));

        }

        private Spin spin;

        private int twoM;

        /// <summary>
        /// Gets the spin value of the spin state.
        /// </summary>
        public double J {
            get {
                return (spin.J);
            }
        }

        /// <summary>
        /// Gets the magnetic substate value of the spin state.
        /// </summary>
        public double M {
            get {
                return (twoM / 2.0);
            }
        }

        internal int TwoJ {
            get {
                return (spin.TwoJ);
            }
        }

        internal int TwoM {
            get {
                return (twoM);
            }
        }

        internal int JPlusM {
            get {
                return ((spin.TwoJ + twoM) / 2);
            }
        }

        internal int JMinusM {
            get {
                return ((spin.TwoJ - twoM) / 2);
            }
        }

        internal SpinState Invert () {
            twoM = -twoM;
            return (this);
        }

        /// <summary>
        /// Gets the spinor representation to which the spin state belongs.
        /// </summary>
        public Spin Representation {
            get {
                return (spin);
            }
        }

        // equality

        private static bool Equals (SpinState s1, SpinState s2) {
            return (s1.spin == s2.spin) && (s1.twoM == s2.twoM);
        }

        /// <summary>
        /// Tests whether another instance equals the current instance.
        /// </summary>
        /// <param name="other">The other <see cref="SpinState"/> instance.</param>
        /// <returns><see langword="true"/> if the current instance is equal to  <paramref name="other"/>, otherwise <see langword="false"/>.</returns>
        public bool Equals (SpinState other) {
            return Equals(this, other);
        }

        /// <summary>
        /// Determines whether two spin states are equal.
        /// </summary>
        /// <param name="s1">The first spin state.</param>
        /// <param name="s2">The second spin state.</param>
        /// <returns><see langword="true"/> if <paramref name="s1"/> and <paramref name="s2"/> are equal, otherwise <see langword="true"/>.</returns>
        public static bool operator == (SpinState s1, SpinState s2) {
            return Equals(s1, s2);
        }

        /// <summary>
        /// Determines whether two spin states are unequal.
        /// </summary>
        /// <param name="s1">The first spin state.</param>
        /// <param name="s2">The second spin state.</param>
        /// <returns><see langword="false"/> if <paramref name="s1"/> and <paramref name="s2"/> are equal, otherwise <see langword="true"/>.</returns>
        public static bool operator != (SpinState s1, SpinState s2) {
            return !Equals(s1, s2);
        }

        /// <summary>
        /// Determines whether the given object represents the same spin state.
        /// </summary>
        /// <param name="obj">The object to compare.</param>
        /// <returns><see langword="true"/> if <paramref name="obj"/> is an equal spin state, otherwise <see langword="false"/>.</returns>
        public override bool Equals (object obj) {
            return (obj is SpinState) && Equals(this, (SpinState) obj);
        }

        /// <summary>
        /// Computes a hash function for the spin state.
        /// </summary>
        /// <returns>An integer which is guaranteed equal for equal spin states, an unlikely to be equal for unequal spin states.</returns>
        public override int GetHashCode () {
            // left-shift 2j by about half the width of an int (which is 32 bits) so 2j and 2m values are unlikely to interfere
            return (TwoJ << 14) ^ (TwoM);
        }

        /// <summary>
        /// Produces a string representation of the spin state.
        /// </summary>
        /// <returns>A string representation of the spin state.</returns>
        public override string ToString () {
            return ToString(CultureInfo.CurrentCulture);
        }

        private string ToString (IFormatProvider provider) {
            NumberFormatInfo numberFormat = (NumberFormatInfo) provider.GetFormat(typeof(NumberFormatInfo));
            StringBuilder text = new StringBuilder(spin.ToString());
            text.Append(",");
            if (twoM > 0) {
                text.Append(numberFormat.PositiveSign);
            } else if (twoM < 0) {
                text.Append(numberFormat.NegativeSign);
            }
            text.Append(SpinString(Math.Abs(twoM)));
            return text.ToString();
        }

        private static string SpinString (int twoJ) {
            if (twoJ % 2 == 0) {
                return ((twoJ / 2).ToString(CultureInfo.CurrentCulture));
            } else {
                return (twoJ.ToString(CultureInfo.CurrentCulture) + "/2");
            }
        }

    }

}
