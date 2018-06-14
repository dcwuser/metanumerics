using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Text;

namespace Meta.Numerics.Functions {

    /// <summary>
    /// Represents a spinor.
    /// </summary>
    /// <remarks>
    /// <para>From a physicist's point of view, a spinor is an object with a particular quantum-mechanical spin. The quantum state of such
    /// an object is represented by a <see cref="SpinState"/> object.</para>
    /// <para>From a mathematician's point of view, a spinor labels an irreducible representation of the SO(3) or SU(2) Lie group.
    /// Individual vectors within each irreducible representation are represented by <see cref="SpinState"/> objects.</para>
    /// </remarks>
    public struct Spin : IEquatable<Spin> {

        // construction

        /// <summary>
        /// Instantiates a new spinor.
        /// </summary>
        /// <param name="j">The spin, which must be an integer or half-integer.</param>
        public Spin (double j) {

            // no negative (or too large) spins
            double tj = 2.0 * j;
            if ((tj < 0.0) || (tj > Int32.MaxValue)) throw new ArgumentOutOfRangeException(nameof(j));

            // spin must be integer or half-integer
            double tjt = Math.Floor(tj);
            if (tjt != tj) throw new ArgumentOutOfRangeException(nameof(j));

            // store 2*j
            twoJ = (int) tjt;
            
        }

        // stored data

        private int twoJ;

        // accessors

        /// <summary>
        /// Gets the spin of the spinor.
        /// </summary>
        public double J {
            get {
                return (twoJ / 2.0);
            }
        }

        internal int TwoJ {
            get {
                return (twoJ);
            }
        }

        // casting

        /*
        public static implicit operator Spin (double j) {
            return (new Spin(j));
        }

        public static implicit operator double (Spin s) {
            return (s.TwoJ / 2.0);
        }
        */

        // pre-defined

        /// <summary>
        /// Gets a spin-0 spinor.
        /// </summary>
        public static Spin SpinZero {
            get {
                return (new Spin(0.0));
            }
        }

        /// <summary>
        /// Gets a spin-1/2 spinor.
        /// </summary>
        public static Spin SpinOneHalf {
            get {
                return (new Spin(0.5));
            }
        }

        /// <summary>
        /// Gets a spin-1 spinor.
        /// </summary>
        public static Spin SpinOne {
            get {
                return (new Spin(1.0));
            }
        }

        // irrep data

        /// <summary>
        /// Gets the dimension of the spinor.
        /// </summary>
        public int Dimension {
            get {
                return (twoJ + 1);
            }
        }

        /// <summary>
        /// Returns the set of spinor states.
        /// </summary>
        /// <returns>An array of spin states that spans the spinor subspace.</returns>
        public SpinState[] States () {
            List<SpinState> states = new List<SpinState>();
            for (int i = -twoJ; i <= twoJ; i += 2) {
                states.Add(new SpinState(twoJ / 2.0, i / 2.0));
            }
            //for (int twoM = -TwoJ; twoM <= 0; twoM = twoM + 2) {
            //    states.Add(new SpinState(twoJ / 2.0, twoM / 2.0));
            //}
            return (states.ToArray());
        }

        // equality

        /// <summary>
        /// Determines whether two spinors are equal.
        /// </summary>
        /// <param name="a">The first spin.</param>
        /// <param name="b">The second spin.</param>
        /// <returns>True if <paramref name="a"/> and <paramref name="b"/> are equal, otherwise false.</returns>
        public static bool Equals (Spin a, Spin b) {
            return(a.TwoJ == b.TwoJ);
        }

        /// <summary>
        /// Determines whether two spinors are equal.
        /// </summary>
        /// <param name="a">The first spin.</param>
        /// <param name="b">The second spin.</param>
        /// <returns>True if <paramref name="a"/> and <paramref name="b"/> are equal, otherwise false.</returns>
        public static bool operator == (Spin a, Spin b) {
            return (Equals(a, b));
        }

        /// <summary>
        /// Determines whether two spinors are unequal.
        /// </summary>
        /// <param name="a">The first spin.</param>
        /// <param name="b">The second spin.</param>
        /// <returns>False if <paramref name="a"/> and <paramref name="b"/> are equal, otherwise true.</returns>
        public static bool operator != (Spin a, Spin b) {
            return (!Equals(a, b));
        }

        /// <summary>
        /// Determines whether the given spinor is equal to this one.
        /// </summary>
        /// <param name="s">The spinor to compare.</param>
        /// <returns>True if <paramref name="s"/> is equal to this one, otherwise false.</returns>
        public bool Equals (Spin s) {
            return (Equals(this, s));
        }

        /// <summary>
        /// Determines whether the given object represents the same spinor.
        /// </summary>
        /// <param name="obj">The object to compare.</param>
        /// <returns>True if <paramref name="obj"/> is an equal spin, otherwise false.</returns>
        public override bool Equals (object obj) {
            return ((obj is Spin) && Equals(this, (Spin)obj));
        }

        /// <summary>
        /// Computes a hash function for the spinor.
        /// </summary>
        /// <returns>An integer which is guaranteed equal for equal spinors, an unlikely to be equal for unequal spinors.</returns>
        public override int GetHashCode () {
            return (twoJ);
        }

        /// <summary>
        /// Produces a string representation of the spinor.
        /// </summary>
        /// <returns>A string representation of the spinor.</returns>
        public override string ToString () {
            if (twoJ % 2 == 0) {
                return ((twoJ / 2).ToString(CultureInfo.CurrentCulture));
            } else {
                return (twoJ.ToString(CultureInfo.CurrentCulture) + "/2");
            }
        }

    }

    

}
