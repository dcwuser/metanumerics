using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Diagnostics;

using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents an experimental data point that is a function of an arbitrary variable.
    /// </summary>
    /// <typeparam name="T">The type of the ordinate (independent variable) characterizing the data point.</typeparam>
    public class UncertainMeasurement<T>  {

        private T x;
        private UncertainValue y;

        /// <summary>
        /// Initializes a new data point with the given values for the ordinate and uncertain abcissa.
        /// </summary>
        /// <param name="x">The ordinate.</param>
        /// <param name="y">The abcissa.</param>
        public UncertainMeasurement (T x, UncertainValue y) {
            this.x = x;
            this.y = y;
        }

        /// <summary>
        /// Initializes a new data point with the given values for the ordinate, abcissa, and uncertainty.
        /// </summary>
        /// <param name="x">The ordinate.</param>
        /// <param name="y">The best estimate of the abcissa.</param>
        /// <param name="dy">The uncertainty in the abcissa.</param>
        public UncertainMeasurement (T x, double y, double dy) {
            this.x = x;
            this.y = new UncertainValue(y, dy);
        }

        /// <summary>
        /// Gets or sets the value of the ordinate (independent variable).
        /// </summary>
        public T X {
            get {
                return (x);
            }
        }

        /// <summary>
        /// Gets or sets the uncertain value of the abcissa (the depdent variable).
        /// </summary>
        public UncertainValue Y {
            get {
                return (y);
            }
        }

        // equality

        private static bool Equals (UncertainMeasurement<T> d1, UncertainMeasurement<T> d2) {
            if (Object.ReferenceEquals(d1, null)) {
                if (Object.ReferenceEquals(d2, null)) {
                    return (true);
                } else {
                    return (false);
                }
            } else {
                if (Object.ReferenceEquals(d2, null)) {
                    return (false);
                } else {
                    return ((d1.X.Equals(d2.X)) && (d1.Y == d2.Y));
                }
            }
        }

        /// <summary>
        /// Determines whether two data points are equal.
        /// </summary>
        /// <param name="d1">The first data point.</param>
        /// <param name="d2">The second data point.</param>
        /// <returns>True if the data points are equal, otherwise false.</returns>
        public static bool operator == (UncertainMeasurement<T> d1, UncertainMeasurement<T> d2) {
            return (Equals(d1, d2));
        }

        /// <summary>
        /// Determines whether two data points are not equal.
        /// </summary>
        /// <param name="d1">The first data point.</param>
        /// <param name="d2">The second data point.</param>
        /// <returns>True if the data points are not equal, otherwise false.</returns>
        public static bool operator != (UncertainMeasurement<T> d1, UncertainMeasurement<T> d2) {
            return (!Equals(d1, d2));
        }

        /// <summary>
        /// Determines whether the object represents the same data point.
        /// </summary>
        /// <param name="obj">The object.</param>
        /// <returns>True if the object represents the same data point, otherwise false.</returns>
        public override bool Equals (object obj) {
            UncertainMeasurement<T> d = obj as UncertainMeasurement<T>;
            if (Object.ReferenceEquals(d, null)) {
                return (false);
            } else {
                return (Equals(this, d));
            }
        }

        /// <summary>
        /// Gets a hash code for the data point.
        /// </summary>
        /// <returns>A hash code for the data point.</returns>
        public override int GetHashCode () {
            return (x.GetHashCode() ^ y.GetHashCode());
        }

    }

}
