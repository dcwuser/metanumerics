using System;

using Meta.Numerics;

namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {

        /// <summary>
        /// Computes the value of a spherical harmonic function.
        /// </summary>
        /// <param name="l">The order, which must be non-negative.</param>
        /// <param name="m">The sub-order, which must lie between -l and l inclusive.</param>
        /// <param name="theta">The azimuthal angle &#x3B8;. This angle is usually expressed as between -&#x3C0;/2 and +&#x3C0;/2, with positive values representing the upper hemisphere and negative values representing the lower hemisphere.</param>
        /// <param name="phi">The cylindrical angle &#x3C6;. This angle is usually expressed as between 0 and 2&#x3C0;, measured counter-clockwise (as seen from above) from the positive x-axis. It is also possible to use negative values to represent clockwise movement. </param>
        /// <returns>The value of Y<sub>l,m</sub>(&#x3B8;,&#x3C6;).</returns>
        public static Complex SphericalHarmonic (int l, int m, double theta, double phi) {
            if (l < 0) throw new ArgumentOutOfRangeException("l");
            if ((m > l) || (m < -l)) throw new ArgumentOutOfRangeException("m");
            if (m < 0) {
                Complex y = SphericalHarmonic(l, -m, theta, phi);
                if ((m % 2) != 0) y = -y;
                return (y.Conjugate);
            }

            double LP = OrthogonalPolynomials.LegendrePe(l, m, Math.Cos(theta));
            double mp = m * phi;
            return (new Complex(LP * AdvancedMath.Cos(mp,0.0), LP * AdvancedMath.Sin(mp,0.0)));

        }
        // Use: Spherically symmetric wave functions in QM; expansion in spherical coordinates
    }
}
