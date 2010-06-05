using System;

namespace Meta.Numerics {

    public static class MoreMath {

        /// <summary>
        /// Rasises an argument to an integer power.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <param name="n">The power.</param>
        /// <returns>The value of x<sup>n</sup>.</returns>
        /// <remarks>
        /// <para>Low integer powers can be computed by optimized algorithms much faster than the general
        /// alrogithm for an arbitrary real power employed by <see cref="System.Math.Pow"/>.</para>
        /// </remarks>
        public static double Pow (double x, int n) {

            if (n < 0) return (1.0 / Pow(x, -n));

            switch(n) {
                case 0:
                    // we follow convention that 0^0 = 1
                    return (1.0);
                case 1:
                    return (x);
                case 2:
                    // 1 multiply
                    return (x * x);
                case 3:
                    // 2 multiplies
                    return (x * x * x);
                case 4: {
                        // 2 multiplies
                        double x2 = x * x;
                        return (x2 * x2);
                    }
                case 5: {
                        // 3 multiplies
                        double x2 = x * x;
                        return (x2 * x2 * x);
                    }
                case 6: {
                        // 3 multiplies
                        double x2 = x * x;
                        return (x2 * x2 * x2);
                    }
                case 7: {
                        // 4 multiplies
                        double x3 = x * x * x;
                        return (x3 * x3 * x);
                    }
                case 8: {
                        // 3 multiplies
                        double x2 = x * x;
                        double x4 = x2 * x2;
                        return (x4 * x4);
                    }
                case 9: {
                        // 4 multiplies
                        double x3 = x * x * x;
                        return (x3 * x3 * x3);
                    }
                case 10: {
                        // 4 multiplies
                        double x2 = x * x;
                        double x4 = x2 * x2;
                        return (x4 * x4 * x2);
                    }
                default:
                    return (Math.Pow(x, n));
            }


        }

        // an internal method for squaring

        internal static double Pow2 (double x) {
            return (x * x);
        }

        /// <summary>
        /// Computes the length of a right triangle's hypotenuse.
        /// </summary>
        /// <param name="x">The length of one side.</param>
        /// <param name="y">The length of another side.</param>
        /// <returns>The length of the hypotenuse, sqrt(x<sup>2</sup> + y<sup>2</sup>).</returns>
        /// <remarks>
        /// <para>The length is computed accurately, even in cases where
        /// x<sup>2</sup> or y<sup>2</sup> would overflow.</para>
        /// </remarks>
        public static double Hypot (double x, double y) {
            if ((x == 0.0) && (y == 0.0)) {
                return (0.0);
            } else {
                double ax = Math.Abs(x);
                double ay = Math.Abs(y);
                if (ax > ay) {
                    double r = y / x;
                    return (ax * Math.Sqrt(1.0 + r * r));
                } else {
                    double r = x / y;
                    return (ay * Math.Sqrt(1.0 + r * r));
                }
            }
        }

    }


}