using System;

namespace Meta.Numerics.Functions {

    /// <summary>
    /// Contains methods that compute advanced functions with real arguments.
    /// </summary>
    public static partial class AdvancedMath {

        // members are defined in other files

        // internal utility trig functions that are accurate for large arguments

        internal static double Sin (double x, double y) {
            return (Math.Sin(Reduce(x,y)));
        }

        internal static double Cos (double x, double y) {
            return (Math.Cos(Reduce(x,y)));
        }

        // reduces an argument to its corresponding argument between -2 Pi < x < 2 Pi

        internal static double Reduce (double x, double y) {

            double t = x + 2.0 * Math.PI * y;
            if ((Math.Abs(t) < 64.0) || (Math.Abs(t) > dmax)) {
                // if the argument is small we don't need the high accurary reduction
                // if the argument is too big, we can't do the high accuracy reduction because it would overflow a decimal vairable
                return (t);
            } else {
                // otherwise, convert to decimal, subtract a multiple of 2 Pi, and return

                // reduce x by factors of 2 Pi
                decimal dx = Convert.ToDecimal(x);
                decimal dn = Decimal.Truncate(dx / dPI2);
                dx = dx - dn * dPI2;

                // reduce y by factors of 1
                decimal dy = Convert.ToDecimal(y);
                decimal dm = Decimal.Truncate(dy / 1.0m);
                dy = dy - dm * 1.0m;

                // form the argument
                decimal dt = dx + dy * dPI2;
                return (Convert.ToDouble(dt));

            }
        }

        private static readonly decimal dPI2 = 2.0m * 3.1415926535897932384626433832795m;

        private static readonly double dmax = Convert.ToDouble(Decimal.MaxValue);

    }

}
