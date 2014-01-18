using System;


namespace Meta.Numerics.Functions {

    // code to invert beta, gamma, and error functions
    // this is unfortunately a big pain
#if PAST
    public static partial class AdvancedMath {

        // the probit function is the inverse CDF of the normal distribution

        internal static double ApproximateProbit (double P) {

            if (P < 0.1) {
                return (-Global.SqrtTwo * ApproximateInverseErfc(2.0 * P));
            } else if (P < 0.9) {
                return (Global.SqrtTwo * ApproximateInverseErf(2.0 * P - 1.0));
            } else {
                return (Global.SqrtTwo * ApproximateInverseErfc(2.0 * (1.0 - P)));
            }

        }
        


    }
#endif

}
