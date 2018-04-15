using System;
using System.Diagnostics;
using System.Collections.Generic;

namespace Meta.Numerics {
    
    internal static class Global {

        // maximum number of iterations of a series
        public const int SeriesMax = 250;

        // double dedicates 10 bits to the magnitude of the the exponent so 2^10 = 1024 is the largest exponent,
        // and 2^1024 is the largest representable double. 2^512 is its square root, which we compare to to
        // decide if we are in danger of overflowing
        public static readonly double SqrtMax = Math.Pow(2.0, 512);

        // double dedicates 52 bits to the magnitude of the mantissa, so 2^-52 is the smallest fraction difference
        // it can detect; in order to avoid any funny effects at the margin, we try for one byte less, 2^-49
        public static readonly double Accuracy = Math.Pow(2.0, -49);

        public static void Swap<T> (ref T a, ref T b) {
            T t = a;
            a = b;
            b = t;
        }

        // pre-calculate some often-used factors
        // these are actually useless; the compiler pre-calculates them for us

        public const double TwoPI = 2.0 * Math.PI;

        public const double HalfPI = Math.PI / 2.0;

        // pre-calculate some often-used square roots

        public static readonly double SqrtTwo = Math.Sqrt(2.0);

        public static readonly double SqrtThree = Math.Sqrt(3.0);

        // Math.Sqrt(Math.PI) is off by 2 x 10^{-16} from nearest representable value
        public const double SqrtPI = 1.7724538509055160273;

        // Math.Sqrt(2.0 * Math.PI) is off by 2 X 10^{-16} from nearest representable value
        public const double SqrtTwoPI = 2.5066282746310005024;

        // Math.Sqrt(Math.PI / 2.0) is off by 2 X 10^{-16} from nearest representable value
        public const double SqrtHalfPI = 1.2533141373155002512;

        // pre-calculate some often-used logs

        public static readonly double LogTwo = Math.Log(2.0);

        // IList defines CopyTo. IReadOnlyList doesn't. So in order to switch from IList to IReadOnlyList,
        // we define an extension method to implement it.
        public static void CopyTo<T> (this IReadOnlyList<T> source, T[] target, int startIndex) {
            Debug.Assert(source != null);
            Debug.Assert(target != null);
            Debug.Assert(startIndex >= 0);
            for (int i = 0; i < source.Count; i++) {
                target[startIndex + i] = source[i];
            }
        }

    }

}
