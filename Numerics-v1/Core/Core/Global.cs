using System;

namespace Meta.Numerics {
    
    internal static class Global {

        // maximum number of iterations of a series
        public const int SeriesMax = 250;

        // double dedicates 10 bits to the magnitude of the the exponent so 2^10 = 1024 is the largest exponent,
        // and 2^1024 is the largest representable double. 2^512 is its square root, which we compare to to
        // decide if we are in danger of overflowing
        public static readonly double SqrtMax = Math.Pow(2.0, 512);

        // double dedicates 52 bits to the magnitude of the mantissa, so 2^-52 is the smallest fraction difference
        // it can detect; in order to avoid any funny effects at the margin, we only try for 2^-50
        public static readonly double Accuracy = Math.Pow(2.0, -50);

        public static void Swap<T> (ref T a, ref T b) {
            T t = a;
            a = b;
            b = t;
        }

        // pre-calculate some often-used factors

        public const double TwoPI = 2.0 * Math.PI;

        public const double HalfPI = Math.PI / 2.0;

        // pre-calculate some often-used square roots

        public static readonly double SqrtTwo = Math.Sqrt(2.0);

        public static readonly double SqrtThree = Math.Sqrt(3.0);

        public static readonly double SqrtPI = Math.Sqrt(Math.PI);

        public static readonly double SqrtTwoPI = Math.Sqrt(2.0 * Math.PI);

        // pre-calculate some often-used logs

        public static readonly double LogTwo = Math.Log(2.0);

    }

}
