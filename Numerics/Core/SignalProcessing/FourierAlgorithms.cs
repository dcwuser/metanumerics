using System;
using System.Collections.Generic;

using Meta.Numerics;

namespace Meta.Numerics.SignalProcessing {

    /// <summary>
    /// Specifies the normalization convention to be used in a Fourier transform.
    /// </summary>
    /// <remarks>
    /// <para>The most common convention in signal processing applications is <see cref="FourierNormalization.None"/>.</para>
    /// </remarks>
    public enum FourierNormalization {
        None, Unitary, Inverse
    }

    /// <summary>
    /// Specifies the sign convention to be used in the exponent of a Fourier transform.
    /// </summary>
    /// <remarks>
    /// <para>The most common convention in signal processing applications is <see cref="FourierSign.Negative"/>.</para>
    /// </remarks>
    public enum FourierSign {

        /// <summary>
        /// The exponent has positive imaginary values.
        /// </summary>
        Positive,

        /// <summary>
        /// The exponent has negative imaginary values.
        /// </summary>
        Negative
    }

    /// <summary>
    /// An engine for performing Fourier transforms on complex series.
    /// </summary>
    /// <seealso href="http://en.wikipedia.org/wiki/Discrete-time_Fourier_transform"/>
    public class FourierTransformer {

        /// <summary>
        /// Initializes a new instance of the Fourier transformer.
        /// </summary>
        /// <param name="size">The series length of the transformer, which must be positive.</param>
        public FourierTransformer (int size) : this(size, FourierSign.Negative, FourierNormalization.None) {
        }

        /// <summary>
        /// Initializes a new instance of the Fourier transformer with the given sign and normalization conventions.
        /// </summary>
        /// <param name="size">The series length of the transformer, which must be positive.</param>
        /// <param name="signConvention">The sign convention of the transformer.</param>
        /// <param name="normalizationConvention">The normalization convention of the transformer.</param>
        public FourierTransformer (int size, FourierSign signConvention, FourierNormalization normalizationConvention) {
            if (size < 1) throw new ArgumentOutOfRangeException("size");
            this.size = size;
            this.factors = Factor.Factorize(size);
            this.signConvention = signConvention;
            this.normalizationConvention = normalizationConvention;
            this.roots = FourierAlgorithms.ComputeRoots(size, +1);
        }

        private int size;
        private List<Factor> factors;
        private FourierNormalization normalizationConvention;
        private FourierSign signConvention;
        private Complex[] roots;

        /// <summary>
        /// The series length for which the transformer is specialized.
        /// </summary>
        public int Length {
            get {
                return (size);
            }
        }

        /// <summary>
        /// Gets the normalization convention used by the transformer.
        /// </summary>
        public FourierNormalization NormalizationConvention {
            get {
                return (normalizationConvention);
            }
        }

        /// <summary>
        /// Gets the normalization convention used by the transformer.
        /// </summary>
        public FourierSign SignConvention {
            get {
                return (signConvention);
            }
        }

        private int GetSign () {
            if (signConvention == FourierSign.Positive) {
                return (+1);
            } else {
                return (-1);
            }
        }

        private static void Normalize (Complex[] x, double f) {
            for (int i = 0; i < x.Length; i++) {
                x[i] = new Complex(f * x[i].Re, f * x[i].Im);
            }
        }


        /// <summary>
        /// Computes the Fourier transform of the given series.
        /// </summary>
        /// <param name="values">The series to transform.</param>
        /// <returns>The discrete Fourier transform of the series.</returns>
        public Complex[] Transform (IList<Complex> values) {
            if (values == null) throw new ArgumentNullException("values");
            if (values.Count != size) throw new InvalidOperationException();

            // copy the original values into a new array
            Complex[] x = new Complex[size];
            values.CopyTo(x, 0);

            // normalize the copy appropriately
            if (normalizationConvention == FourierNormalization.Unitary) {
                Normalize(x, 1.0 / Math.Sqrt(size));
            } else if (normalizationConvention == FourierNormalization.Inverse) {
                Normalize(x, 1.0 / size);
            }

            // create a scratch array
            Complex[] y = new Complex[size];

            // do the FFT
            FourierAlgorithms.Fft(values.Count, factors, ref x, ref y, roots, GetSign());

            return (x);

        }

        /// <summary>
        /// Computes the inverse Fourier transform of the given series.
        /// </summary>
        /// <param name="values">The series to invert.</param>
        /// <returns>The inverse discrete Fourier transform of the series.</returns>
        public Complex[] InverseTransform (IList<Complex> values) {
            if (values == null) throw new ArgumentNullException("values");
            if (values.Count != size) throw new InvalidOperationException();

            // copy the original values into a new array
            Complex[] x = new Complex[size];
            values.CopyTo(x, 0);

            // normalize the copy appropriately
            if (normalizationConvention == FourierNormalization.None) {
                Normalize(x, 1.0 / size);
            } else if (normalizationConvention == FourierNormalization.Unitary) {
                Normalize(x, 1.0 / Math.Sqrt(size));
            }

            // create a scratch array
            Complex[] y = new Complex[size];

            // do the FFT
            FourierAlgorithms.Fft(values.Count, factors, ref x, ref y, roots, -GetSign());

            return (x);

        }

    }

    internal static class FourierAlgorithms {

        // computes the Nth roots of unity, which are the factors in a length-N Fourier transform

        public static Complex[] ComputeRoots (int N, int sign) {
            Complex[] u = new Complex[N + 1];
            double t = sign * Global.TwoPI / N;
            u[0] = 1.0;
            for (int r = 1; r < N; r++) {
                double rt = r * t;
                u[r] = new Complex(Math.Cos(rt), Math.Sin(rt));
            }
            u[N] = 1.0;
            return (u);
        }

        // does a FFT on x
        // the length N is factored as factors
        // y is scratch space, u are the Nth roots of unity, sign is the sign convention for the exponent

        public static void Fft (int N, List<Factor> factors, ref Complex[] x, ref Complex[] y, Complex[] u, int sign) {
            // keep track of 
            int Ns = 1;
            // loop over factors
            foreach (Factor factor in factors) {
                for (int m = 0; m < factor.Multiplicity; m++) {
                    FftPass(N, factor.Value, Ns, x, y, u, sign);
                    Ns *= factor.Value;
                    // exchange x and y, so result is now in x and y is ready to be next target
                    // this is the reason that x and y have the ref modifier; if they did not the switch would just be internal
                    // and the caller would have to figure out which array was the final target to get the result
                    Complex[] t = x; x = y; y = t;
                }
            }
        }

        private static void FftPass (int N, int R, int Ns, Complex[] x, Complex[] y, Complex[] u, int sign) {
            Complex[] v = new Complex[R];
            int dx = (N / R);
            for (int j = 0; j < N / R; j++) {
                // the j-loop can be parallelized (if the v buffer is not shared)
                int xi = j;
                int ui = 0; if (sign < 0) ui = N;
                int du = (dx / Ns) * (j % Ns); if (sign < 0) du = -du;
                    v[0] = x[xi]; // this could be in the r-loop, but by pulling it out we avoid a 6-op complex multiply by 1
                    for (int r = 1; r < R; r++) {
                        xi += dx;
                        ui += du;
                        v[r] = x[xi] * u[ui];
                    }
                    int y0 = Expand(j, Ns, R);
                    Fft(R, v, y, y0, Ns, sign);
            }
        }

        private static int Expand (int idxL, int N1, int N2) {
            return ((idxL / N1) * N1 * N2 + (idxL % N1));
        }

        private static void Fft (int R, Complex[] v, Complex[] y, int y0, int dy, int sign) {
            switch (R) {
                case 2:
                    FftKernel2(v, y, y0, dy);
                    break;
                case 3:
                    FftKernel3(v, y, y0, dy, sign);
                    break;
                case 5:
                    FftKernel5(v, y, y0, dy, sign);
                    break;
                default:
                    FftKernel(R, v, y, y0, dy, sign);
                    break;
            }
            /*
            if (R == 2) {
                FftKernel2(v, y, y0, dy);
            } else if (R == 3) {
                FftKernel3(v, y, y0, dy, sign);
            } else if (R == 5) {
                //FftKernel(R, v, y, y0, dy, sign);
                FftKernel5(v, y, y0, dy, sign);
            } else {
                FftKernel(R, v, y, y0, dy, sign);
            }
            */
        }

        private static void FftKernel2 (Complex[] v, Complex[] y, int y0, int dy) {
            double a0 = v[0].Re; double b0 = v[0].Im;
            double a1 = v[1].Re; double b1 = v[1].Im;
            y[y0] = new Complex(a0 + a1, b0 + b1);
            y[y0 + dy] = new Complex(a0 - a1, b0 - b1);
            // for some reason, this looks to be faster than using the complex add and subtract; i don't see why
            //u[u0] = v[0] + v[1];
            //u[u0 + du] = v[0] - v[1];
        }

        // This length-3 FFT kernel requires just 12 additions and 4 multiplications, or 16 operations total.
        // This is half as many operations as is required by the naive implementation. 

        private static void FftKernel3 (Complex[] v, Complex[] y, int y0, int dy, int sign) {
            double a12p = v[1].Re + v[2].Re;
            double b12p = v[1].Im + v[2].Im;
            double sa = v[0].Re + r31.Re * a12p;
            double sb = v[0].Im + r31.Re * b12p;
            double ta = r31.Im * (v[1].Re - v[2].Re);
            double tb = r31.Im * (v[1].Im - v[2].Im);
            if (sign < 0) { ta = -ta; tb = -tb; }
            y[y0] = new Complex(v[0].Re + a12p, v[0].Im + b12p);
            int y1 = y0 + dy;
            y[y1] = new Complex(sa - tb, sb + ta);
            int y2 = y1 + dy;
            y[y2] = new Complex(sa + tb, sb - ta);
        }

        private static readonly Complex r31 = new Complex(-1.0 / 2.0, Math.Sqrt(3.0) / 2.0);

        // This length-5 kernel requires just 48 operations, about a third as many as are required by the naive implementation.

        private static void FftKernel5 (Complex[] v, Complex[] y, int y0, int dy, int sign) {
            // first set of combinations
            double a14p = v[1].Re + v[4].Re;
            double a14m = v[1].Re - v[4].Re;
            double a23p = v[2].Re + v[3].Re;
            double a23m = v[2].Re - v[3].Re;
            double b14p = v[1].Im + v[4].Im;
            double b14m = v[1].Im - v[4].Im;
            double b23p = v[2].Im + v[3].Im;
            double b23m = v[2].Im - v[3].Im;
            // second set of combinations, for v[1] and v[4]
            double s14a = v[0].Re + r51.Re * a14p + r52.Re * a23p;
            double s14b = v[0].Im + r51.Re * b14p + r52.Re * b23p;
            double t14a = r51.Im * a14m + r52.Im * a23m;
            double t14b = r51.Im * b14m + r52.Im * b23m;
            // second set of combinations, for v[2] and v[3]
            double s23a = v[0].Re + r52.Re * a14p + r51.Re * a23p;
            double s23b = v[0].Im + r52.Re * b14p + r51.Re * b23p;
            double t23a = r52.Im * a14m - r51.Im * a23m;
            double t23b = r52.Im * b14m - r51.Im * b23m;
            // take care of sign
            if (sign < 0) { t14a = -t14a; t14b = -t14b; t23a = -t23a; t23b = -t23b; }
            // bring together results
            y[y0] = new Complex(v[0].Re + a14p + a23p, v[0].Im + b14p + b23p);
            y[y0 + dy] = new Complex(s14a - t14b, s14b + t14a);
            y[y0 + 2 * dy] = new Complex(s23a - t23b, s23b + t23a);
            y[y0 + 3 * dy] = new Complex(s23a + t23b, s23b - t23a);
            y[y0 + 4 * dy] = new Complex(s14a + t14b, s14b - t14a);
        }

        private static readonly double S5 = Math.Sqrt(5.0);
        private static readonly Complex r51 = new Complex((S5 - 1.0) / 4.0, Math.Sqrt((5.0 + S5) / 8.0));
        private static readonly Complex r52 = new Complex(-(S5 + 1.0) / 4.0, Math.Sqrt((5.0 - S5) / 8.0));

        // this is an O(N^2) kernel which is used when no Winograd kernel is available

        private static void FftKernel (int n, Complex[] v, Complex[] y, int y0, int dy, int sign) {

            // change this to use pre-computed roots
            Complex[] u = ComputeRoots(n, sign);

            int yi = y0;
            y[yi] = 0.0;
            for (int j = 0; j < n; j++) {
                y[yi] += v[j];
            }
            for (int i = 1; i < n; i++) {
                yi += dy;
                y[yi] = v[0];
                int ui = 0;
                for (int j = 1; j < n; j++) {
                    ui += i;
                    if (ui >= n) ui -= n;
                    y[yi] += v[j] * u[ui];
                }
            }

        }

    }

    // factorization logic

    internal struct Factor {

        public Factor (int value, int multiplicity) {
            this.Value = value;
            this.Multiplicity = multiplicity;
        }

        public int Value;
        public int Multiplicity;

        public static List<Factor> Factorize (int n) {
            List<Factor> factors = new List<Factor>();
            foreach (int p in primes) {
                int m = 0;
                while ((n % p) == 0) {
                    m++;
                    n = n / p;
                }
                if (m > 0) factors.Add(new Factor(p, m));
                if (n == 1) return (factors);
            }
            factors.Add(new Factor(n, 1));
            return (factors);
        }

        private static readonly int[] primes = new int[] { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31 };

    }
}
