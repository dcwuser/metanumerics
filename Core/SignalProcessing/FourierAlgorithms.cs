using System;
using System.Collections.Generic;

using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Meta.Numerics.SignalProcessing {

    // Given a length-N DFT, we decompose N = R1 R2 ... Rn into prime factors R. The total DFT can then be expressed as a
    // series of length-R DFTs, where each length-R DFT is repeated N/R times. (Each R need not actually be prime, only
    // co-prime to the other factors.)

    // If a length-N DFT is O(N^2) and N = R1 R2, then a naive implementation would be order N^2 = R1^2 R2^2. But the
    // decomposed work is order (N / R1) R1^2 + (N / R2) R2^2 = R1 R2 (R1 + R2), which is less. We handle large prime
    // factors with the Bluestein algorithm.

    // Each length-R DFT is handled by a transformlet. We have a general transformlet for arbitrary R (the Transformlet class),
    // speicalized dedicated transformlets for small values of R (LengthTwoTransformlet, LengthThreeTransformlet, etc.) and
    // a BluesteinTransformlet for larger R.

    internal class Transformlet {

        // R is the radix, N the total length, and u contains the Nth complex roots of unity

        public Transformlet (int R, int N, Complex[] u) {
            this.R = R;
            this.N = N;
            this.u = u;
            this.dx = N / R;
        }

        protected int R;
        protected int N;
        protected Complex[] u;
        int dx;

        public int Radix {
            get {
                return (R);
            }
        }

        // we don't use the Multiplicity in any Transformlet methods, but it is used to execute the whole plan, and since 
        // there is one per transformlet we save ourselves from creating an additonal class by adding it

        public int Multiplicity { get; internal set; }

        public virtual void FftPass (Complex[] x, Complex[] y, int Ns, int sign) {
            Complex[] v = new Complex[R];
            int dx = N / R;
            for (int j = 0; j < dx; j++) {
                // note: the j-loop could be parallelized, if the v-buffer is not shared between threads
                int xi = j;
                int ui = 0; if (sign < 0) ui = N;
                int du = (dx / Ns) * (j % Ns); if (sign < 0) du = -du;
                // basically, we need to copy x[j + r * dx] * u[r * du] into v[r]
                // we do this in a complicated-looking way in order to avoid unnecessary multiplication when u = 1
                // such a complex multiply requires 6 flops which results in a no-op; we have shown this to be a measurable time-saver
                if (false) {
                    // all u-factors are 1, so we can just reference x directly without copying into v
                    // to do this, we need to change FftKernel to accept an offset and stride for x
                } else {
                    v[0] = x[xi]; // the first u is guaranteed to be 1
                    for (int r = 1; r < R; r++) {
                        xi += dx;
                        ui += du;
                        v[r] = x[xi] * u[ui];
                    }
                }
                int y0 = Expand(j, Ns, R);
                FftKernel(v, y, y0, Ns, sign);
            }
        }

        // This is an O(R^2) DFT. The only thing we have done to speed it up is special-case the u = 1 case to avoid
        // unnecessary complex multiplications. For good performance, override this with a custom kernel for each radix.

        // I am a little worried that this call is virtual. It's not in the innermost loop (which is inside it), but it
        // is in the next loop out. But the whole transformlet architecture gives a significant performance boost over
        // our last architecture, so it's a price I'm willing to pay for now.

        public virtual void FftKernel (Complex[] x, Complex[] y, int y0, int dy, int sign) {

            int yi = y0;
            y[yi] = 0.0;
            for (int j = 0; j < R; j++) {
                y[yi] += x[j];
            }
            for (int i = 1; i < R; i++) {
                yi += dy;
                y[yi] = x[0];
                int ui = 0; if (sign < 0) ui = N;
                int du = dx * i; if (sign < 0) du = -du;
                for (int j = 1; j < R; j++) {
                    ui += du;
                    if (ui >= N) ui -= N; if (ui < 0) ui += N;
                    y[yi] += x[j] * u[ui];
                }
            }

        }

        protected static int Expand (int idxL, int N1, int N2) {
            return ((idxL / N1) * N1 * N2 + (idxL % N1));
        }

    }

    internal class RadixTwoTransformlet : Transformlet {

        public RadixTwoTransformlet (int N, Complex[] u) : base(2, N, u) { }

        public override void FftPass (Complex[] x, Complex[] y, int Ns, int sign) {
            int dx = N / 2;
            for (int j = 0; j < dx; j++) {
                int du = (dx / Ns) * (j % Ns);
                int y0 = Expand(j, Ns, 2);
                if (sign < 0) {
                    FftKernel(x[j], x[j + dx] * u[N - du], out y[y0], out y[y0 + Ns]);
                } else {
                    FftKernel(x[j], x[j + dx] * u[du], out y[y0], out y[y0 + Ns]);
                }
            }
        }

        // performs a length-2 FFT

        private static void FftKernel (Complex x0, Complex x1, out Complex y0, out Complex y1) {
            double a0 = x0.Re; double b0 = x0.Im;
            double a1 = x1.Re; double b1 = x1.Im;
            y0 = new Complex(a0 + a1, b0 + b1);
            y1 = new Complex(a0 - a1, b0 - b1);
            // for some reason, this looks to be faster than using the complex add and subtract; i don't see why
        }

    }

    internal class RadixThreeTransformlet : Transformlet {

        public RadixThreeTransformlet (int N, Complex[] u) : base(3, N, u) { }

        public override void FftPass (Complex[] x, Complex[] y, int Ns, int sign) {
            int dx = N / 3;
            for (int j = 0; j < dx; j++) {
                int du = (dx / Ns) * (j % Ns);
                int y0 = Expand(j, Ns, 3);
                if (sign < 0) {
                    FftKernel(x[j], x[j + dx] * u[N - du], x[j + 2 * dx] * u[N - 2 * du], out y[y0], out y[y0 + Ns], out y[y0 + 2 * Ns] , -1);
                } else {
                    FftKernel(x[j], x[j + dx] * u[du], x[j + 2 * dx] * u[2 * du], out y[y0], out y[y0 + Ns], out y[y0 + 2 * Ns], 1);
                }
            }
        }

        private static void FftKernel (Complex x0, Complex x1, Complex x2, out Complex y0, out Complex y1, out Complex y2, int sign) {
            double a12p = x1.Re + x2.Re;
            double b12p = x1.Im + x2.Im;
            double sa = x0.Re + r31.Re * a12p;
            double sb = x0.Im + r31.Re * b12p;
            double ta = r31.Im * (x1.Re - x2.Re);
            double tb = r31.Im * (x1.Im - x2.Im);
            if (sign < 0) { ta = -ta; tb = -tb; }
            y0 = new Complex(x0.Re + a12p, x0.Im + b12p);
            y1 = new Complex(sa - tb, sb + ta);
            y2 = new Complex(sa + tb, sb - ta);
        }

        private static readonly Complex r31 = new Complex(-1.0 / 2.0, Math.Sqrt(3.0) / 2.0);

    }

    internal class RadixFiveTransformlet : Transformlet {

        public RadixFiveTransformlet (int N, Complex[] u) : base(5, N, u) { }

        public override void FftKernel (Complex[] x, Complex[] y, int y0, int dy, int sign) {
            // first set of combinations
            double a14p = x[1].Re + x[4].Re;
            double a14m = x[1].Re - x[4].Re;
            double a23p = x[2].Re + x[3].Re;
            double a23m = x[2].Re - x[3].Re;
            double b14p = x[1].Im + x[4].Im;
            double b14m = x[1].Im - x[4].Im;
            double b23p = x[2].Im + x[3].Im;
            double b23m = x[2].Im - x[3].Im;
            // second set of combinations, for v[1] and v[4]
            double s14a = x[0].Re + r51.Re * a14p + r52.Re * a23p;
            double s14b = x[0].Im + r51.Re * b14p + r52.Re * b23p;
            double t14a = r51.Im * a14m + r52.Im * a23m;
            double t14b = r51.Im * b14m + r52.Im * b23m;
            // second set of combinations, for v[2] and v[3]
            double s23a = x[0].Re + r52.Re * a14p + r51.Re * a23p;
            double s23b = x[0].Im + r52.Re * b14p + r51.Re * b23p;
            double t23a = r52.Im * a14m - r51.Im * a23m;
            double t23b = r52.Im * b14m - r51.Im * b23m;
            // take care of sign
            if (sign < 0) { t14a = -t14a; t14b = -t14b; t23a = -t23a; t23b = -t23b; }
            // bring together results
            y[y0] = new Complex(x[0].Re + a14p + a23p, x[0].Im + b14p + b23p);
            y[y0 + dy] = new Complex(s14a - t14b, s14b + t14a);
            y[y0 + 2 * dy] = new Complex(s23a - t23b, s23b + t23a);
            y[y0 + 3 * dy] = new Complex(s23a + t23b, s23b - t23a);
            y[y0 + 4 * dy] = new Complex(s14a + t14b, s14b - t14a);
        }

        private static readonly double S5 = Math.Sqrt(5.0);
        private static readonly Complex r51 = new Complex((S5 - 1.0) / 4.0, Math.Sqrt((5.0 + S5) / 8.0));
        private static readonly Complex r52 = new Complex(-(S5 + 1.0) / 4.0, Math.Sqrt((5.0 - S5) / 8.0));

    }

    internal class RadixSevenTransformlet : Transformlet {

        public RadixSevenTransformlet (int N, Complex[] u) : base(7, N, u) { }

        public override void FftKernel (Complex[] x, Complex[] y, int y0, int dy, int sign) {
            // relevent sums and differences
            double a16p = x[1].Re + x[6].Re;
            double a16m = x[1].Re - x[6].Re;
            double a25p = x[2].Re + x[5].Re;
            double a25m = x[2].Re - x[5].Re;
            double a34p = x[3].Re + x[4].Re;
            double a34m = x[3].Re - x[4].Re;
            double b16p = x[1].Im + x[6].Im;
            double b16m = x[1].Im - x[6].Im;
            double b25p = x[2].Im + x[5].Im;
            double b25m = x[2].Im - x[5].Im;
            double b34p = x[3].Im + x[4].Im;
            double b34m = x[3].Im - x[4].Im;
            // combinations used in y[1] and y[6]
            double s16a = x[0].Re + r71.Re * a16p + r72.Re * a25p + r73.Re * a34p;
            double s16b = x[0].Im + r71.Re * b16p + r72.Re * b25p + r73.Re * b34p;
            double t16a = r71.Im * a16m + r72.Im * a25m + r73.Im * a34m;
            double t16b = r71.Im * b16m + r72.Im * b25m + r73.Im * b34m;
            // combinations used in y[2] and y[5]
            double s25a = x[0].Re + r71.Re * a34p + r72.Re * a16p + r73.Re * a25p;
            double s25b = x[0].Im + r71.Re * b34p + r72.Re * b16p + r73.Re * b25p;
            double t25a = r71.Im * a34m - r72.Im * a16m + r73.Im * a25m;
            double t25b = r71.Im * b34m - r72.Im * b16m + r73.Im * b25m;
            // combinations used in y[3] and y[4]
            double s34a = x[0].Re + r71.Re * a25p + r72.Re * a34p + r73.Re * a16p;
            double s34b = x[0].Im + r71.Re * b25p + r72.Re * b34p + r73.Re * b16p;
            double t34a = r71.Im * a25m - r72.Im * a34m - r73.Im * a16m;
            double t34b = r71.Im * b25m - r72.Im * b34m - r73.Im * b16m;
            // if sign is negative, invert t's
            if (sign < 0) {
                t16a = -t16a; t16b = -t16b;
                t25a = -t25a; t25b = -t25b;
                t34a = -t34a; t34b = -t34b;
            }
            // combine to get results
            y[y0] = new Complex(x[0].Re + a16p + a25p + a34p, x[0].Im + b16p + b25p + b34p);
            y[y0 + dy] = new Complex(s16a - t16b, s16b + t16a);
            y[y0 + 2 * dy] = new Complex(s25a + t25b, s25b - t25a);
            y[y0 + 3 * dy] = new Complex(s34a + t34b, s34b - t34a);
            y[y0 + 4 * dy] = new Complex(s34a - t34b, s34b + t34a);
            y[y0 + 5 * dy] = new Complex(s25a - t25b, s25b + t25a);
            y[y0 + 6 * dy] = new Complex(s16a + t16b, s16b - t16a);
        }

        // seventh roots of unity
        // a la Gauss, these are not expressible in closed form using rationals and rational roots

        private static readonly Complex r71 = new Complex(0.62348980185873353053, 0.78183148246802980871);
        private static readonly Complex r72 = new Complex(-0.22252093395631440429, 0.97492791218182360702);
        private static readonly Complex r73 = new Complex(-0.90096886790241912624, 0.43388373911755812048);

    }


    // The Bluestein technique works as follows. Given the length-N FT
    //  \tilde{x}_m = \sum_{n=0}^{N-1} x_n \exp{i \pm 2 \pi m n / N}
    // use m n = \frac{m^2 + n^2 - (m - n)^2}{2} to turn this into
    //   \tilde{x}_m = \exp{i \pm \pi m^2 / N} \sum_{n=0}^{N-1} x_n \exp{i \pm \pi n^2 / N} \exp{i \mp \pi (m - n)^2 / N}
    // The summed expression is a convolution of
    //   a_n = x_n \exp{i \pm \pi n^2 / N}
    //   b_n = \exp{i \mp \pi n^2 / N}
    // A convolution can be done via an FT of any length larger than 2N-1. The 2N is necessary so that a_0 can be multiplied
    // by b_{-N} and a_N can be multiplied by b_0. This thus the sequences to be convolved are
    //   0 0  0       0   0  a_0 a_1 a_2 ... a_n 0 0 0
    //   0 0 b_n ... b_2 b_1 b_0 b_1 b_2 ... b_n 0 0 0
    // Since this is a convolution, it doesn't matter how far out we zero-pad. We pick an M >= 2N-1 that is composed of
    // small prime factors, so we won't need the Bluestein technique to do the convolution itself.

    internal class BluesteinTransformlet : Transformlet {

        public BluesteinTransformlet (int R, int N, Complex[] u) : base(R, N, u) {

            // figure out the right Bluestein length and create a transformer for it
            Nb = SetBluesteinLength(2 * R - 1);
            ft = new FourierTransformer(Nb);

            // compute the Bluestein coeffcients and compute the FT of the filter based on them
            b = ComputeBluesteinCoefficients(R);
            Complex[] c = new Complex[Nb];
            c[0] = 1.0;
            for (int i = 1; i < R; i++) {
                c[i] = b[i].Conjugate;
                c[Nb - i] = c[i];
            }
            bt = ft.Transform(c);

        }

        // Length of convolution transform
        private int Nb;

        // Fourier transform for convolution transform
        private FourierTransformer ft;

        // R Bluestein coefficients
        private Complex[] b;

        // Nb-length Fourier transform of the symmetric Bluestein coefficient filter
        private Complex[] bt;

        // This method computes b_n = \exp{i \pi n^2 / N}. If we do this naively, by computing sin and cos of \pi n^2 / N, then
        // the argument can get large, up to N \pi, and the inaccuracy of trig methods for large arguments will hit us
        // To avoid this, note that the difference n^2 - (n-1)^2 = 2n-1. So we can add 2n-1 each time and take the result mod 2N
        // to keep the argument less than 2 \pi.

        private static Complex[] ComputeBluesteinCoefficients (int R) {

            Complex[] b = new Complex[R];
            double t = Math.PI / R;
            b[0] = 1.0;
            int s = 0;
            int TwoR = 2 * R;
            for (int i = 1; i < R; i++) {
                s += (2 * i - 1); if (s >= TwoR) s -= TwoR;
                double ts = t * s;
                b[i] = new Complex(Math.Cos(ts), Math.Sin(ts));
            }
            return (b);
        }

        public override void FftKernel (Complex[] x, Complex[] y, int y0, int dy, int sign) {

            // all we have to do here is convolve (b x) with b-star
            // to do this convolution, we need to multiply the DFT of (b x) with the DFT of b-star, the IDFT the result back
            // we have already stored the DFT of b-star

            // create c = b x and transform it into Fourier space
            Complex[] c = new Complex[Nb];
            if (sign > 0) {
                for (int i = 0; i < R; i++) c[i] = b[i] * x[i];
            } else {
                for (int i = 0; i < R; i++) c[i] = b[i] * x[i].Conjugate;
            }
            Complex[] ct = ft.Transform(c);

            // multiply b-star and c = b x in Fourier space, and inverse transform the product back into configuration space
            for (int i = 0; i < Nb; i++) {
                ct[i] = ct[i] * bt[i];
            }
            c = ft.InverseTransform(ct);

            // read off the result
            if (sign > 0) {
                for (int i = 0; i < R; i++) y[y0 + i * dy] = b[i] * c[i];
            } else {
                for (int i = 0; i < R; i++) y[y0 + i * dy] = b[i].Conjugate * c[i].Conjugate;
            }

            // for the sign < 0 case, we have used the fact that the convolution of (b-star x) with b is
            // just the convolution of (b x-star) with b-star, starred

        }

        // This is all about determining a good value to use for the Bluestein length. We choose a length based on powers
        // of two and three, since those give very fast Fourier transforms. Our method is heuristic and not optomized.

        private static int SetBluesteinLength (int N) {

            // try the next power of two
            int t = NextPowerOfBase(N, 2);
            int M = t;

            // see if replacing factors of 4 by 3, which shortens the length, will still be long enough
            while (M % 4 == 0) {
                t = (M / 4) * 3;
                if (t < N) break;
                if (t < M) M = t;
            }

            // try the next power of three
            t = NextPowerOfBase(N, 3);
            if ((t > 0) && (t < M)) M = t;

            return (M);

        }

        private static int NextPowerOfBase (int n, int b) {
            for (int m = b; m <= Int32.MaxValue; m *= b) {
                if (m >= n) return (m);
            }
            return (-1);
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

    }

}
