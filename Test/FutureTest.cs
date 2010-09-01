using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;

namespace Test {

    [TestClass()]
    public class FutureTest {

        private TestContext testContextInstance;

        /// <summary>
        ///Gets or sets the test context which provides
        ///information about and functionality for the current test run.
        ///</summary>
        public TestContext TestContext {
            get {
                return testContextInstance;
            }
            set {
                testContextInstance = value;
            }
        }


        public static void DFT2 (Complex[] x, int[] p) {

            Complex f0 = x[p[0]];
            Complex f1 = x[p[1]];
            x[p[0]] = f0 + f1;
            x[p[1]] = f0 - f1;

        }

        public static void DFT3 (Complex[] x, int[] p) {

            Complex f0 = x[p[0]];
            Complex f1 = x[p[1]];
            Complex f2 = x[p[2]];

            x[p[0]] = f0 + f1 + f2;
            x[p[1]] = f0 + W31 * f1 + W32 * f2;
            x[p[2]] = f0 + W32 * f1 + W31 * f2;

        }

        private static readonly Complex W31 = new Complex(-1.0 / 2.0, -Math.Sqrt(3.0) / 2.0);
        private static readonly Complex W32 = new Complex(-1.0 / 2.0, +Math.Sqrt(3.0) / 2.0);

        public static void DFT4 (Complex[] x, int[] p) {

            Complex f0 = x[p[0]];
            Complex f1 = x[p[1]];
            Complex f2 = x[p[2]];
            Complex f3 = x[p[3]];

            Complex f0p2 = f0 + f2;
            Complex f0m2 = f0 - f2;
            Complex f1p3 = f1 + f3;
            Complex f1m3 = f1 - f3;

            x[p[0]] = f0p2 + f1p3;
            x[p[1]] = f0m2 - ComplexMath.I * f1m3;
            x[p[2]] = f0p2 - f1p3;
            x[p[3]] = f0m2 + ComplexMath.I * f1m3;

        }

        public static void ProductDFT (Complex[] x, int n1, int n2) {

            int n = x.Length;
            if (n1 * n2 != n) throw new InvalidOperationException();

            for (int k1 = 0; k1 < n1; k1++) {
                // set up mapping for n2-length transform
                int[] p = new int[n2];
                for (int k2 = 0; k2 < n2; k2++) {
                    int k = (n2 * k1 + n1 * k2) % n;
                    p[k2] = k;
                }
                // do n2-length transform
                switch (n2) {
                    case 2:
                        DFT2(x, p);
                        break;
                    case 3:
                        DFT3(x, p);
                        break;
                    case 4:
                        DFT4(x, p);
                        break;
                    default:
                        throw new NotImplementedException();
                }
            }

            for (int k2 = 0; k2 < n2; k2++) {
                // set up mapping for n1-length transform
                int[] p = new int[n1];
                for (int k1 = 0; k1 < n1; k1++) {
                    int k = (n2 * k1 + n1 * k2) % n;
                    p[k1] = k;
                }
                // do n1-length transform
                switch (n1) {
                    case 2:
                        DFT2(x, p);
                        break;
                    case 3:
                        DFT3(x, p);
                        break;
                    case 4:
                        DFT4(x, p);
                        break;
                    default:
                        throw new NotImplementedException();
                }
            }

        }

        public static Complex[] DumbDFT (Complex[] x) {

            // extract length
            int N = x.Length;

            // compute twidle factors
            double t = 2.0 * Math.PI / N;
            Complex[] W = new Complex[N];
            W[0] = 1.0;
            W[1] = new Complex(Math.Cos(t), -Math.Sin(t));
            for (int iw = 2; iw < N; iw++) {
                W[iw] = W[iw - 1] * W[1];
            }

            // compute each frequency component
            Complex[] y = new Complex[N];
            for (int iy = 0; iy < N; iy++) {
                Complex y0 = 0.0;
                int iw = 0;
                for (int ix = 0; ix < N; ix++) {
                    y0 = y0 + x[ix] * W[iw];
                    iw = (iw + iy) % N;
                }
                y[iy] = y0;
            }

            return (y);
        }

        //public static Complex[] DFT (Complex[] x, int[] F) {

        //    int N = x.Length;

        //    // check length
        //    int NN = 1;
        //    for (int k = 0; k < F.Length; k++) {
        //        NN = NN * F[k];
        //    }
        //    if (NN != N) throw new InvalidOperationException();

        //    Complex[] y = new Complex[x.Length];

        //    // loop over factors
        //    for (int k = 0; k < F.Length; k++) {
        //        int N1 = F[k];
        //        int N2 = N / F[k];

        //        int L = 1;
        //        int N3 = N2 - N1 * (N2 / N1);

        //        int[] LP;
        //        for (int j = 0; j < N1; j++) {
        //            LP[j] = L;
        //            L = L + N3;
        //            if (L >= N1) L = L - N1;
        //        }

        //        int[] I, IP;
        //        for (int j = 0; j < N; j += N1) {

        //            int it = j;
        //            for (int l = 0; l < N1; l++) {
        //                I[L] = it;
        //                IP[LP[L]] = it;
        //                it = it + N2;
        //                if (it >= N) it = it - N;
        //            }

        //            if (N1 == 2) {
        //                Complex r1 = x[I[1]];
        //                x[I[1]] = r1 + x[I[2]];
        //                x[I[2]] = r1 - x[I[2]];
        //            } else {
        //                throw new NotImplementedException();
        //            }
        //            /*
        //            // populate y from x
        //            int i1 = j;
        //            for (int l = 0; l < N1; l++) {
        //                y[l] = x[i1];
        //                i1 = (i1 + N2) % N;
        //            }
        //            // do explicit FFT of order N1
        //            if (N1 == 2) {
        //                Complex y0 = y[0];
        //                Complex y1 = y[1];
        //                y[0] = y0 + y1;
        //                y[1] = y0 - y1;
        //            } else {
        //                throw new NotImplementedException();
        //            }
        //            // populate x from y
        //            int i2 = j;
        //            for (int l = 0; l < N1; l++) {
        //                x[i2] = y[l];
        //                i2 = (i2 + N2) % N;
        //            }
        //            */
        //        }
        //    }
        //    // unscramble
        //    return (x);

        //}

        [TestMethod]
        public void TestDFT () {

            /*
            Complex[] x = new Complex[] {
                new Complex(1.0,0.0), new Complex(1.0,0.0)
            };
            int[] F = new int[] { 2 };
            */


            Complex[] x1 = new Complex[] {
                new Complex(1.0,0.0), new Complex(2.0,0.0), new Complex(3.0,0.0),
                new Complex(4.0,0.0), new Complex(5.0,0.0), new Complex(6.0,0.0)
            };
            Complex[] x2 = new Complex[x1.Length];
            Array.Copy(x1, x2, x1.Length);

            Complex[] y1 = DumbDFT(x1);
            for (int i = 0; i < y1.Length; i++) {
                Console.WriteLine("{0} {1}", i, y1[i]);
            }

            //DFT4(x2, new int[] { 0, 1, 2, 3 });
            ProductDFT(x2, 2, 3);
            for (int i = 0; i < x2.Length; i++) {
                Console.WriteLine("{0} {1}", i, x2[i]);
            }


        }

    }
}
