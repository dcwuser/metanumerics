using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Statistics;

namespace FutureTest {

    // proposed matrix class hierarchy

    public abstract class RectangularMatrixBase {

        public abstract double this[int r, int c] { get; set; }

        public abstract int RowCount { get; }

        public abstract int ColumnCount { get; }

        // we an implement some operations, but many will be slow because they do not have access to the underlying storage
        // override them with faster implementations

        public virtual double OneNorm () {
            throw new NotImplementedException();
        }

        public virtual double InfinityNorm () {
            throw new NotImplementedException();
        }

        public static RectangularMatrix operator + (RectangularMatrixBase M1, RectangularMatrixBase M2) {
            if (M1 == null) throw new ArgumentNullException("M1");
            if (M2 == null) throw new ArgumentNullException("M2");
            throw new NotImplementedException();
        }

        public static RectangularMatrix operator * (RectangularMatrixBase M1, RectangularMatrixBase M2) {
            throw new NotImplementedException();
        }

        public static bool operator == (RectangularMatrixBase M1, RectangularMatrixBase M2) {
            throw new NotImplementedException();
        }

        public static bool operator != (RectangularMatrixBase M1, RectangularMatrixBase M2) {
            return (!(M1 == M2));
        }

        public override bool Equals (object obj) {
            RectangularMatrixBase M = obj as RectangularMatrix;
            if (obj == null) {
                return (false);
            } else {
                return ((this == M));
            }
        }

        public override int GetHashCode () {
            throw new NotImplementedException();
        }
    }

    public sealed class RectangularMatrix : RectangularMatrixBase {

        public override double this[int r, int c] {
            get {
                throw new NotImplementedException();
            }
            set {
                throw new NotImplementedException();
            }
        }

        public override int RowCount {
            get { throw new NotImplementedException(); }
        }

        public override int ColumnCount {
            get { throw new NotImplementedException(); }
        }

        public static RectangularMatrix operator * (RectangularMatrix M1, RectangularMatrix M2) {
            // this is faster than the base operator, because it knows about the underlying structure
            throw new NotImplementedException();
        }

        public static ColumnVector operator * (RectangularMatrix M1, ColumnVector v1) {
            // this is faster than the base operator, because it knows about the underlying structure
            throw new NotImplementedException();
        }

    }

    public abstract class SquareMatrixBase : RectangularMatrixBase {

        public abstract int Dimension { get; }

        public override int RowCount {
            get { return (Dimension); }
        }

        public override int ColumnCount {
            get { return (Dimension); }
        }

        public virtual double Trace () {
            throw new NotImplementedException();
        }

        public static SquareMatrix operator + (SquareMatrixBase M1, SquareMatrixBase M2) {
            throw new NotImplementedException();
        }

        public static SquareMatrix operator * (SquareMatrixBase M1, SquareMatrixBase M2) {
            throw new NotImplementedException();
        }

    }

    public sealed class SquareMatrix : SquareMatrixBase {

        public override double this[int r, int c] {
            get {
                throw new NotImplementedException();
            }
            set {
                throw new NotImplementedException();
            }
        }

        public override int Dimension {
            get { throw new NotImplementedException(); }
        }

        public SquareMatrix Inverse () {
            throw new NotImplementedException();
        }

    }

    public sealed class SymmetricMatrix : SquareMatrixBase {

        public override double this[int r, int c] {
            get {
                throw new NotImplementedException();
            }
            set {
                throw new NotImplementedException();
            }
        }

        public override int Dimension {
            get { throw new NotImplementedException(); }
        }

        public SymmetricMatrix Inverse () {
            throw new NotImplementedException();
        }

        public static SymmetricMatrix operator + (SymmetricMatrix M1, SymmetricMatrix M2) {
            throw new NotImplementedException();
        }

    }

    public sealed class TridiagonalMatrix : SquareMatrixBase {

        public override double this[int r, int c] {
            get {
                throw new NotImplementedException();
            }
            set {
                throw new NotImplementedException();
            }
        }

        public override int Dimension {
            get { throw new NotImplementedException(); }
        }

        public SquareMatrix Inverse () {
            throw new NotImplementedException();
        }

        public static TridiagonalMatrix operator + (TridiagonalMatrix M1, TridiagonalMatrix M2) {
            throw new NotImplementedException();
        }

        public static SquareMatrix operator * (TridiagonalMatrix M1, TridiagonalMatrix M2) {
            throw new NotImplementedException();
        }

    }

    public abstract class VectorBase : RectangularMatrixBase {

        public virtual double this[int n] {
            get {
                throw new NotImplementedException();
            }
            set {
                throw new NotImplementedException();
            }
        }

        public virtual int Dimension {
            get {
                throw new NotImplementedException();
            }
        }

        public override int RowCount {
            get { return (Dimension); }
        }

        public override int ColumnCount {
            get { return (Dimension); }
        }

        public virtual double TwoNorm () {
            throw new NotImplementedException();
        }

    }

    public sealed class ColumnVector : VectorBase {

        public override double this[int r, int c] {
            get { throw new NotImplementedException(); }
            set { throw new NotImplementedException(); }
        }

        public static ColumnVector operator + (ColumnVector v1, ColumnVector v2) {
            throw new NotImplementedException();
        }

        public static ColumnVector operator * (RectangularMatrixBase M1, ColumnVector v2) {
            throw new NotImplementedException();
        }

    }

    public sealed class RowVector : VectorBase {

        public override double this[int r, int c] {
            get { throw new NotImplementedException(); }
            set { throw new NotImplementedException(); }
        }

        public static RowVector operator + (RowVector v1, RowVector v2) {
            throw new NotImplementedException();
        }

        public static RowVector operator * (RowVector v1, RectangularMatrixBase M2) {
            throw new NotImplementedException();
        }

        public static double operator * (RowVector v1, ColumnVector v2) {
            throw new NotImplementedException();
        }

    }



    [TestClass()]
    public class FutureTest {

        //[TestMethod]
        public void MatrixOps () {
            RectangularMatrix RT = new RectangularMatrix();
            SquareMatrix SQ = new SquareMatrix();
            SymmetricMatrix SY = new SymmetricMatrix();
            TridiagonalMatrix TR = new TridiagonalMatrix();

            ColumnVector CV = new ColumnVector();
            RowVector RV = new RowVector();

            RectangularMatrixBase R;
            R = RT * RT;
            R = RT * SQ;
            R = RT * SY;
            R = RT * TR;
            R = SQ * RT;
            R = SQ * SQ;
            R = SQ * SY;
            R = SQ * TR;
            R = SY * RT;
            R = SY * SQ;
            R = SY * SY;
            R = SY * TR;
            R = TR * RT;
            R = TR * SQ;
            R = TR * SY;
            R = TR * TR;

            ColumnVector V;
            V = RT * CV;
            V = SQ * CV;
            V = SY * CV;
            V = TR * CV;

        }

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

        public static void StockhamFFTIteration (Complex[] x, Complex[] y, int j, int N, int R, int Ns) {
            StockhamFFTIteration(x, y, j, N, R, Ns, 1);
        }

        public static void StockhamFFTIteration (Complex[] x, Complex[] y, int j, int N, int R, int Ns, int q) {
            Complex[] v = new Complex[R];
            double angle = -2.0 * Math.PI * (j % Ns) / (Ns * R);
            for (int r = 0; r < R; r++) {
                Complex tw = new Complex(Math.Cos(r * angle), Math.Sin(r * angle));
                Console.WriteLine("R={0} Ns={1} j={2} r={3} tw={4}", R, Ns, j, r, tw);
                Console.WriteLine("{0} -> {1}", j + r * N / R, ((j / Ns) * Ns * R + (j % Ns) + r * Ns * q) % N);
                v[r] = x[j + r * N / R] * tw;
            }
            if (R == 2) {
                Complex v0 = v[0];
                v[0] = v0 + v[1];
                v[1] = v0 - v[1];
            } else if (R == 3) {
                Complex t0 = v[0] + v[1] + v[2];
                Complex t1 = v[0] + R3.Conjugate * v[1] + R3 * v[2];
                Complex t2 = v[0] + R3 * v[1] + R3.Conjugate * v[2];
                v[0] = t0; v[1] = t1; v[2] = t2;
            } else {
                throw new NotImplementedException();
            }
            int idxD = (j / Ns) * Ns * R + (j % Ns);
            for (int r = 0; r < R; r++) {
                y[(idxD + r * Ns * q) % N] = v[r];
            }
        }

        public static void FFT (Complex[] v, int p) {
            switch (p) {
                case 2:
                    break;
                case 3:
                    break;
                default:
                    throw new NotImplementedException();
            }
        }

        public static void StockhamFFT (Complex[] x, Complex[] y, int N, int R) {
            int p = 1;
            for (int Ns = 1; Ns < N; Ns *= R) {
                for (int j = 0; j < N / R; j++) {
                    StockhamFFTIteration(x, y, j, N, R, Ns);
                }
                Complex[] t = x; x = y; y = t; p *= -1;
            }
            if (p < 0) {
                for (int i = 0; i < x.Length; i++) {
                    y[i] = x[i];
                }
            }
        }

        // Basic implementation

        public static void FftIteration (int j, int N, int R, int Ns, Complex[] data0, Complex[] data1, Complex[] factors) {
            Complex[] v = new Complex[R];
            int idxS = j;
            double angle = -2.0 * Math.PI * (j % Ns) / (Ns * R);
            int step = (N / Ns / R) * (j % Ns);
            for (int r = 0; r < R; r++) {
                if (factors == null) {
                    v[r] = data0[idxS + r * N / R] * new Complex(Math.Cos(r * angle), Math.Sin(r * angle));
                } else {
                    v[r] = data0[idxS + r * N / R] * factors[r * step];
                }
                //if (r == 1) Console.WriteLine("{0} {1} {2} {3} {4}", N, Ns, R, j, angle / (-2.0 * Math.PI));
            }
            Fft(R, v, -1);
            int idxD = Expand(j, Ns, R);
            for (int r = 0; r < R; r++) {
                //Console.WriteLine("j={0} r={1} {2} -> {3}", j, r, idxS + r * N / R, idxD + r * Ns); 
                data1[idxD + r * Ns] = v[r];
            }
        }

        public static void FftKernel2 (Complex[] x, int x0, int dx, Complex[] y, int y0, int dy) {
            int x1 = x0 + dx;
            int y1 = y0 + dy;
            y[y0] = x[x0] + x[x1];
            y[y1] = x[x0] - x[x1];
        }

        // This length-3 FFT kernel requires just 12 additions and 4 multiplications, or 16 operations total.
        // This is half as many operations as is required by the naive implementation. 

        public static void Fft3 (Complex[] v, int sign) {
            double a12p = v[1].Re + v[2].Re;
            double b12p = v[1].Im + v[2].Im;
            double sa = v[0].Re + r31.Re * a12p;
            double sb = v[0].Im + r31.Re * b12p;
            double ta = r31.Im * (v[1].Re - v[2].Re);
            double tb = r31.Im * (v[1].Im - v[2].Im);
            if (sign < 0) { ta = -ta; tb = -tb; }
            v[0] = new Complex(v[0].Re + a12p, v[0].Im + b12p);
            v[1] = new Complex(sa - tb, sb + ta);
            v[2] = new Complex(sa + tb, sb - ta);
        }

        private static readonly Complex r31 = new Complex(-1.0 / 2.0, Math.Sqrt(3.0) / 2.0);

        // This length-5 kernel requires just 48 operations, about a third as many as are required by the naive implementation.

        public static void Fft5 (Complex[] v, int sign) {
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
            v[0] = new Complex(v[0].Re + a14p + a23p, v[0].Im + b14p + b23p);
            v[1] = new Complex(s14a - t14b, s14b + t14a);
            v[2] = new Complex(s23a - t23b, s23b + t23a);
            v[3] = new Complex(s23a + t23b, s23b - t23a);
            v[4] = new Complex(s14a + t14b, s14b - t14a);
        }

        private static readonly double S5 = Math.Sqrt(5.0);
        private static readonly Complex r51 = new Complex((S5 - 1.0) / 4.0, Math.Sqrt((5.0 + S5) / 8.0));
        private static readonly Complex r52 = new Complex(-(S5 + 1.0) / 4.0, Math.Sqrt((5.0 - S5) / 8.0));

        public static void Fft (Complex[] v, int sign) {
            Complex[] w = new Complex[v.Length];
            Array.Copy(v, w, v.Length);
            Complex[] u = ComputePowers(v.Length, sign);
            for (int i = 0; i < v.Length; i++) {
                v[i] = 0.0;
                int c = 0;
                for (int j = 0; j < w.Length; j++) {
                    v[i] += u[c] * w[j];
                    c += i;
                    if (c >= v.Length) c -= v.Length;
                }
            }
        }

        public static void Fft (int R, Complex[] v, int sign) {
            if (R == 2) {
                double a0 = v[0].Re; double b0 = v[0].Im;
                double a1 = v[1].Re; double b1 = v[1].Im;
                v[0] = new Complex(a0 + a1, b0 + b1);
                v[1] = new Complex(a0 - a1, b0 - b1);
            } else if (R == 3) {
                Fft3(v, sign);
            } else if (R == 5) {
                Fft5(v, sign);
            } else {
                Fft(v, sign);
            }
        }

        private static readonly Complex R3 = new Complex(-1.0 / 2.0, Math.Sqrt(3.0) / 2.0);

        public static int Expand (int idxL, int N1, int N2) {
            return ((idxL / N1) * N1 * N2 + (idxL % N1));
        }

        // Radix-R entry point

        public static Complex[] Fft (int N, int R, Complex[] data0) {
            Complex[] data1 = new Complex[data0.Length];
            for (int Ns = 1; Ns < N; Ns *= R) {
                for (int j = 0; j < N / R; j++) {
                    FftIteration(j, N, R, Ns, data0, data1, null);
                }
                Complex[] temp = data0; data0 = data1; data1 = temp;
            }
            return (data0);
        }

        // Mixed-radix entry point

        public static Complex[] Fft (int N, int R, int p, Complex[] data0) {
            Complex[] data1 = new Complex[data0.Length];
            int Ns = 1;
            for (int t = 0; t < p; t++) {
                for (int j = 0; j < N / R; j++) {
                    FftIteration(j, N, R, Ns, data0, data1, null);
                }
                Ns *= R;
                Complex[] temp = data0; data0 = data1; data1 = temp;
            }
            return (data0);
        }

        internal static Complex[] Fft (int N, List<Factor> factors, Complex[] x, int sign) {
            Complex[] powers = ComputePowers(N, sign);
            Complex[] y = new Complex[x.Length];
            int Ns = 1;
            foreach (Factor factor in factors) {
                for (int t = 0; t < factor.Multiplicity; t++) {
                    FftCore(N, factor.Value, Ns, x, y, powers, sign);
                    Ns *= factor.Value;
                    Complex[] temp = x; x = y; y = temp;
                }
            }
            return (x);
        }

        public static void FftCore (int N, int R, int Ns, Complex[] data0, Complex[] data1, Complex[] factors, int sign) {
            Complex[] v = new Complex[R];
            for (int j = 0; j < N / R; j++) {
                // the j-loop can be parallelized (if the v buffer is not shared)
                int idxS = j;
                int dStep = (N / R);
                int fStep = (N / Ns / R) * (j % Ns);
                v[0] = data0[idxS];
                for (int r = 1; r < R; r++) {
                    v[r] = data0[idxS + r * dStep] * factors[r * fStep];
                }
                Fft(R, v, sign);
                int idxD = Expand(j, Ns, R);
                for (int r = 0; r < R; r++) {
                    data1[idxD + r * Ns] = v[r];
                }
            }
        }

        public Complex[] ComputeFactors (int N) {
            Complex[] u = new Complex[N];
            double t = -2.0 * Math.PI / N;
            for (int r = 0; r < N; r++) {
                double rt = r * t;
                u[r] = new Complex(Math.Cos(rt), Math.Sin(rt));
            }
            return (u);
        }

        // Test director method

        [TestMethod]
        public void TestFft () {

            /*

            //Complex[] x = { 0.0, 1.0 };
            //Complex[] x = { 0.0, 0.0, 0.0, 1.0 };
            //Complex[] x = { 0.0, 0.0, 1.0 };
            Complex[] x = { 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

            Complex[] z = Fft(9, 3, x);

            for (int i = 0; i < z.Length; i++) {
                Console.WriteLine("{0}", z[i]);
            }
            */

            //Complex[] x = { 0.0, 1.0, 0.0, 0.0, 0.0 };
            //Complex[] x = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 };
            //Complex[] x = { 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

            Complex[] x = new Complex[25];
            x[1] = 1.0;

            Stopwatch sw = Stopwatch.StartNew();
            //Complex[] u = ComputeFactors(x.Length);
            List<Factor> factors = Factor.Factorize(x.Length);
            Complex[] y = Fft(x.Length, factors, x, -1);
            sw.Stop();
            Console.WriteLine("{0} ms", sw.ElapsedMilliseconds);
            /*
            Complex[] y = new Complex[x.Length];
            for (int j = 0; j < 6; j++) {
                FftIteration(j, 12, 2, 1, x, y);
            }
            for (int j = 0; j < 6; j++) {
                FftIteration(j, 12, 2, 2, y, x);
            }
            for (int i = 0; i < y.Length; i++) {
                Console.WriteLine(y[i]);
            }

            for (int j = 0; j < 4; j++) {
                FftIteration(j, 12, 3, 4, x, y);
            }
            */
            
            Console.WriteLine(x.Length);
            for (int i = 0; i < x.Length; i++) {
                Console.WriteLine(y[i]);
            }
            Console.WriteLine("-");
            Complex[] z = Fft(y.Length, factors, y, +1);
            for (int i = 0; i < x.Length; i++) {
                Console.WriteLine(z[i]);
            }
            Console.WriteLine("-");
            double theta = - 2.0 * Math.PI / x.Length;
            Complex z1 = new Complex(Math.Cos(theta), Math.Sin(theta));
            Console.WriteLine("{0} {1}", z1, z1 * z1);
            
            
            /*
            Complex[] x = { 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

            
            //Complex[] y = Fft(6, 2, 1, x);
            Complex[] z = Fft(9, 3, 2, x);

            for (int i = 0; i < z.Length; i++) {
                Console.WriteLine(z[i]);
            }
            */

        }

        public static void StockhamFFT (Complex[] x, Complex[] y, int N, int R, int p, int q) {
            int pi = 1;
            int Ns = 1;
            for (int pt = 0; pt < p; pt++) {
                for (int j = 0; j < N / R; j++) {
                    StockhamFFTIteration(x, y, j, N, R, Ns, q);
                }
                Complex[] t = x; x = y; y = t; pi *= -1;
                Ns *= R;
            }
            if (pi < 0) {
                for (int i = 0; i < x.Length; i++) {
                    y[i] = x[i];
                }
            }
        }

        public static Complex[] ComputePowers (int N, int sign) {
            Complex[] u = new Complex[N];
            double alpha = 2.0 * Math.PI / N;
            if (sign < 0) alpha = -alpha;
            for (int i = 0; i < N; i++) {
                u[i] = new Complex(Math.Cos(i * alpha), Math.Sin(i * alpha));
            }
            return (u);
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
        public void TestSFFT () {

            // Powers of 2
            //Complex[] x = { 0.0, 1.0 };
            //Complex[] x = { 1.0, 1.0, 1.0, 1.0 };
            //Complex[] x = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
            //Complex[] x = { 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            //Complex[] x = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
            //Complex[] x = { 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0 };
            // Powers of 3
            //Complex[] x = { 0.0, 0.0, 1.0 };
            //Complex[] x = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
            // 6 = 2 X 3
            Complex[] x = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            // 12 = 2^2 X 3
            //Complex[] x = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

            Complex[] y = new Complex[x.Length];
            //StockhamFFT(x, y, x.Length, 2);

            StockhamFFT(x, y, x.Length, 2, 1, 1);
            StockhamFFT(x, y, x.Length, 3, 1, 1);

            //StockhamFFT(x, y, x.Length, 3, 1, 1.0);
            //StockhamFFT(x, y, x.Length, 2, 1, ComplexMath.Exp(-2.0 * Math.PI * ComplexMath.I * 3.0 / 6.0));

            for (int i = 0; i < x.Length; i++) {
                Console.WriteLine("{0}    {1}", x[i], y[i]);
            }

        }

        //[TestMethod]
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

        [TestMethod]
        public void TestFactorize () {

            int nm = 0;
            int ns = 0;
            for (int n = 2; n < 64; n++) {

                List<Factor> factors = Factor.Factorize(n);
                foreach (Factor factor in factors) {
                    Console.Write("{0}^{1} ", factor.Value, factor.Multiplicity);
                }
                Console.WriteLine(" = {0}", n);
                if (factors[factors.Count - 1].Value > 16) ns++;
                nm++;

            }

            Console.WriteLine(((double) ns) / nm);

        }

    }

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
