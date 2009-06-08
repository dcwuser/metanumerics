using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Matrices {

#if FUTURE

    /// <summary>
    /// Represents a tridiagonal matrix.
    /// </summary>
    public sealed class TridiagonalMatrix : IMatrix {

        /// <summary>
        /// Initializes a new tridiagonal matrix of the given dimension.
        /// </summary>
        /// <param name="dimension"></param>
        public TridiagonalMatrix (int dimension) {
            if (dimension < 1) throw new ArgumentOutOfRangeException("dimension");
            this.dimension = dimension;
            U = new double[dimension - 1];
            D = new double[dimension];
            L = new double[dimension - 1];
        }

        int dimension;
        private double[] U;
        private double[] D;
        private double[] L;

        /// <summary>
        /// Gets the dimension of the matrix.
        /// </summary>
        public int Dimension {
            get {
                return (dimension);
            }
        }

        int IMatrix.RowCount {
            get {
                return (dimension);
            }
        }

        int IMatrix.ColumnCount {
            get {
                return (dimension);
            }
        }

        internal double GetDiagonalElement (int i) {
            return (D[i]);
        }

        internal double GetSuperdiagonalElement (int i) {
            return (U[i]);
        }

        internal double GetSubdiagonalElement (int i) {
            return (L[i]);
        }

        internal void SetDiagonalElement (int i, double value) {
            D[i] = value;
        }

        internal void SetSuperdiagonalElement (int i, double value) {
            U[i] = value;
        }

        internal void SetSubdiagonalElement (int i, double value) {
            L[i] = value;
        }

        private int BoundsCheck (int r, int c) {
            if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException("r");
            if ((c < 0) || (c >= dimension)) throw new ArgumentOutOfRangeException("c");
            return (c - r);
        }

        /// <summary>
        /// Gets or sets a matrix element.
        /// </summary>
        /// <param name="r">The (zero-based) row number of the element.</param>
        /// <param name="c">The (zero-based) column number of the element.</param>
        /// <returns>M<sub>r,c</sub></returns>
        /// <remarks>
        /// <para>Elements on the tridiagonal strip can be set and gotten normally. Other elements
        /// will always have the value zero, and any attempt to set them to a non-zero value will
        /// result in an <see cref="InvalidOperationException"/>.
        /// </para>
        /// </remarks>
        public double this[int r, int c] {
            get {
                switch (BoundsCheck(r, c)) {
                    case -1:
                        return (GetSubdiagonalElement(c));
                        break;
                    case 0:
                        return (GetDiagonalElement(r));
                        break;
                    case +1:
                        return (GetSuperdiagonalElement(r));
                        break;
                    default:
                        return (0.0);
                }
            }
            set {
                switch (BoundsCheck(r, c)) {
                    case -1:
                        SetSubdiagonalElement(c, value);
                        break;
                    case 0:
                        SetDiagonalElement(r, value);
                        break;
                    case +1:
                        SetSuperdiagonalElement(r, value);
                        break;
                    default:
                        if (value != 0.0) throw new InvalidOperationException();
                        break;
                }
            }
        }

        public TridiagonalMatrix Clone () {
            TridiagonalMatrix C = new TridiagonalMatrix(dimension);
            for (int i = 0; i < (dimension-1); i++) {
                C.SetDiagonalElement(i, GetDiagonalElement(i));
                C.SetSubdiagonalElement(i, GetSubdiagonalElement(i));
                C.SetSuperdiagonalElement(i, GetSuperdiagonalElement(i));
            }
            C.SetDiagonalElement(dimension - 1, GetDiagonalElement(dimension - 1));
            return (C);
        }

        IMatrix IMatrix.Clone () {
            return (Clone());
        }

        public TridiagonalMatrix Transpose () {
            TridiagonalMatrix T = new TridiagonalMatrix(dimension);
            for (int i = 0; i < (dimension - 1); i++) {
                T.SetDiagonalElement(i, GetDiagonalElement(i));
                T.SetSubdiagonalElement(i, GetSuperdiagonalElement(i));
                T.SetSuperdiagonalElement(i, GetSubdiagonalElement(i));
            }
            T.SetDiagonalElement(dimension - 1, GetDiagonalElement(dimension - 1));
            return (T);
        }

        IMatrix IMatrix.Transpose () {
            return (Transpose());
        }

        public double Determinant () {

            // the determinant is just the determinant of T, which can be computed by the recursion
            // det A_{n,n} = a_{n, n} det A_{n-1,n-1} - a_{n, n-1} a_{n-1, n} det A_{n-2,n-2}

            int n = Dimension;

            double a0 = GetDiagonalElement(0);

            if (n == 1) return (a0);

            double a1 = GetDiagonalElement(1) * GetDiagonalElement(0) - GetSubdiagonalElement(0) * GetSuperdiagonalElement(0);

            for (int i = 2; i < n; i++) {
                double a2 = GetDiagonalElement(i) * a1 - GetSubdiagonalElement(i - 1) * GetSuperdiagonalElement(i - 1) * a0;
                a0 = a1;
                a1 = a2;
            }

            return (a1);

        }

        public TridiagonalLUDecomposition LUDecompose () {

            double[] LC = (double[]) L.Clone();
            double[] DC = (double[]) D.Clone();
            double[] UC = (double[]) U.Clone();

            double[] V;
            int[] P;
            int pi;

            TridiagonalLUDecompose(LC, DC, UC, out V, out P, out pi);

            return (new TridiagonalLUDecomposition(LC, DC, UC, V, P, pi));

        }

        private static void TridiagonalLUDecompose (double[] L, double[] D, double[] U, out double[] V, out int[] P, out int pi) {

            int n = D.Length;

            Debug.Assert(L.Length == n - 1);
            Debug.Assert(D.Length == n);
            Debug.Assert(U.Length == n - 1);

            // keep track of the super-super-diagonal, which may be populated by row exchanges
            V = new double[n - 2];

            // keep track of row exchanges
            P = new int[n];
            for (int i = 0; i < n; i++) P[i] = i;
            pi = 1;

            // iterate across the columns
            for (int i = 0; i < (n - 1); i++) {

                Console.WriteLine("i = {0}", i);

                if (Math.Abs(L[i]) > Math.Abs(D[i])) {
                    // the subdiagonal element is larger, so permute it onto the diagonal for use as a pivot
                    double t;
                    t = D[i];
                    D[i] = L[i];
                    L[i] = t;
                    t = U[i];
                    U[i] = D[i + 1];
                    D[i + 1] = t;
                    if (i < (n - 2)) {
                        V[i] = U[i + 1];
                        U[i + 1] = 0.0;
                    }

                    P[i] = i + 1;
                    //int tp;
                    //tp = P[i];
                    //P[i] = P[i + 1];
                    //P[i + 1] = tp;

                }

                Console.WriteLine("exchanged");

                // Here is what we are working on...
                // (     D_{i-1}  U_{i-1}  V_{i-1}                       )
                // (     0        D_{i}    U_{i}    V_{i}                )  <- pivot row
                // (              L_{i}    D_{i+1}  U_{i+1}              )
                // (                       L_{i+1}  D_{i+2}  U_{i+2}     )
                // The last iteration left L_{i-1}=0. In this iteration, we will make L_{i}=0.
                // We use L_{j} for j < i to store the factor of the pivot used to reduce the subdiagonal elements

                // zero the subdiagonal element, and keep track of the effect on other elements of that row
                double f = L[i] / D[i];
                L[i] = 0.0;
                D[i + 1] = D[i + 1] - f * U[i];
                if (i < (n - 2)) U[i + 1] = U[i + 1] - f * V[i];

                // store the reduction factor
                L[i] = f;

                for (int j = 0; j < n; j++) {
                    Console.Write("  {0}", P[j]);
                }
                Console.WriteLine();

                for (int j = 0; j < n - 2; j++) {
                    Console.Write("  {0}", V[j]);
                }
                Console.WriteLine();

                for (int j = 0; j < n - 1; j++) {
                    Console.Write("  {0}", U[j]);
                }
                Console.WriteLine();

                for (int j = 0; j < n; j++) {
                    Console.Write("  {0}", D[j]);
                }
                Console.WriteLine();

                for (int j = 0; j < n - 1; j++) {
                    Console.Write("  {0}", L[j]);
                }
                Console.WriteLine();

            }

        }

        // tridiagonal matrix arithmetic

        // ( X X       ) ( Y Y       )
        // ( X X X     ) ( Y Y Y     )
        // (   X X X   ) (   Y Y Y   )
        // (     X X X ) (     Y Y Y )
        // (       X X ) (       Y Y )

        public static TridiagonalMatrix operator + (TridiagonalMatrix T1, TridiagonalMatrix T2) {

            int n = T1.Dimension;
            if (T2.Dimension != n) throw new DimensionMismatchException();

            TridiagonalMatrix T = new TridiagonalMatrix(n);
            for (int i = 0; i < n; i++) {
                T[i,i] = T1[i,i-1] * T2[i-1,i] + T1[i,i] * T2[i,i] + T1[i,i+1] * T2[i+1,i];
            }

            return (T);

        }


    }

    /// <summary>
    /// Represents the LU decomposition of a tridiagonal matrix.
    /// </summary>
    public class TridiagonalLUDecomposition {

        private int n;

        private int[] P;
        int parity;

        private double[] L;
        private double[] D;
        private double[] U;
        private double[] V;

        internal TridiagonalLUDecomposition (double[] L, double[] D, double[] U, double[] V, int[] P, int parity) {

            this.n = D.Length;

            this.L = L;
            this.D = D;
            this.U = U;
            this.V = V;
            this.P = P;
            this.parity = parity;

            Debug.Assert(L.Length == n - 1);
            Debug.Assert(U.Length == n - 1);

            Debug.Assert(P.Length == n);

        }

        public int Dimension {
            get {
                return (P.Length);
            }
        }

        public double Determinant () {
            double det = parity;
            for (int i = 0; i < n; i++) {
                det = det * D[i];
            }
            return (det);
        }

        public ColumnVector Solve (IList<double> rhs) {

            if (rhs.Count != n) throw new DimensionMismatchException();

            Console.WriteLine("fill in");

            ColumnVector x = new ColumnVector(n);
            for (int i = 0; i < n; i++) {
                x[i] = rhs[i];
            }

            PrintVector(x);

            Console.WriteLine("forward sub");

            // forward-substitute to solve Ly=z
            for (int i = 0; i < (n-1); i++) {

                double t = x[i+1-P[i]+i] - L[i] * x[P[i]];
                x[i] = x[P[i]];
                x[i+1] = t;


                //x[i] = x[i] - L[i - 1] * x[i - 1];
            }

            PrintVector(x);


            Console.WriteLine("back sub");

            // back-substitute to solve Ux=y
            x[n - 1] = x[n - 1] / D[n - 1];
            if (n > 1) {
                x[n - 2] = (x[n - 2] - U[n - 2] * x[n - 1]) / D[n - 2];
                for (int i = n - 3; i >= 0; i--) {
                    x[i] = (x[i] - U[i] * x[i + 1] - V[i] * x[i + 2]) / D[i];
                }
            }

            PrintVector(x);


            return (x);

        }

        private void PrintVector (IList<double> v) {
            for (int i = 0; i < v.Count; i++) {
                Console.Write("  {0}", v[i]);
            }
            Console.WriteLine();
        }

        public SquareMatrix Inverse () {

            SquareMatrix TI = new SquareMatrix(n);
            for (int c = 0; c < n; c++) {

                ColumnVector e = new ColumnVector(n);
                e[c] = 1.0;

                ColumnVector v = Solve(e);

                for (int r = 0; r < n; r++) {
                    TI[r, c] = v[r];
                }

            }

            return (TI);

        }

        public SquareMatrix PMatrix () {
            SquareMatrix PM = new SquareMatrix(n);
            for (int r = 0; r < n; r++) {
                PM[r, P[r]] = 1.0;
            }
            return (PM);
        }

        public SquareMatrix LMatrix () {
            SquareMatrix LM = new SquareMatrix(n);
            for (int i = 0; i < (n - 1); i++) {
                LM[i, i] = 1.0;
                LM[i + 1, i] = L[i];
            }
            LM[n - 1, n - 1] = 1.0;
            return (LM);
        }

        public SquareMatrix UMatrix () {
            SquareMatrix UM = new SquareMatrix(n);
            for (int i = 0; i < n; i++) {
                UM[i, i] = D[i];
                if (i < (n-1)) UM[i, i+1] = U[i];
                if (i < (n-2)) UM[i, i+2] = V[i];
            }
            return (UM);
        }

    }

#endif

}
