using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

namespace Meta.Numerics.Matrices {


    /// <summary>
    /// Represents a tri-diagonal matrix.
    /// </summary>
    public sealed class TridiagonalMatrix : AnySquareMatrix {

        /// <summary>
        /// Initializes a new tri-diagonal matrix of the given dimension.
        /// </summary>
        /// <param name="dimension">The dimension of the matrix, which must be positive.</param>
        public TridiagonalMatrix (int dimension) {
            if (dimension < 1) throw new ArgumentOutOfRangeException(nameof(dimension));
            this.dimension = dimension;
            superDiag = new double[dimension - 1];
            diag = new double[dimension];
            subDiag = new double[dimension - 1];
        }

        internal TridiagonalMatrix (int dimension, double[] superDiag, double[] diag, double[] subDiag, bool isReadOnly) : base(isReadOnly) {
            Debug.Assert(dimension > 0);
            Debug.Assert(superDiag != null);
            Debug.Assert(diag != null);
            Debug.Assert(subDiag != null);
            Debug.Assert(superDiag.Length == dimension - 1);
            Debug.Assert(diag.Length == dimension);
            Debug.Assert(subDiag.Length == dimension - 1);
            this.dimension = dimension;
            this.superDiag = superDiag;
            this.diag = diag;
            this.subDiag = subDiag;
        }

        private readonly int dimension;
        private readonly double[] superDiag;
        private readonly double[] diag;
        private readonly double[] subDiag;

        /// <summary>
        /// Gets the dimension of the matrix.
        /// </summary>
        public override int Dimension {
            get {
                return (dimension);
            }
        }

        internal double GetDiagonalElement (int i) {
            return (diag[i]);
        }

        internal double GetSuperdiagonalElement (int i) {
            return (superDiag[i]);
        }

        internal double GetSubdiagonalElement (int i) {
            return (subDiag[i]);
        }

        internal void SetDiagonalElement (int i, double value) {
            diag[i] = value;
        }

        internal void SetSuperdiagonalElement (int i, double value) {
            superDiag[i] = value;
        }

        internal void SetSubdiagonalElement (int i, double value) {
            subDiag[i] = value;
        }

        private int BoundsCheck (int r, int c) {
            if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException(nameof(r));
            if ((c < 0) || (c >= dimension)) throw new ArgumentOutOfRangeException(nameof(c));
            return (c - r);
        }

        /// <summary>
        /// Gets or sets a matrix element.
        /// </summary>
        /// <param name="r">The (zero-based) row number of the element.</param>
        /// <param name="c">The (zero-based) column number of the element.</param>
        /// <returns>M<sub>r,c</sub></returns>
        /// <remarks>
        /// <para>Get and set operates normally for elements on the tri-diagonal strip. Other elements
        /// will always have the value zero, and any attempt to set them to a non-zero value will
        /// result in an <see cref="InvalidOperationException"/>.
        /// </para>
        /// </remarks>
        public override double this[int r, int c] {
            get {
                switch (BoundsCheck(r, c)) {
                    case -1:
                        return (GetSubdiagonalElement(c));
                    case 0:
                        return (GetDiagonalElement(r));
                    case +1:
                        return (GetSuperdiagonalElement(r));
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

        /// <summary>
        /// Copies the matrix.
        /// </summary>
        /// <returns>An independent copy of the matrix.</returns>
        public TridiagonalMatrix Copy () {
            return (new TridiagonalMatrix(dimension, CopyArray(superDiag), CopyArray(diag), CopyArray(subDiag), false));
        }

        private double[] CopyArray (double[] array) {
            double[] copy = new double[array.Length];
            Array.Copy(array, copy, array.Length);
            return (copy);
        }

        /// <summary>
        /// Creates a transpose of the matrix.
        /// </summary>
        /// <value>The matrix transpose M<sup>T</sup>.</value>
        /// <remarks>
        /// <para>The returned transpose matrix is not independent of the original matrix.
        /// Instead, it is a read-only view into the same storage as the original matrix with row an column indices reversed.
        /// This has the advantage that is can be produced with almost no time and memory cost.
        /// It has the disadvantage that any subsequent changes to the original matrix will be reflected in the returned transpose,
        /// which might not be what you expected. If you want an independent, write-able transpose matrix, call <see cref="Copy"/>
        /// on the returned matrix.
        /// </para>
        /// </remarks>
        public TridiagonalMatrix Transpose {
            get {
                // switch super-diagonal and sub-diagonal
                return (new TridiagonalMatrix(dimension, subDiag, diag, superDiag, true));
            }
        }


        /// <summary>
        /// Computes the determinant of the matrix.
        /// </summary>
        /// <returns>The determinant det M.</returns>
        /// <remarks>
        /// <para>Computing the determinant of a tri-diagonal matrix is an O(N) operation.</para>
        /// </remarks>
        public double Determinant () {

            // the determinant is just the determinant of T, which can be computed by the recursion
            // det A_{n,n} = a_{n, n} det A_{n-1,n-1} - a_{n, n-1} a_{n-1, n} det A_{n-2,n-2}

            double a0 = GetDiagonalElement(0);
            if (dimension == 1) return (a0);
            double a1 = GetDiagonalElement(1) * GetDiagonalElement(0) - GetSubdiagonalElement(0) * GetSuperdiagonalElement(0);

            for (int i = 2; i < dimension; i++) {
                double a2 = GetDiagonalElement(i) * a1 - GetSubdiagonalElement(i - 1) * GetSuperdiagonalElement(i - 1) * a0;
                a0 = a1;
                a1 = a2;
            }

            return (a1);

        }

        /// <summary>
        /// Computes the LU decomposition of the matrix.
        /// </summary>
        /// <returns>The LU decomposition of the matrix.</returns>
        /// <remarks>
        /// <para>Computing the LU decomposition of a tri-diagonal matrix is an O(N) operation.</para>
        /// </remarks>
        public TridiagonalLUDecomposition LUDecomposition () {

            double[] LC = CopyArray(subDiag);
            double[] DC = CopyArray(diag);
            double[] UC = CopyArray(superDiag);

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

                if (Math.Abs(L[i]) > Math.Abs(D[i])) {
                    // the sub-diagonal element is larger, so permute it onto the diagonal for use as a pivot
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
                    pi = -pi;

                }

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

            }

        }

        // tri- diagonal matrix arithmetic

        // ( X X       ) ( Y Y       )
        // ( X X X     ) ( Y Y Y     )
        // (   X X X   ) (   Y Y Y   )
        // (     X X X ) (     Y Y Y )
        // (       X X ) (       Y Y )

        /// <summary>
        /// Adds two tri-diagonal matrices.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>The matrix sum A + B.</returns>
        public static TridiagonalMatrix operator + (TridiagonalMatrix A, TridiagonalMatrix B) {

            if (A == null) throw new ArgumentNullException(nameof(A));
            if (B == null) throw new ArgumentNullException(nameof(B));

            int n = A.Dimension;
            if (B.Dimension != n) throw new DimensionMismatchException();

            TridiagonalMatrix T = new TridiagonalMatrix(n);

            // add the elements
            // this is a glorified for-next loop structured so as to minimize the
            // integer arithmetic required to keep track of two numbers one unit apart
            // and avoid if-tests inside the loop
            T.SetDiagonalElement(0, A.GetDiagonalElement(0) + B.GetDiagonalElement(0));
            int i0 = 0;
            int i1 = 1;
            while (i1 < n) {
                T.SetDiagonalElement(i1, A.GetDiagonalElement(i1) + B.GetDiagonalElement(i1));
                T.SetSubdiagonalElement(i0, A.GetSubdiagonalElement(i0) + B.GetSubdiagonalElement(i0));
                T.SetSuperdiagonalElement(i0, A.GetSuperdiagonalElement(i0) + B.GetSuperdiagonalElement(i0));
                i0 = i1;
                i1++;
            }

            return (T);

        }

        /// <summary>
        /// Subtracts two tri-diagonal matrices.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>The difference matrix A - B.</returns>
        public static TridiagonalMatrix operator - (TridiagonalMatrix A, TridiagonalMatrix B) {

            if (A == null) throw new ArgumentNullException(nameof(A));
            if (B == null) throw new ArgumentNullException(nameof(B));

            int n = A.Dimension;
            if (B.Dimension != n) throw new DimensionMismatchException();

            TridiagonalMatrix T = new TridiagonalMatrix(n);

            // add the elements
            // this is a glorified for-next loop structured so as to minimize the
            // integer arithmetic required to keep track of two numbers one unit apart
            // and avoid if-tests inside the loop
            T.SetDiagonalElement(0, A.GetDiagonalElement(0) - B.GetDiagonalElement(0));
            int i0 = 0;
            int i1 = 1;
            while (i1 < n) {
                T.SetDiagonalElement(i1, A.GetDiagonalElement(i1) - B.GetDiagonalElement(i1));
                T.SetSubdiagonalElement(i0, A.GetSubdiagonalElement(i0) - B.GetSubdiagonalElement(i0));
                T.SetSuperdiagonalElement(i0, A.GetSuperdiagonalElement(i0) - B.GetSuperdiagonalElement(i0));
                i0 = i1;
                i1++;
            }

            return (T);

        }

        /// <summary>
        /// Multiplies a tri-diagonal matrix by a real constant.
        /// </summary>
        /// <param name="f">The constant.</param>
        /// <param name="T">The matrix.</param>
        /// <returns>The product matrix.</returns>
        public static TridiagonalMatrix operator * (double f, TridiagonalMatrix T) {

            if (T == null) throw new ArgumentNullException(nameof(T));

            int n = T.Dimension;
            TridiagonalMatrix fT = new TridiagonalMatrix(n);

            fT.SetDiagonalElement(0, f * T.GetDiagonalElement(0));

            int i0 = 0;
            int i1 = 1;
            while (i1 < n) {
                fT.SetDiagonalElement(i1, f * T.GetDiagonalElement(i1));
                fT.SetSubdiagonalElement(i0, f * T.GetSubdiagonalElement(i0));
                fT.SetSuperdiagonalElement(i0, f * T.GetSuperdiagonalElement(i0));
                i0 = i1;
                i1++;
            }

            return (fT);

        }

        // equality testing
        /*
        private static bool Equals (TridiagonalMatrix T1, TridiagonalMatrix T2) {
            if (((object) T1) == null) {
                if (((object) T2) == null) {
                    return (true);
                } else {
                    return (false);
                }
            } else {
                if (((object) T2) == null) {
                    return (false);
                } else {
                    int n = T1.Dimension;
                    if (T2.Dimension != n) throw new DimensionMismatchException();

                    if (T1.GetDiagonalElement(0) != T2.GetDiagonalElement(0)) return (false);
                    // this is a glorified for-next loop structured so as to minimize the
                    // integer arithmetic required to keep track of two numbers one unit apart
                    int i0 = 0;
                    int i1 = 1;
                    while (i1 < n) {
                        if (T1.GetDiagonalElement(i1) != T2.GetDiagonalElement(i1)) return (false);
                        if (T1.GetSubdiagonalElement(i0) != T2.GetSubdiagonalElement(i0)) return (false);
                        if (T1.GetSuperdiagonalElement(i0) != T2.GetSuperdiagonalElement(i0)) return (false);
                        i0 = i1;
                        i1++;
                    }
                    return (true);
                }
            }
        }

        /// <summary>
        /// Determines whether two tridiagonal matrices are equal.
        /// </summary>
        /// <param name="T1">The first matrix.</param>
        /// <param name="T2">The second matrix.</param>
        /// <returns>True if <paramref name="T1"/> and <paramref name="T2"/> are equal, otherwise false.</returns>
        public static bool operator == (TridiagonalMatrix T1, TridiagonalMatrix T2) {
            return (TridiagonalMatrix.Equals(T1, T2));
        }

        /// <summary>
        /// Determines whether two tridiagonal matrices are not equal.
        /// </summary>
        /// <param name="T1">The first matrix.</param>
        /// <param name="T2">The second matrix.</param>
        /// <returns>False if <paramref name="T1"/> and <paramref name="T2"/> are equal, otherwise true.</returns>
        public static bool operator != (TridiagonalMatrix T1, TridiagonalMatrix T2) {
            return (!TridiagonalMatrix.Equals(T1, T2));
        }

        /// <summary>
        /// Determines whether the given object is an equal matrix.
        /// </summary>
        /// <param name="obj">The object to compare.</param>
        /// <returns>True if <paramref name="obj"/> is an equal matrix, otherwise false.</returns>
        public override bool Equals (object obj) {
            return (Matrix.Equals(this, obj as IMatrix));
        }

        /// <inheritdoc />
        public override int GetHashCode () {
            return base.GetHashCode();
        }
        

#if SHO
        [Obsolete]
        public string __repr__ () {
            StringWriter writer = new StringWriter();
            Matrix.WriteMatrix(this, writer);
            return (writer.ToString());
        }
#endif
        */

    }

    /// <summary>
    /// Represents the LU decomposition of a tri-diagonal matrix.
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

        /// <summary>
        /// Gets the dimension of the original matrix.
        /// </summary>
        public int Dimension {
            get {
                return (P.Length);
            }
        }

        /// <summary>
        /// Computes the determinant of the original matrix.
        /// </summary>
        /// <returns>The determinant det M.</returns>
        public double Determinant () {
            double det = parity;
            for (int i = 0; i < n; i++) {
                det = det * D[i];
            }
            return (det);
        }

        /// <summary>
        /// Solves a tri-diagonal system of linear equations.
        /// </summary>
        /// <param name="rhs">The right-hand side vector b.</param>
        /// <returns>A vector x which satisfies Ax = b.</returns>
        public ColumnVector Solve (IList<double> rhs) {

            if (rhs == null) throw new ArgumentNullException(nameof(rhs));
            if (rhs.Count != n) throw new DimensionMismatchException();

            ColumnVector x = new ColumnVector(n);
            for (int i = 0; i < n; i++) {
                x[i] = rhs[i];
            }

            // forward-substitute to solve Ly=z
            for (int i = 0; i < (n-1); i++) {

                double t = x[i+1-P[i]+i] - L[i] * x[P[i]];
                x[i] = x[P[i]];
                x[i+1] = t;
            }

            // back-substitute to solve Ux=y
            x[n - 1] = x[n - 1] / D[n - 1];
            if (n > 1) {
                x[n - 2] = (x[n - 2] - U[n - 2] * x[n - 1]) / D[n - 2];
                for (int i = n - 3; i >= 0; i--) {
                    x[i] = (x[i] - U[i] * x[i + 1] - V[i] * x[i + 2]) / D[i];
                }
            }

            return (x);

        }

        /// <summary>
        /// Computes the inverse of the original matrix.
        /// </summary>
        /// <returns>The matrix M<sup>-1</sup>.</returns>
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

        /*
        ISquareMatrix ISquareDecomposition.Inverse () {
            return (Inverse());
        }
        */

        /*
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


        /// <summary>
        /// The matrix U in the decomposition PA = LU.
        /// </summary>
        /// <returns></returns>
        public SquareMatrix UMatrix () {
            SquareMatrix UM = new SquareMatrix(n);
            for (int i = 0; i < n; i++) {
                UM[i, i] = D[i];
                if (i < (n-1)) UM[i, i+1] = U[i];
                if (i < (n-2)) UM[i, i+2] = V[i];
            }
            return (UM);
        }
        */

    }

}
