using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents a vector of numbers.
    /// </summary>
    /// <typeparam name="T">The type of numbers the vector contins.</typeparam>
    /// <remarks><para>For matrix operations, you will usually want to define a ColumnVector or RowVector instead.</para></remarks>
    /// <seealso cref="ColumnVector"/>
    /// <seealso cref="RowVector" />
    public class Vector<T> : IList<T>, ICollection<T>, IEnumerable<T>, IEnumerable where T : struct {

        private T[] v;

        /// <summary>
        /// Instantiates a new vector of the given dimension.
        /// </summary>
        /// <param name="dimension">The dimension of the vector.</param>
        public Vector (int dimension) {
            if (dimension < 1) throw new ArgumentOutOfRangeException("dimension");
            v = new T[dimension];
        }


        /// <summary>
        /// Instantiates a new vector containing the given elements.
        /// </summary>
        /// <param name="list">The component elements.</param>
        public Vector (IList<T> list) {
            if (list == null) throw new ArgumentNullException("list");
            v = new T[list.Count];
            list.CopyTo(v, 0);
        }

        /// <summary>
        /// Gets or sets the specified vector component.
        /// </summary>
        /// <param name="index">The (zero-based) component index.</param>
        /// <returns>The value of the specified vector component.</returns>
        public T this[int index] {
            get {
                return (v[index]);
            }
            set {
                v[index] = value;
            }

        }

        /// <summary>
        /// Gets the dimension of the vector.
        /// </summary>
        public int Dimension {
            get {
                return (v.Length);
            }
        }

        /// <summary>
        /// Converts the vector to an array.
        /// </summary>
        /// <returns>An array of the vector components.</returns>
        public T[] ToArray () {
            T[] u = new T[v.Length];
            v.CopyTo(u, 0);
            return (u);
        }

        // operators

        // equality

        internal static bool Equals (Vector<T> v1, Vector<T> v2) {

            bool v1null = (((object) v1) == null);
            bool v2null = (((object) v2) == null);

            // this is elegant, but too elegant to be clear
            // if (v1null && v2null) return (true);
            // if (v1null ^ v2null) return (false);
            // neither v1 nor v2 is null

            if (v1null) {
                if (v2null) {
                    return (true);
                } else {
                    return (false);
                }
            } else {
                if (v2null) {
                    return (false);
                } else {
                    if (v1.Dimension != v2.Dimension) throw new DimensionMismatchException();
                    for (int i = 0; i < v1.Dimension; i++) {
                        //if (v1[i] != v2[i]) return (false);
                        if (!Object.Equals(v1[i],v2[i])) return(false);
                    }
                    return (true);
                }
            }

        }

        internal static int GetHashCode (Vector<T> v) {
            int hash = 1;
            for (int i = 0; i < v.Dimension; i++) {
                hash = hash ^ (v[i].GetHashCode() + i);
            }
            return (hash);
        }

        // IList methods

        int IList<T>.IndexOf (T item) {
            return (((IList<T>) v).IndexOf(item));
        }

        void IList<T>.Insert (int index, T item) {
            ((IList<T>) v).Insert(index, item);
        }

        void IList<T>.RemoveAt (int index) {
            ((IList<T>) v).RemoveAt(index);
        }

        // ICollection methods

        void ICollection<T>.Add (T item) {
            ((ICollection<T>) v).Add(item);
        }

        void ICollection<T>.Clear () {
            ((ICollection<T>) v).Clear();
        }

        bool ICollection<T>.Contains (T item) {
            return (((ICollection<T>) v).Contains(item));
        }

        void ICollection<T>.CopyTo (T[] array, int arrayIndex) {
            ((ICollection<T>) v).CopyTo(array, arrayIndex);
        }

        bool ICollection<T>.Remove (T item) {
            return (((ICollection<T>) v).Remove(item));
        }

        int ICollection<T>.Count {
            get {
                return(((ICollection<T>) v).Count);
            }
        }

        bool ICollection<T>.IsReadOnly {
            get {
                return (((ICollection<T>) v).IsReadOnly);
            }
        }

        // IEnumerable methods

        IEnumerator<T> IEnumerable<T>.GetEnumerator () {
            return (((IEnumerable<T>) v).GetEnumerator());
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (((IEnumerable) v).GetEnumerator());
        }
    }

    /// <summary>
    /// Represents a column vector of real numbers.
    /// </summary>
    /// <remarks><para>An N-dimensional column vector is an NX1-dimensional matrix.</para></remarks>
    public sealed class ColumnVector : Vector<double>, IMatrix {

        /// <summary>
        /// Initializes a column vector of a given dimension.
        /// </summary>
        /// <param name="dimension">The dimension of the vector.</param>
        public ColumnVector (int dimension) : base(dimension) {
        }

        /// <summary>
        /// Initializes a column vector with the components given by a list.
        /// </summary>
        /// <param name="list">A list of the initial components.</param>
        /// <remarks>
        /// <para>Many ordered collections of reals, including double[] and List&lt;double&gt;, implement
        /// IList&lt;double&gt; and thus can be used with this constructor.</para>
        /// </remarks>
        public ColumnVector (IList<double> list) : base(list) {
        }

        /// <summary>
        /// Returns the transpose of the column vector.
        /// </summary>
        /// <returns>An independent row vector with the same components as the column vector.</returns>
        public RowVector Transpose () {
            RowVector u = new RowVector(Dimension);
            for (int i = 0; i < Dimension; i++) {
                u[i] = this[i];
            }
            return (u);
        }

        /// <summary>
        /// Returns a clone of the column vector.
        /// </summary>
        /// <returns>An independent column vector with the same components as the original.</returns>
        public ColumnVector Clone () {
            ColumnVector u = new ColumnVector(Dimension);
            for (int i = 0; i < Dimension; i++) {
                u[i] = this[i];
            }
            return (u);
        }

        // operators

        // equality

        /// <summary>
        /// Determines whether two column vectors are equal.
        /// </summary>
        /// <param name="v1">The first column vector.</param>
        /// <param name="v2">The second column vector.</param>
        /// <returns>True if <paramref name="v1"/> and <paramref name="v2"/> are equal, otherwise false.</returns>
        public static bool operator == (ColumnVector v1, ColumnVector v2) {
            return (Vector<double>.Equals(v1, v2));
        }

        /// <summary>
        /// Determines whether two column vectors are not equal.
        /// </summary>
        /// <param name="v1">The first column vector.</param>
        /// <param name="v2">The second column vector.</param>
        /// <returns>False if <paramref name="v1"/> and <paramref name="v2"/> are equal, otherwise true.</returns>
        public static bool operator != (ColumnVector v1, ColumnVector v2) {
            return (!Vector<double>.Equals(v1, v2));
        }

        /// <summary>
        /// Determine whether the supplied object is an equal column vector.
        /// </summary>
        /// <param name="obj">An object.</param>
        /// <returns>True if <paramref name="obj"/> is an equal column vector, otherwise false.</returns>
        public override bool Equals (object obj) {
            return (Vector<double>.Equals(this, obj as ColumnVector));
        }

        /// <inheritdoc />
        public override int GetHashCode () {
            return (Vector<double>.GetHashCode(this));
        }

        // arithmatic

        /// <summary>
        /// Computes the sum of two column vectors.
        /// </summary>
        /// <param name="v1">The first column vector.</param>
        /// <param name="v2">The second column vector.</param>
        /// <returns>The sum of <paramref name="v1"/> and <paramref name="v2"/>.</returns>
        public static ColumnVector operator + (ColumnVector v1, ColumnVector v2) {
            if (v1.Dimension != v2.Dimension) throw new DimensionMismatchException();
            int d = v1.Dimension;
            ColumnVector u = new ColumnVector(d);
            for (int i = 0; i < d; i++) {
                u[i] = v1[i] + v2[i];
            }
            return (u);
        }

        /// <summary>
        /// Computes the difference of two column vectors.
        /// </summary>
        /// <param name="v1">The first column vector.</param>
        /// <param name="v2">The second column vector.</param>
        /// <returns>The difference of <paramref name="v1"/> and <paramref name="v2"/>.</returns>
        public static ColumnVector operator - (ColumnVector v1, ColumnVector v2) {
            if (v1 == null) throw new ArgumentNullException("v1");
            if (v2 == null) throw new ArgumentNullException("v2");
            if (v1.Dimension != v2.Dimension) throw new DimensionMismatchException();
            int d = v1.Dimension;
            ColumnVector u = new ColumnVector(d);
            for (int i = 0; i < d; i++) {
                u[i] = v1[i] - v2[i];
            }
            return (u);
        }

        /// <summary>
        /// Computes the product of a real number and a column vector.
        /// </summary>
        /// <param name="a">The real number.</param>
        /// <param name="v">The column vector.</param>
        /// <returns>The product of <paramref name="a"/> and <paramref name="v"/>.</returns>
        public static ColumnVector operator * (double a, ColumnVector v) {
            int d = v.Dimension;
            ColumnVector u = new ColumnVector(d);
            for (int i = 0; i < d; i++) {
                u[i] = a * v[i];
            }
            return (u);
        }

        /*
        public static ColumnVector operator * (ColumnVector v, double a) {
            return (a * v);
        }
        */

        /*
        public static ColumnVector operator / (ColumnVector v, double a) {
            int d = v.Dimension;
            ColumnVector u = new ColumnVector(d);
            for (int i = 0; i < d; i++) {
                u[i] = v[i] / a;
            }
            return (u);
        }
        */

        /// <summary>
        /// Computes the additive inverse of a column vector.
        /// </summary>
        /// <param name="v">The column vector.</param>
        /// <returns>The additive inverse of the column vector.</returns>
        public static ColumnVector operator - (ColumnVector v) {
            int d = v.Dimension;
            ColumnVector u = new ColumnVector(d);
            for (int i = 0; i < d; i++) {
                u[i] = -v[i];
            }
            return (u);
        }

        // dot product

        /// <summary>
        /// Computes the inner product of a row and a column vector.
        /// </summary>
        /// <param name="v1">The row vector.</param>
        /// <param name="v2">The column vector.</param>
        /// <returns>The (dot) inner product v1 * v2.</returns>
        public static double operator * (RowVector v1, ColumnVector v2) {
            if (v1.Dimension != v2.Dimension) throw new DimensionMismatchException();
            int d = v1.Dimension;
            double u = 0.0;
            for (int i = 0; i < d; i++) {
                u += v1[i] * v2[i];
            }
            return (u);
        }

        /// <summary>
        /// Computes the product of a matrix and a column vector.
        /// </summary>
        /// <param name="M">The matrix.</param>
        /// <param name="v">The column vector.</param>
        /// <returns>The product Mv.</returns>
        public static ColumnVector operator * (IMatrix M, ColumnVector v) {
            if (M.ColumnCount != v.Dimension) throw new DimensionMismatchException();
            ColumnVector u = new ColumnVector(M.RowCount);
            for (int i = 0; i < u.Dimension; i++) {
                u[i] = 0.0;
                for (int j = 0; j < v.Dimension; j++) {
                    u[i] += M[i, j] * v[j];

                }
            }
            return(u);
        }

        // additional IMatrix functions

        int IMatrix.ColumnCount {
            get {
                return (1);
            }
        }

        int IMatrix.RowCount {
            get {
                return (Dimension);
            }
        }

        IMatrix IMatrix.Transpose () {
            return (this.Transpose());
        }

        IMatrix IMatrix.Clone () {
            return (this.Clone());
        }

        double IMatrix.this[int r, int c] {
            get {
                if ((r < 0) || (r > Dimension)) throw new ArgumentOutOfRangeException("r");
                if (c != 0) throw new ArgumentOutOfRangeException("c");
                return (this[r]);
            }
        }

#if SHO
        /// <summary>
        /// Produces the representation of the vector for the Python interactive console.
        /// </summary>
        /// <returns>A string representation of the vector.</returns>
        public string __repr__ () {
            StringWriter writer = new StringWriter();
            Matrix.WriteMatrix(this, writer);
            return (writer.ToString());
        }
#endif

    }

    /// <summary>
    /// Represents a row vector of real numbers.
    /// </summary>
    /// <remarks><para>An N-dimensional row vector is a 1XN-dimensional matrix.</para></remarks>
    public sealed class RowVector : Vector<double>, IMatrix {

        /// <summary>
        /// Initializes a row vector of a given dimension.
        /// </summary>
        /// <param name="dimension">The dimension of the vector.</param>
        public RowVector (int dimension) : base(dimension) {
        }

        /// <summary>
        /// Initializes a row vector with the components given by a list.
        /// </summary>
        /// <param name="list">A list of the initial components.</param>
        /// <remarks>
        /// <para>Many ordered collections of reals, including double[] and List&lt;double&gt;, implement
        /// IList&lt;double&gt; and thus can be used with this constructor.</para>
        /// </remarks>
        public RowVector (IList<double> list) : base(list) {
        }

        /// <summary>
        /// Returns the transpose of the column vector.
        /// </summary>
        /// <returns>An independent column vector with the same components as the row vector.</returns>
        public ColumnVector Transpose () {
            ColumnVector u = new ColumnVector(Dimension);
            for (int i = 0; i < Dimension; i++) {
                u[i] = this[i];
            }
            return (u);
        }

        /// <summary>
        /// Returns a clone of the row vector.
        /// </summary>
        /// <returns>An independent row vector with the same components as the original.</returns>
        public RowVector Clone () {
            RowVector u = new RowVector(Dimension);
            for (int i = 0; i < Dimension; i++) {
                u[i] = this[i];
            }
            return (u);
        }


        // additional IMatrix functions

        int IMatrix.ColumnCount {
            get {
                return (Dimension);
            }
        }

        int IMatrix.RowCount {
            get {
                return (1);
            }
        }

        IMatrix IMatrix.Transpose () {
            return (this.Transpose());
        }

        IMatrix IMatrix.Clone () {
            return (this.Clone());
        }

        double IMatrix.this[int r, int c] {
            get {
                if (r != 0) throw new ArgumentOutOfRangeException("c");
                if ((c < 0) || (c > Dimension)) throw new ArgumentOutOfRangeException("r");
                return (this[c]);
            }
        }

        // operators

        // equality

        /// <summary>
        /// Determines whether two row vectors are equal.
        /// </summary>
        /// <param name="v1">The first row vector.</param>
        /// <param name="v2">The second row vector.</param>
        /// <returns>True if <paramref name="v1"/> and <paramref name="v2"/> are equal, otherwise false.</returns>
        public static bool operator == (RowVector v1, RowVector v2) {
            return (Vector<double>.Equals(v1, v2));
        }

        /// <summary>
        /// Determines whether two row vectors are not equal.
        /// </summary>
        /// <param name="v1">The first row vector.</param>
        /// <param name="v2">The second row vector.</param>
        /// <returns>False if <paramref name="v1"/> and <paramref name="v2"/> are equal, otherwise true.</returns>
        public static bool operator != (RowVector v1, RowVector v2) {
            return (!Vector<double>.Equals(v1, v2));
        }

        /// <summary>
        /// Determine whether the supplied object is an equal row vector.
        /// </summary>
        /// <param name="obj">An object.</param>
        /// <returns>True if <paramref name="obj"/> is an equal row vector, otherwise false.</returns>
        public override bool Equals (object obj) {
            return (Vector<double>.Equals(this, obj as RowVector));
        }

        /// <inheritdoc />
        public override int GetHashCode () {
            return (Vector<double>.GetHashCode(this));
        }

        // arithmetic

        /// <summary>
        /// Computes the sum of two row vectors.
        /// </summary>
        /// <param name="v1">The first row vector.</param>
        /// <param name="v2">The second row vector.</param>
        /// <returns>The sum of <paramref name="v1"/> and <paramref name="v2"/>.</returns>
        public static RowVector operator + (RowVector v1, RowVector v2) {
            if (v1.Dimension != v2.Dimension) throw new DimensionMismatchException();
            int d = v1.Dimension;
            RowVector u = new RowVector(d);
            for (int i = 0; i < d; i++) {
                u[i] = v1[i] + v2[i];
            }
            return (u);
        }

        /// <summary>
        /// Computes the difference of two row vectors.
        /// </summary>
        /// <param name="v1">The first row vector.</param>
        /// <param name="v2">The second row vector.</param>
        /// <returns>The difference of <paramref name="v1"/> and <paramref name="v2"/>.</returns>
        public static RowVector operator - (RowVector v1, RowVector v2) {
            if (v1.Dimension != v2.Dimension) throw new DimensionMismatchException();
            int d = v1.Dimension;
            RowVector u = new RowVector(d);
            for (int i = 0; i < d; i++) {
                u[i] = v1[i] - v2[i];
            }
            return (u);
        }

        /// <summary>
        /// Computes the product of a real number and a row vector.
        /// </summary>
        /// <param name="a">The real number.</param>
        /// <param name="v">The row vector.</param>
        /// <returns>The product of <paramref name="a"/> and <paramref name="v"/>.</returns>
        public static RowVector operator * (double a, RowVector v) {
            int d = v.Dimension;
            RowVector u = new RowVector(d);
            for (int i = 0; i < d; i++) {
                u[i] = a * v[i];
            }
            return (u);
        }

        // outer product

        /// <summary>
        /// Computes the outer product of a column vector and a row vector.
        /// </summary>
        /// <param name="v">The column vector.</param>
        /// <param name="u">The row vector.</param>
        /// <returns>The outer product <paramref name="v"/> * <paramref name="u"/>.</returns>
        /// <remarks>
        /// <para>The outer product of vectors is a matrix with elements M<sub>i,j</sub> = v<sub>i</sub>u<sub>j</sub>.</para>
        /// </remarks>
        public static Matrix operator * (ColumnVector v, RowVector u) {
            Matrix M = new Matrix(v.Dimension, u.Dimension);
            for (int r = 0; r < v.Dimension; r++) {
                for (int c = 0; c < u.Dimension; c++) {
                    M[r, c] = v[r] * u[c];
                }
            }
            return (M);
        }

        // multliplication by a matrix

        /// <summary>
        /// Computes the product of a row vector and a matrix.
        /// </summary>
        /// <param name="v">The row vector.</param>
        /// <param name="M">The matrix.</param>
        /// <returns>The product vM.</returns>
        public static RowVector operator * (RowVector v, IMatrix M) {
            if (v.Dimension != M.RowCount) throw new DimensionMismatchException();
            RowVector u = new RowVector(M.ColumnCount);
            for (int i = 0; i < u.Dimension; i++) {
                u[i] = 0.0;
                for (int j = 0; j < v.Dimension; j++) {
                    u[i] += v[j] * M[j, i];

                }
            }
            return (u);
        }

#if SHO
        /// <summary>
        /// Produces the representation of the vector for the Python interactive console.
        /// </summary>
        /// <returns>A string representation of the vector.</returns>        
        public string __repr__ () {
            StringWriter writer = new StringWriter();
            Matrix.WriteMatrix(this, writer);
            return (writer.ToString());
        }
#endif

    }

}
