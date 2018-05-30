using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Meta.Numerics.Matrices {

    // we use the storage method described in Knuth
    // at five values per entry, this is not the most storage-efficient format
    // but it is conceptually straightforward and allows easy update
    // other formats i have seen, e.g. Boeing-Harwell, require that the matrix entries be supplied in a particular, non-random order
    // and do not allow easy update

    internal class SparseMatrixElement {

        public SparseMatrixElement (int row, int column, double value) {
            Row = row;
            Column = column;
            Value = value;
        }

        public int Row;
        public int Column;
        public double Value;
        public SparseMatrixElement NextInRow;
        public SparseMatrixElement NextInColumn;


    }

    /// <summary>
    /// Represents a sparse, square matrix.
    /// </summary>
    /// <remarks>
    /// <para>Many applications give rise to very large matrices which consist mostly of zero elements.</para>
    /// <para>When working with sparse matrices, it is important to keep in mind that many operations do not respect sparsity.
    /// For example, the product of two sparse matrices is not necessarily sparse, nor is the inverse of a sparse matrix.</para>
    /// </remarks>
    public sealed class SparseSquareMatrix : AnySquareMatrix {

        /// <summary>
        /// Initializes a new sparse, square matrix.
        /// </summary>
        /// <param name="dimension">The dimension of the matrix, which must be positive.</param>
        public SparseSquareMatrix (int dimension) {
            if (dimension < 1) throw new ArgumentOutOfRangeException(nameof(dimension));
            this.dimension = dimension;
            rows = new SparseMatrixElement[dimension];
            columns = new SparseMatrixElement[dimension];
            fill = 0;
        }

        internal SparseSquareMatrix (int dimension, SparseMatrixElement[] rows, SparseMatrixElement[] columns, int fill) {
            this.dimension = dimension;
            this.rows = rows;
            this.columns = columns;
            this.fill = fill;
        }

        private int dimension;
        private SparseMatrixElement[] rows;
        private SparseMatrixElement[] columns;
        private int fill;

        /// <inheritdoc />
        public override int Dimension {
            get {
                return (dimension);
            }
        }

        /// <inheritdoc />
        public override double this[int r, int c] {
            get {
                if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException(nameof(r));
                if ((c < 0) || (c >= dimension)) throw new ArgumentOutOfRangeException(nameof(c));
                SparseMatrixElement element = GetElement(r, c);
                if (element == null) {
                    return(0.0);
                } else {
                    return(element.Value);
                }
            }
            set {
                if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException(nameof(r));
                if ((c < 0) || (c >= dimension)) throw new ArgumentOutOfRangeException(nameof(c));
                SetElement(r, c, value);
            }
        }

        private SparseMatrixElement GetElement (int r, int c) {
            SparseMatrixElement current = rows[r];
            while (current != null) {
                if (current.Column == c) return (current);
                current = current.NextInRow;
            }
            return (null);
        }

        private SparseMatrixElement SetElement (int r, int c, double value) {

            // if the row is empty: previousInRow == null, nextInRow == null
            // if the new element is the leftmost in the row: previousInRow == null, nextInRow != null
            // if the new element is intermediate in the row: previousInRow != null, nextInRow != null
            // if the new element is the rightmost in the row: previousInRow != null, nextInRow == null
            // if the element already exists: nextInRow is the element element, previousInRow is null if it is the first in the row
            SparseMatrixElement previousInRow = null;
            SparseMatrixElement nextInRow = rows[r];
            while ((nextInRow != null) && (nextInRow.Column < c)) {
                previousInRow = nextInRow;
                nextInRow = nextInRow.NextInRow;
            }

            SparseMatrixElement previousInColumn = null;
            SparseMatrixElement nextInColumn = columns[c];
            while ((nextInColumn != null) && (nextInColumn.Row < r)) {
                previousInColumn = nextInColumn;
                nextInColumn = nextInColumn.NextInColumn;
            }

            if ((nextInRow != null) && (nextInRow.Row == r) && (nextInRow.Column == c)) {
                // the element already exists
                if (value != 0.0) {
                    // if the value is non-zero, just update it
                    nextInRow.Value = value;
                    return (nextInRow);
                } else {
                    // if the value is zero, delete the element
                    if (previousInRow == null) {
                        rows[r] = nextInRow.NextInRow;
                    } else {
                        previousInRow.NextInRow = nextInRow.NextInRow;
                    }
                    if (previousInColumn == null) {
                        columns[c] = nextInColumn.NextInColumn;
                    } else {
                        previousInColumn.NextInColumn = nextInColumn.NextInColumn;
                    }
                    fill--;
                    return (null);
                }
            } else {
                // the element does not yet exist
                if (value == 0.0) {
                    // if the value is zero, we don't need to do anything
                    return (null);
                } else {
                    // otherwise, we need to create an element
                    SparseMatrixElement element = new SparseMatrixElement(r, c, value);
                    if (previousInRow == null) {
                        rows[r] = element;
                    } else {
                        previousInRow.NextInRow = element;
                    }
                    element.NextInRow = nextInRow;
                    if (previousInColumn == null) {
                        columns[c] = element;
                    } else {
                        previousInColumn.NextInColumn = element;
                    }
                    element.NextInColumn = nextInColumn;
                    fill++;
                    return (element);
                }
            }

            
        }

        /// <summary>
        /// Gets the number of non-zero matrix entries.
        /// </summary>
        public int FillCount {
            get {
                return (fill);
            }
        }

        /// <summary>
        /// Gets the fraction of matrix entries that are non-zero.
        /// </summary>
        public double FillFraction {
            get {
                return (1.0 * fill / dimension / dimension);
            }
        }

        /// <inheritdoc />
        public override RowVector Row (int r) {
            if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException(nameof(r));
            double[] store = new double[dimension];
            SparseMatrixElement element = rows[r];
            while (element != null) {
                store[element.Column] = element.Value;
                element = element.NextInRow;
            }
            return (new RowVector(store, dimension));
        }

        /// <inheritdoc />
        public override ColumnVector Column (int c) {
            if ((c < 0) || (c >= dimension)) throw new ArgumentOutOfRangeException(nameof(c));
            double[] store = new double[dimension];
            SparseMatrixElement element = columns[c];
            while (element != null) {
                store[element.Row] = element.Value;
                element = element.NextInColumn;
            }
            return (new ColumnVector(store, dimension));
        }

        /// <summary>
        /// Copies the matrix.
        /// </summary>
        /// <returns>An independent copy of the matrix.</returns>
        public SparseSquareMatrix Copy () {

            // a simple, slow solution is to create an empty sparse matrix with the same dimension and call
            // set element for each entry; but this moves along each linked list many times

            // we would like to just move along, say, each row, just once, copying the elements we find
            // but that isn't simple, because we need to store pointers to the next in row and next in column,
            // which we don't have on hand as we create the element

            // dealing with the next in row is easy: just keep a reference to the last one we do and
            // set its pointer to the next in row when we create the next element in the row

            // dealing with the next in column is harder; we need to know what is above it in the column,
            // but that was created a long time ago when we were dealing with a previous row.
            // to solve this problem, we use auxiluary storage: a N-element array that stores that last
            // element created in each column. when we create a new element, we hook it up to the previous
            // element stored in that array, then put the new element in the array


            SparseMatrixElement[] copyRows = new SparseMatrixElement[dimension];
            SparseMatrixElement[] copyColumns = new SparseMatrixElement[dimension];
            SparseMatrixElement[] lastInColumn = new SparseMatrixElement[dimension];

            for (int r = 0; r < dimension; r++) {
                SparseMatrixElement element = rows[r];
                SparseMatrixElement lastInRow = null;
                while (element != null) {
                    // create a copy of the element
                    SparseMatrixElement copyElement = new SparseMatrixElement(element.Row, element.Column, element.Value);
                    // hook it up to the previous one in the row (and store it for the next one)
                    if (lastInRow != null) {
                        lastInRow.NextInRow = copyElement;
                    } else {
                        copyRows[r] = copyElement;
                    }
                    lastInRow = copyElement;
                    // hook it up to the previous one in the column (and store it for the next one)
                    if (lastInColumn[element.Column] != null) {
                        lastInColumn[element.Column].NextInColumn = copyElement;
                    } else {
                        copyColumns[element.Column] = copyElement;
                    }
                    lastInColumn[element.Column] = copyElement;
                    // move to the next element in the row
                    element = element.NextInRow;
                }
            }

            SparseSquareMatrix copy = new SparseSquareMatrix(dimension, copyRows, copyColumns, fill);
            return (copy);

        }

        /// <summary>
        /// Multiplies a sparse matrix by a real scalar.
        /// </summary>
        /// <param name="alpha">The scalar value.</param>
        /// <param name="A">The sparse matrix.</param>
        /// <returns>The product sparse matrix.</returns>
        public static SparseSquareMatrix operator * (double alpha, SparseSquareMatrix A) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            SparseSquareMatrix aA = A.Copy();
            if (alpha == 0.0) return (aA);
            for (int r = 0; r < aA.dimension; r++) {
                SparseMatrixElement element = aA.rows[r];
                while (element != null) {
                    element.Value = alpha * element.Value;
                    // handle the case where alpha != 0, but alpha * element.Value underflows to zero
                    element = element.NextInRow;
                }
            }
            return(aA);
        }

        /// <summary>
        /// Multiplies a column vector by a sparse matrix.
        /// </summary>
        /// <param name="A">The matrix.</param>
        /// <param name="v">The column vector.</param>
        /// <returns>The product column vector.</returns>
        public static ColumnVector operator * (SparseSquareMatrix A, ColumnVector v) {

           if (A == null) throw new ArgumentNullException(nameof(A));
           if (v == null) throw new ArgumentNullException(nameof(v));
           if (A.Dimension != v.Dimension) throw new DimensionMismatchException();

           ColumnVector Av = new ColumnVector(A.Dimension);
           for (int i = 0; i < A.Dimension; i++) {
               SparseMatrixElement element = A.rows[i];
               while (element != null) {
                   Av[i] += element.Value * v[element.Column];
                   element = element.NextInRow;
               }
           }
           return (Av);

       }


        /// <summary>
        /// Multiplies a sparse matrix by a row vector.
        /// </summary>
        /// <param name="A">The matrix.</param>
        /// <param name="v">The row vector.</param>
        /// <returns>The product row vector.</returns>
        public static RowVector operator * (RowVector v, SparseSquareMatrix A) {

            if (v == null) throw new ArgumentNullException(nameof(v));
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (v.Dimension != A.Dimension) throw new DimensionMismatchException();

            RowVector vA = new RowVector(A.Dimension);
            for (int i = 0; i < A.Dimension; i++) {
                SparseMatrixElement element = A.columns[i];
                while (element != null) {
                    vA[i] += element.Value * v[element.Row];
                    element = element.NextInColumn;
                }
            }
            return (vA);

        }

        /// <summary>
        /// Solves Ax = b using iterative methods.
        /// </summary>
        /// <param name="rhs">The right-hand side vector b.</param>
        /// <returns>The solution vector x.</returns>
        /// <remarks>
        /// <para>In general, neither the inverse nor any decomposition of a sparse matrix is itself sparse. Therefore,
        /// to solve large, sparse linear systems, iterative methods are employed. An iterative method begins with
        /// an approximate or guessed solution vector and progresses toward an improved solution. Iterative methods
        /// are often successful at converging to a sufficiently accurate solution vector, but this is not guaranteed.
        /// If this method fails to converge, it throws a <see cref="NonconvergenceException"/>.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="rhs"/> is null.</exception>
        /// <exception cref="DimensionMismatchException"><paramref name="rhs"/>'s dimension does not equal the matrix's dimension.</exception>
        /// <exception cref="NonconvergenceException">The method did not converge to a solution.</exception>
        public ColumnVector Solve (ColumnVector rhs) {

            if (rhs == null) throw new ArgumentNullException(nameof(rhs));
            if (rhs.Dimension != dimension) throw new DimensionMismatchException();

            // accuracy and iteration limits are very huristic, revisit later
            double accuracyGoal = Global.Accuracy * 128.0;
            int iterationMax = 256 + dimension / 2;

            // we use the stabilized biconjugate gradient algorithm

            // choose an initial guess
            ColumnVector x = new ColumnVector(dimension);
            for (int i = 0; i < x.Dimension; i++) x[i] = 1.0;

            // 
            RowVector rt = x.Transpose;

            // r is the deviation vector that we are trying to drive to zero
            ColumnVector r = rhs - this * x;

            double rho0 = 1.0;
            double a = 1.0;
            double omega = 1.0;

            ColumnVector v = new ColumnVector(dimension);
            ColumnVector p = new ColumnVector(dimension);

            for (int i = 1; i < iterationMax; i++) {

                double rho1 = rt * r;
                double beta = (rho1 / rho0) * (a / omega);
                p = r + beta * (p - omega * v);
                v = this * p;
                a = rho1 / (rt * v);
                ColumnVector s = r - a * v;
                ColumnVector t = this * s;

                omega = (t.Transpose * s) / (t.Transpose * t);

                x = x + a * p + omega * s;

                r = s - omega * t;

                if (r.FrobeniusNorm() <= accuracyGoal * rhs.FrobeniusNorm()) return (x);

                // prepare for next iteration
                rho0 = rho1;
            }

            throw new NonconvergenceException();

        }
    }

}
