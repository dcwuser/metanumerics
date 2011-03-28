using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Meta.Numerics.Matrices {

    // we use the storage method described in Knuth
    // at five values per entry, this is not the most storage-efficient format
    // but it is conceptually straightforward and allows easy update
    // other formats i have seen, e.g. Boeing-Harwell, require that the matrix entries be supplied in a particular, non-random order
    // and do nog allow easy update

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
    /// </remarks>
    public sealed class SparseSquareMatrix : AnySquareMatrix {

        /// <summary>
        /// Initializes a new sparse, square matrix.
        /// </summary>
        /// <param name="dimension">The dimension of the matrix, which must be positive.</param>
        public SparseSquareMatrix (int dimension) {
            if (dimension < 1) throw new ArgumentOutOfRangeException("dimension");
            this.dimension = dimension;
            rows = new SparseMatrixElement[dimension];
            columns = new SparseMatrixElement[dimension];
            fill = 0;
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
                if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException("r");
                if ((c < 0) || (c >= dimension)) throw new ArgumentOutOfRangeException("c");
                SparseMatrixElement element = GetElement(r, c);
                if (element == null) {
                    return(0.0);
                } else {
                    return(element.Value);
                }
            }
            set {
                if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException("r");
                if ((c < 0) || (c >= dimension)) throw new ArgumentOutOfRangeException("c");
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
        public override ColumnVector Column (int c) {
            if ((c < 0) || (c >= dimension)) throw new ArgumentOutOfRangeException("c");
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
            // this is not a very fast copy, it would be nice to improve it
            SparseSquareMatrix copy = new SparseSquareMatrix(dimension);
            for (int r = 0; r < dimension; r++) {
                SparseMatrixElement element = rows[r];
                while (element != null) {
                    copy[element.Row, element.Column] = element.Value;
                    element = element.NextInRow;
                }
            }
            return (copy);
        }

        public static SparseSquareMatrix operator * (double alpha, SparseSquareMatrix A) {

            if (A == null) throw new ArgumentNullException("A");
            SparseSquareMatrix aA = new SparseSquareMatrix(A.dimension);
            for (int r = 0; r < A.dimension; r++) {
                SparseMatrixElement element = A.rows[r];
                while (element != null) {

                }
            }
            throw new NotImplementedException();
        }

        /// <summary>
        /// Multiplies a column vector by a sparse matrix.
        /// </summary>
        /// <param name="A">The matrix.</param>
        /// <param name="v">The column vector.</param>
        /// <returns>The product column vector.</returns>
        public static ColumnVector operator * (SparseSquareMatrix A, ColumnVector v) {

           if (A == null) throw new ArgumentNullException("A");
           if (v == null) throw new ArgumentNullException("v");
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
        /// Multiplies a sparse matrix by a row vector..
        /// </summary>
        /// <param name="A">The matrix.</param>
        /// <param name="v">The row vector.</param>
        /// <returns>The product row vector.</returns>
        public static RowVector operator * (RowVector v, SparseSquareMatrix A) {

            if (v == null) throw new ArgumentNullException("v");
            if (A == null) throw new ArgumentNullException("A");
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
    }

}
