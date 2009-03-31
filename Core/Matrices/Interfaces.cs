using System;
using System.Collections.Generic;
using System.Text;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// The contract fufilled by any real matrix.
    /// </summary>
    public interface IMatrix {

        /// <summary>
        /// Gets the number of rows in the matrix.
        /// </summary>
        int RowCount { get; }

        /// <summary>
        /// Gets the number of columns in the matrix.
        /// </summary>
        int ColumnCount { get; }

        /// <summary>
        /// Gets the value of a matrix entry.
        /// </summary>
        /// <param name="r">The (zero-based) row number.</param>
        /// <param name="c">The (zero-based) column number.</param>
        /// <returns>The value of the matrix entry.</returns>
        double this[int r, int c] { get; }

        /// <summary>
        /// Returns a clone of the matrix.
        /// </summary>
        /// <returns></returns>
        IMatrix Clone ();

        /// <summary>
        /// Returns a transpose of the matrix.
        /// </summary>
        /// <returns></returns>
        IMatrix Transpose ();

    }

    /// <summary>
    /// The contract fufilled by any real, square matrix.
    /// </summary>
    public interface ISquareMatrix : IMatrix {

        /// <summary>
        /// Gets the dimension of the matrix.
        /// </summary>
        int Dimension { get; }

        /// <summary>
        /// Computes the trace of the matrix.
        /// </summary>
        /// <returns></returns>
        double Trace ();

        //ISquareDecomposition Decomposition ();

    }


    /// <summary>
    /// The contract fufilled by any real, symmetric matrix.
    /// </summary>
    public interface ISymmetricMatrix : ISquareMatrix {
    }


    /// <summary>
    /// The contract fufilled by a real LU decomposition.
    /// </summary>
    public interface ISquareDecomposition {

        /// <summary>
        /// Gets the dimension of the system.
        /// </summary>
        int Dimension { get; }

        /// <summary>
        /// Computes the solution to a system of equations.
        /// </summary>
        /// <param name="rhs"></param>
        /// <returns></returns>
        ColumnVector Solve (IList<double> rhs);

        /// <summary>
        /// Computes the inverse of the original matrix.
        /// </summary>
        /// <returns></returns>
        ISquareMatrix Inverse ();

        /// <summary>
        /// Computes the determinant of the original matrix.
        /// </summary>
        /// <returns></returns>
        double Determinant ();

    }

}
