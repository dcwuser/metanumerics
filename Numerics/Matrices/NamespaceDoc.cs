using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace Meta.Numerics.Matrices {


    /// <summary>
    /// Contains types that store and operate on matrices.
    /// </summary>
    /// <remarks>
    /// <para>You can create matrices of real values using the <see cref="RectangularMatrix"/> and <see cref="SquareMatrix"/>
    /// classes, and vectors of real values using the <see cref="ColumnVector"/> and <see cref="RowVector"/> classes.
    /// Operator overloads are defined that allow you to perform allowed arithmetic operations, such as adding
    /// two vectors or matrices, or multiplying a vector by a scalar, vector, or matrix. Each type defines
    /// methods corresponding to common linear algebra operations, such as inversion (<see cref="SquareMatrix.Inverse"/>),
    /// finding eigenvalues and eigenvectors (<see cref="SquareMatrix.Eigenvalues"/> and <see cref="SquareMatrix.Eigendecomposition"/>),
    /// and decompositions (<see cref="RectangularMatrix.QRDecomposition"/> and <see cref="RectangularMatrix.SingularValueDecomposition"/>).</para>
    /// <para>The fastest way to solve a linear system A x = b is to form the <see cref="SquareMatrix.LUDecomposition"/> of A
    /// and call <see cref="LUDecomposition.Solve(IReadOnlyList{double})"/> with the right-hand-side b.</para>
    /// <para>There are several additional matrix containers that support smaller storage requirements and faster operations for
    /// matrices with particular structures, such as <see cref="DiagonalMatrix"/>, <see cref="TridiagonalMatrix"/>,
    /// <see cref="SymmetricMatrix"/>, and <see cref="SparseSquareMatrix"/>.</para>
    /// <para>Where possible, we quickly return new matrix objects that implement a new view of existing stored values,
    /// without copying or otherwise disturbing the original values. Examples include <see cref="SquareMatrix.Transpose"/> and <see cref="SquareMatrix.Row(int)"/>.
    /// For read-only purposes, this is much faster and requires less memory that computing and storing new values. The returned matrix objects are,
    /// however, necessarily read-only. Whether an matrix object is read-only can be determined from <see cref="AnyMatrix{T}.IsReadOnly"/>. If
    /// you want to modify a read-only matrix object, simply make a copy using the object's Copy method (e.g. <see cref="SquareMatrix.Copy"/> or
    /// <see cref="RowVector.Copy"/>).</para>
    /// </remarks>
    [CompilerGenerated]
    internal class NamespaceDoc {
    }
}
