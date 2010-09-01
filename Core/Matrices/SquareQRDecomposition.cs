using System;
using System.Collections.Generic;

namespace Meta.Numerics.Matrices {

    public sealed class SquareQRDecomposition : ISquareDecomposition {

        private double[] qtStore;
        private double[] rStore;
        private int dimension;

        internal SquareQRDecomposition (double[] qtStore, double[] rStore, int dimension) {
            this.qtStore = qtStore;
            this.rStore = rStore;
            this.dimension = dimension;
        }

        /// <summary>
        /// The orthogonal matrix Q.
        /// </summary>
        /// <returns>The orthogonal matrix Q.</returns>
        public SquareMatrix QMatrix () {
            double[] qStore = MatrixAlgorithms.Transpose(qtStore, dimension, dimension);
            return(new SquareMatrix(qStore, dimension));
        }

        /// <summary>
        /// The upper-right triangular matrix R.
        /// </summary>
        /// <returns>The upper-right triangular matrix R.</returns>
        public SquareMatrix RMatrix () {
            double[] store = MatrixAlgorithms.Clone(rStore, dimension, dimension);
            return (new SquareMatrix(store, dimension));
        }

        public ColumnVector Solve (IList<double> rhs) {

            if (rhs == null) throw new ArgumentNullException("rhs");
            if (rhs.Count != dimension) throw new DimensionMismatchException();

            double[] y = new double[rhs.Count];
            rhs.CopyTo(y, 0);

            y = MatrixAlgorithms.Multiply(qtStore, dimension, dimension, y, dimension, 1);

            return (new ColumnVector(y));
        }

        public double Determinant () {
            double det = 1.0;
            for (int i = 0; i < dimension; i++) {
                det *= MatrixAlgorithms.GetEntry(rStore, dimension, dimension, i, i);
            }
            return (det);
        }

        public SquareMatrix Inverse () {
            throw new NotImplementedException();
        }

        ISquareMatrix ISquareDecomposition.Inverse () {
            return (Inverse());
        }


        public int Dimension {
            get {
                return(dimension);
            }
        }
        
    }
}
