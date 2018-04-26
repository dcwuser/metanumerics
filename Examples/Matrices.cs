using System;
using Meta.Numerics.Matrices;

namespace Examples {

    public static class Matrices {

        public static void PrintMatrix(string name, AnyRectangularMatrix M) {
            Console.WriteLine($"{name}=");
            for (int r = 0; r < M.RowCount; r++) {
                for (int c = 0; c < M.ColumnCount; c++) {
                    Console.Write("{0,10:g4}", M[r, c]);
                }
                Console.WriteLine();
            }
        }

        [ExampleMethod]
        public static void VectorsAndMatrices () {

            ColumnVector v = new ColumnVector(0.0, 1.0, 2.0);
            ColumnVector w = new ColumnVector(new double[] {1.0, -0.5, 1.5});
            
            SquareMatrix A = new SquareMatrix(new double[,] {
               {1, -2, 3},
               {2, -5, 12},
               {0, 2, -10}
            });

            RowVector u = new RowVector(4);
            for (int i = 0; i < u.Dimension; i++) u[i] = i;

            Random rng = new Random(1);
            RectangularMatrix B = new RectangularMatrix(4, 3);
            for (int r = 0; r < B.RowCount; r++) {
                for (int c = 0; c < B.ColumnCount; c++) {
                    B[r, c] = rng.NextDouble();
                }
            }

            SquareMatrix AI = A.Inverse();
            PrintMatrix("A * AI", A * AI);

            PrintMatrix("v + 2.0 * w", v + 2.0 * w);
            PrintMatrix("Av", A * v);
            PrintMatrix("B A", B * A);

            PrintMatrix("v^T", v.Transpose);
            PrintMatrix("B^T", B.Transpose);

            Console.WriteLine($"|v| = {v.Norm()}");
            Console.WriteLine($"sqrt(v^T v) = {Math.Sqrt(v.Transpose * v)}");

            UnitMatrix I = UnitMatrix.OfDimension(3);
            PrintMatrix("IA", I * A);

            Console.WriteLine(v == w);
            Console.WriteLine(I * A == A);


        }

        [ExampleMethod]
        public static void QRDecomposition () {

            SquareMatrix A = new SquareMatrix(new double[,] {
                { 1, -2, 3 },
                { 2, -5, 12 },
                { 0, 2, -10 }
            });

            ColumnVector b = new ColumnVector(2, 8, -4);

                        SquareQRDecomposition qrd = A.QRDecomposition();
            ColumnVector x = qrd.Solve(b);
            PrintMatrix("x", x);

            SquareMatrix Q = qrd.QMatrix;
            SquareMatrix R = qrd.RMatrix;
            PrintMatrix("QR", Q * R);

            PrintMatrix("Q  Q^T", Q.MultiplySelfByTranspose());

            SquareMatrix AI = qrd.Inverse();
            PrintMatrix("A^{-1}", AI);
            PrintMatrix("A^{-1} A", AI * A);

        }

        public static void LUDecomposition () {

            SquareMatrix A = new SquareMatrix(new double[,] {
                { 1, -2, 3 },
                { 2, -5, 12 },
                { 0, 2, -10 }
            });

            ColumnVector b = new ColumnVector(2, 8, -4);

            LUDecomposition lud = A.LUDecomposition();
            ColumnVector x = lud.Solve(b);
            PrintMatrix("x", x);
            PrintMatrix("Ax", A * x);

            SquareMatrix L = lud.LMatrix();
            SquareMatrix U = lud.UMatrix();
            SquareMatrix P = lud.PMatrix();
            PrintMatrix("LU", L * U);
            PrintMatrix("PA", P * A);

            SquareMatrix AI = lud.Inverse();
            PrintMatrix("A * AI", A * AI);

            Console.WriteLine($"det(a) = {lud.Determinant()}");

        }

        [ExampleMethod]
       public static void SingularValueDecomposition () {

            SquareMatrix A = new SquareMatrix(new double[,] {
                { 1, -2, 3 },
                { 2, -5, 12 },
                { 0, 2, -10 }
            });

            SingularValueDecomposition svdA = A.SingularValueDecomposition();
            Console.WriteLine($"A has rank = {svdA.Rank}, condition number = {svdA.ConditionNumber}.");

            SquareMatrix V = svdA.LeftTransformMatrix;
            SquareMatrix W = svdA.RightTransformMatrix;

            PrintMatrix("V^T V", V.MultiplyTransposeBySelf());
            PrintMatrix("W^T W", W.MultiplyTransposeBySelf());
            PrintMatrix("V^T A W", V.Transpose * A * W);

            ColumnVector b = new ColumnVector(0.0, 1.0, 2.0);
            ColumnVector x = svdA.Solve(b);
            PrintMatrix("Ax", A * x);

            SquareMatrix S = new SquareMatrix(new double[,] {
                { 11, 9, 0 },
                { 9, 11, 0 },
                { 0,  0, 0 }
            });
            SingularValueDecomposition svdS = S.SingularValueDecomposition();
            Console.WriteLine($"S has rank = {svdS.Rank}, condition number = {svdS.ConditionNumber}.");

            RectangularMatrix R = new RectangularMatrix(svdS.RowCount, svdS.ColumnCount);
            foreach (SingularValueContributor contributor in svdS.Contributors) {
                Console.WriteLine($"s * L * R^T with s = {contributor.SingularValue}");
                PrintMatrix("L^T", contributor.LeftSingularVector.Transpose);
                PrintMatrix("R^T", contributor.RightSingularVector.Transpose);
                R += contributor.SingularValue * contributor.LeftSingularVector * contributor.RightSingularVector.Transpose;
            }
            PrintMatrix("R", R);

            SingularValueContributorCollection contributors = svdS.Contributors;
            for (int i = 0; i < contributors.Count; i++) {
                for (int j = 0; j <= i; j++) {
                    Console.WriteLine($"l_{i}^T l_{j} = " +
                        $"{contributors[i].LeftSingularVector.Transpose * contributors[j].LeftSingularVector}");
                    Console.WriteLine($"r_{i}^T r_{j} = " +
                        $"{contributors[i].RightSingularVector.Transpose * contributors[j].RightSingularVector}");
                }
            }

            ColumnVector t = new ColumnVector(1.0, -1.0, 0.0);
            ColumnVector y = svdS.Solve(t);
            PrintMatrix("Sy", S * y);

        }

        [ExampleMethod]
        public static void EigenvaluesAndEigenvectors () {

            SymmetricMatrix H = new SymmetricMatrix(3);
            for (int r = 0; r < H.Dimension; r++) {
                for (int c = 0; c <= r; c++) {
                    H[r, c] = 1.0 / (r + c + 1);
                }
            }

            RealEigendecomposition ed = H.Eigendecomposition();

            SquareMatrix V = ed.TransformMatrix;
            PrintMatrix ("V^T V", V.Transpose * V);
            PrintMatrix("V D V^T", V * ed.DiagonalizedMatrix * V.Transpose);

            foreach(RealEigenpair pair in ed.Eigenpairs) {
                Console.WriteLine($"Eigenvalue {pair.Eigenvalue}");
                PrintMatrix("Hv", H * pair.Eigenvector);
                PrintMatrix("ev", pair.Eigenvalue * pair.Eigenvector);
            }

            ed.Eigenpairs.Sort(OrderBy.MagnitudeDescending);
            Console.WriteLine($"Largest eigenvalue {ed.Eigenpairs[0].Eigenvalue}");
            ed.Eigenpairs.Sort(OrderBy.ValueAscending);
            Console.WriteLine($"Least eigenvalue {ed.Eigenpairs[0].Eigenvalue}");

            double[] eigenvalues = H.Eigenvalues();

            double sum = 0.0;
            double product = 1.0;
            foreach (double eigenvalue in eigenvalues) {
                sum += eigenvalue;
                product *= eigenvalue;
            }
            Console.WriteLine($"sum(e) = {sum},  tr(H) = {H.Trace()}");
            Console.WriteLine($"prod(e) = {product}, det(H) = {H.CholeskyDecomposition().Determinant()}");

            SquareMatrix G1 = new SquareMatrix(4);
            G1[0, 3] = 1.0;
            G1[1, 2] = 1.0;
            G1[2, 1] = -1.0;
            G1[3, 0] = -1.0;

            ComplexEigendecomposition ced = G1.Eigendecomposition();
            foreach (ComplexEigenpair pair in ced.Eigenpairs) {
                Console.WriteLine(pair.Eigenvalue);
            }

        }

    }

}