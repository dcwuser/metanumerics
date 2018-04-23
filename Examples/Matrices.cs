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
       public static void SingularValueDecomposition () {

            SquareMatrix A = new SquareMatrix(new double[,] {
                { 1, -2, 3 },
                { 2, -5, 12 },
                { 0, 2, -10 }
            });

            SingularValueDecomposition SVD = A.SingularValueDecomposition();
            Console.WriteLine($"A has rank = {SVD.Rank}, condition number = {SVD.ConditionNumber}.");

            ColumnVector b = new ColumnVector(0.0, 1.0, 2.0);
            ColumnVector x = SVD.Solve(b);
            PrintMatrix("Ax", A * x);

            SquareMatrix S = new SquareMatrix(new double[,] {
                { 11, 9, 0 },
                { 9, 11, 0 },
                { 0,  0, 0 }
            });
            SingularValueDecomposition svdS = S.SingularValueDecomposition();
            Console.WriteLine($"D has rank = {svdS.Rank}, condition number = {svdS.ConditionNumber}.");

            RectangularMatrix R = new RectangularMatrix(svdS.RowCount, svdS.ColumnCount);
            foreach (SingularValueContributor contributor in svdS.Contributors) {
                Console.WriteLine($"s * L * R^T with s = {contributor.SingularValue}");
                PrintMatrix("L", contributor.LeftSingularVector);
                PrintMatrix("R^T", contributor.RightSingularVector.Transpose);
                R += contributor.SingularValue * contributor.LeftSingularVector * contributor.RightSingularVector.Transpose;
            }
            PrintMatrix("R", R);

            ColumnVector t = new ColumnVector(1.0, -1.0, 0.0);
            ColumnVector y = svdS.Solve(t);
            PrintMatrix("Sy", S * y);

        }

    }

}