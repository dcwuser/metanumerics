using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Analysis {

#if FUTURE

    internal class Vertex {

        public double[] X;

        public double Y;

    }

    public static partial class MultiFunctionMath {

        const double alpha = 1.0;

        const double sigma = 0.5;

        public static void FindExtremum_Amobea (Func<IList<double>, double> function, IList<double> x0, EvaluationSettings settings) {
            MultiFunctor f = new MultiFunctor(function);

            int d = x0.Count;

            Vertex[] vertexes = new Vertex[d + 1];
            for (int i = 0; i < d; i++) {
                double[] x = new double[d];
                x0.CopyTo(x, 0);
                x[i] = x[i] + 1.0 + Math.Abs(x[i]);
                double y = f.Evaluate(x);
                vertexes[i] = new Vertex() { X = x, Y = y };
            }
            double[] x00 = new double[d];
            x0.CopyTo(x00, 0);
            double y00 = f.Evaluate(x00);
            vertexes[d] = new Vertex() { X = x00, Y = y00 };

            FindExtremum_Amobea(f, vertexes, settings);
        }

        private static void FindExtremum_Amobea (MultiFunctor f, Vertex[] vertexes, EvaluationSettings settings) {

            int d = vertexes.Length - 1;

            while (f.EvaluationCount < settings.EvaluationBudget) {

                // Identify the best and worst vertexes.
                int minVertex = 0; double minY = vertexes[0].Y;
                int maxVertex = 0; double maxY = vertexes[0].Y;
                int nextMaxVertex = 0; double nextMaxY = vertexes[0].Y;
                for (int i = 1; i < vertexes.Length; i++) {
                    double y = vertexes[i].Y;
                    if (y < minY) {
                        minVertex = i; minY = y;
                    }
                    if (y > nextMaxY) {
                        if (y > maxY) {
                            nextMaxVertex = maxVertex; nextMaxY = maxY;
                            maxVertex = i; maxY = y;
                        } else {
                            nextMaxVertex = i; nextMaxY = y;
                        }
                    }
                }

                // Terminate based on spread between vertexes.
                if ((maxY - minY) <= Math.Abs(maxY) * settings.RelativePrecision) {
                    Debug.WriteLine(minY);
                    return;
                }

                // Produce a new candidate vertex by reflecting the worst vertex through the opposite face.
                double[] centroid = new double[d];
                for (int i = 0; i < vertexes.Length; i++) {
                    if (i != maxVertex) {
                        for (int j = 0; j < d; j++) {
                            centroid[j] += vertexes[i].X[j] / d;
                        }
                    }
                }
                double[] newX = new double[d];
                for (int j = 0; j < d; j++) {
                    newX[j] = centroid[j] + alpha * (centroid[j] - vertexes[maxVertex].X[j]);
                }
                double newY = f.Evaluate(newX);

                if (newY < nextMaxY) {
                    // As long as the new point is not terrible, we are going to replace the worst point with it.
                    vertexes[maxVertex] = new Vertex() { X = newX, Y = newY };
                    Debug.WriteLine("Reflect");

                    if (newY < minY) {
                        // If the new point was very good, we will try to extend the simplex further in that direction.
                        double[] extendedX = new double[d];
                        for (int j = 0; j < d; j++) {
                            extendedX[j] = centroid[j] + 2.0 * (centroid[j] - vertexes[maxVertex].X[j]);
                        }
                        double extendedY = f.Evaluate(extendedX);
                        if (extendedY < minY) {
                            // If the extension is also very good, we replace the second worst point too.
                            vertexes[maxVertex] = new Vertex() { X = extendedX, Y = extendedY };
                            Debug.WriteLine("No, Extend");
                        }
                    }
                } else {
                    // The reflected point was pretty terrible, so we will try to produce a new candidate
                    // point by contracting the worst point toward the centroid instead.
                    for (int j = 0; j < d; j++) {
                        newX[j] = centroid[j] + (vertexes[maxVertex].X[j] - centroid[j]) / 2.0;
                    }
                    newY = f.Evaluate(newX);

                    if (newY < nextMaxY) {
                        // If that candidate is not terrible, accept it.
                        vertexes[maxVertex] = new Vertex() { X = newX, Y = newY };
                        Debug.WriteLine("Contract");
                    } else {
                        // Otherwise, we give up and simply shrink our simplex down toward the minimum.
                        for (int i = 0; i < vertexes.Length; i++) {
                            if (i != minVertex) {
                                double[] shrunkX = new double[d];
                                for (int j = 0; j < d; j++) {
                                    shrunkX[j] = vertexes[minVertex].X[j] + (vertexes[i].X[j] - vertexes[minVertex].X[j]) / 2.0;
                                }
                                double shrunkY = f.Evaluate(shrunkX);
                                vertexes[i] = new Vertex() { X = shrunkX, Y = shrunkY };
                            }
                        }
                        Debug.WriteLine("Shrink");
                    }


                }

            }

        }

    }

#endif

}
