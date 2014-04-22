using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Functions {

    public static partial class MultiFunctionMath {

        private const double crossoverProbability = 0.625; //0.5;
        private const double mutationFactor = 0.625; //0.5;

        public static MultiExtremum FindGlobalMinimum (Func<IList<double>, double> function, IList<Interval> volume) {
            if (function == null) throw new ArgumentNullException("function");
            if (volume == null) throw new ArgumentNullException("volume");
            int d = volume.Count;
            EvaluationSettings settings = new EvaluationSettings() {
                RelativePrecision = 1.0E-4,
                AbsolutePrecision = 1.0E-6,
                EvaluationBudget = (8 * d) * (8 * d) * (8 * d)
            };
            return (FindGlobalMinimum(function, volume, settings));
        }

        public static MultiExtremum FindGlobalMinimum (Func<IList<double>, double> function, IList<Interval> volume, EvaluationSettings settings) {
            MultiFunctor f = new MultiFunctor(function);
            return (FindGlobalMinimum(f, volume, settings));
        }

        private static MultiExtremum FindGlobalMinimum (MultiFunctor f, IList<Interval> volume, EvaluationSettings settings) {

            int d = volume.Count;

            // Choose a number of agents that increases with dimension and required precision.
            int m = 8 * d;
            Debug.WriteLine("d={0} m={1}", d, m);

            Random rng = new Random(1);

            double[][] points = new double[m][];
            double[] values = new double[m];
            for (int i = 0; i < m; i++) {
                points[i] = new double[d];
                for (int j = 0; j < d; j++) points[i][j] = volume[j].LeftEndpoint + rng.NextDouble() * volume[j].Width;
                values[i] = f.Evaluate(points[i]);
            }

            while (f.EvaluationCount < settings.EvaluationBudget) {

                //double mutationFactor = 0.5 + 0.5 * rng.NextDouble();

                double[][] newPoints = new double[m][];
                double[] newValues = new double[m];

                for (int i = 0; i < m; i++) {
                    
                    // Mutation
                    // construct donor vector
                    int a = i; while (a == i) a = rng.Next(m);
                    int b = i; while ((b == i) || (b == a)) b = rng.Next(m);
                    int c = i; while ((c == i) || (c == b) || (c == a)) c = rng.Next(m);
                    double[] donor = new double[d];
                    for (int j = 0; j < d; j++) {
                        donor[j] = points[a][j] + mutationFactor * (points[b][j] - points[c][j]);
                        if (donor[j] < volume[j].LeftEndpoint) donor[j] = volume[j].LeftEndpoint;
                        if (donor[j] > volume[j].RightEndpoint) donor[j] = volume[j].RightEndpoint;
                    }

                    // Recombination
                    double[] trial = new double[d];
                    int k = rng.Next(d);
                    for (int j = 0; j < d; j++) {
                        if ((j == k) || (rng.NextDouble() < crossoverProbability)) {
                            trial[j] = donor[j];
                        } else {
                            trial[j] = points[i][j];
                        }
                    }

                    // Selection
                    double value = f.Evaluate(trial);
                    if (value <= values[i]) {
                        newPoints[i] = trial;
                        newValues[i] = value;
                    } else {
                        newPoints[i] = points[i];
                        newValues[i] = values[i];
                    }

                }

                points = newPoints;
                values = newValues;

                // Check termination criteria
                int minIndex = -1;
                double minValue = Double.MaxValue;
                double maxValue = Double.MinValue;
                for (int i = 0; i < m; i++) {
                    if (values[i] < minValue) {
                        minValue = values[i];
                        minIndex = i;
                    }
                    if (values[i] > maxValue) maxValue = values[i];
                }
                double range = maxValue - minValue;
                if ((range <= settings.AbsolutePrecision) || (range <= Math.Abs(maxValue) * settings.RelativePrecision)) {
                    MultiExtremum result = new MultiExtremum(points[minIndex], values[minIndex], null, f.EvaluationCount);
                    return (result);
                }

                settings.OnUpdate(new MultiExtremum(points[minIndex], values[minIndex], null, f.EvaluationCount));

            }

            throw new NonconvergenceException();

        }

    }


}
