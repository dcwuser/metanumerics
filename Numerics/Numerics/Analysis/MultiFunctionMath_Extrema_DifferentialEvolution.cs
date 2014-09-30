using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Analysis {

    internal class DifferentialEvolutionSettings : EvaluationSettings {

        public double CrossoverProbability { get; set; }

        public int Population { get; set; }

    }

    public static partial class MultiFunctionMath {

        //private const double crossoverProbability =  1.0 - 1.0 / 8.0 - 1.0 / 14.0; //0.9286; //0.625; //0.5;
        private const double mutationFactor = 0.625; //0.5;

        public static MultiExtremum FindGlobalMinimum (Func<IList<double>, double> function, IList<Interval> volume) {
            if (function == null) throw new ArgumentNullException("function");
            if (volume == null) throw new ArgumentNullException("volume");
            return (FindGlobalMinimum(function, volume, null));
        }

        private static DifferentialEvolutionSettings GetDefaultSettings (EvaluationSettings settings, int d) {

            DifferentialEvolutionSettings deSettings = new DifferentialEvolutionSettings();
            deSettings.Population = 8 * d + 4;
            deSettings.CrossoverProbability = 1.0 - 1.0 / 8.0 - 1.0 / d;

            if (settings == null) {
                deSettings.RelativePrecision = 1.0E-4;
                deSettings.AbsolutePrecision = 1.0E-8;
                deSettings.EvaluationBudget = 128 * d * d * d * d;
            } else {
                deSettings.RelativePrecision = settings.RelativePrecision;
                deSettings.AbsolutePrecision = settings.AbsolutePrecision;
                deSettings.EvaluationBudget = settings.EvaluationBudget;
            }

            return (deSettings);
        }

        public static MultiExtremum FindGlobalMinimum (Func<IList<double>, double> function, IList<Interval> volume, EvaluationSettings settings) {
            if (function == null) throw new ArgumentNullException("function");
            if (volume == null) throw new ArgumentNullException("volume");
            DifferentialEvolutionSettings deSettings = GetDefaultSettings(settings, volume.Count);
            MultiFunctor f = new MultiFunctor(function);
            return (FindGlobalExtremum(f, volume, deSettings));
        }

        public static MultiExtremum FindGlobalMaximum (Func<IList<double>, double> function, IList<Interval> volume) {
            if (function == null) throw new ArgumentNullException("function");
            if (volume == null) throw new ArgumentNullException("volume");
            return (FindGlobalMaximum(function, volume, null));
        }

        public static MultiExtremum FindGlobalMaximum (Func<IList<double>, double> function, IList<Interval> volume, EvaluationSettings settings) {
            if (function == null) throw new ArgumentNullException("function");
            if (volume == null) throw new ArgumentNullException("volume");
            DifferentialEvolutionSettings deSettings = GetDefaultSettings(settings, volume.Count);
            MultiFunctor f = new MultiFunctor(function, true);
            MultiExtremum maximum = FindGlobalExtremum(f, volume, deSettings);
            return (maximum);
        }

        private static MultiExtremum FindGlobalExtremum (MultiFunctor f, IList<Interval> volume, DifferentialEvolutionSettings settings) {

            int d = volume.Count;

            // Choose a number of agents that increases with dimension and required precision.
            int m = settings.Population;
            Debug.WriteLine("d={0} m={1}", d, m);

            Random rng = new Random(3);

            // Start with random points in the allowed region.
            double[][] points = new double[m][];
            double[] values = new double[m];
            for (int i = 0; i < m; i++) {
                points[i] = new double[d];
                for (int j = 0; j < d; j++) {
                    points[i][j] = volume[j].LeftEndpoint + rng.NextDouble() * volume[j].Width;
                }
                values[i] = f.Evaluate(points[i]);
            }

            //int[] indexes = new int[m];
            //for (int i = 0; i < indexes.Length; i++) indexes[i] = i;

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
                        if ((j == k) || (rng.NextDouble() < settings.CrossoverProbability)) {
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
                    MultiExtremum result = new MultiExtremum(f.EvaluationCount, settings, points[minIndex], f.IsNegated ? -values[minIndex] : values[minIndex], null);
                    return (result);
                }

                //settings.OnUpdate(new MultiExtremum(points[minIndex], values[minIndex], null, f.EvaluationCount));

            }

            throw new NonconvergenceException();

        }

    }


}
