using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Analysis {

    internal class DifferentialEvolutionSettings : MultiExtremumSettings {

        public double CrossoverProbability { get; set; }

        public int Population { get; set; }

    }

    public static partial class MultiFunctionMath {

        private const double mutationFactor = 0.625; //0.5;

        /// <summary>
        /// Finds the minimum of a function within the given volume.
        /// </summary>
        /// <param name="function">The function.</param>
        /// <param name="volume">The volume to search.</param>
        /// <returns>The global minimum.</returns>
        public static MultiExtremum FindGlobalMinimum (Func<IReadOnlyList<double>, double> function, IReadOnlyList<Interval> volume) {
            return (FindGlobalExtremum(function, volume, null, false));
        }

        private static DifferentialEvolutionSettings GetDefaultSettings (MultiExtremumSettings settings, int d) {

            if (settings == null) settings = new MultiExtremumSettings();

            DifferentialEvolutionSettings deSettings = new DifferentialEvolutionSettings();
            deSettings.Population = 8 * d + 4;
            deSettings.CrossoverProbability = 1.0 - 1.0 / 8.0 - 1.0 / d;

            deSettings.RelativePrecision = (settings.RelativePrecision < 0.0) ? Math.Pow(10.0, -(2.0 + 4.0 / d)) : settings.RelativePrecision;
            deSettings.AbsolutePrecision = (settings.AbsolutePrecision < 0.0) ? Math.Pow(10.0, -(4.0 + 8.0 / d)) : settings.AbsolutePrecision;
            deSettings.EvaluationBudget = (settings.EvaluationBudget < 0) ? 128 * d * d * d * d : settings.EvaluationBudget;

            deSettings.Listener = settings.Listener;

            return (deSettings);
        }

        /// <summary>
        /// Finds the minimum of a function within the given volume, subject to the given evaluation constraints.
        /// </summary>
        /// <param name="function">The function.</param>
        /// <param name="volume">The volume to search.</param>
        /// <param name="settings">The evaluation constraints to apply.</param>
        /// <returns>The global minimum.</returns>
        /// <remarks>
        /// <para>This algorithm attempts to find the global minimum of the given function within the entire given hyper-cube. It generally
        /// requires many more function evaluations than <see cref="FindGlobalMinimum(Func{IReadOnlyList{double}, double}, IReadOnlyList{Interval}, MultiExtremumSettings)"/>,
        /// but is much more likely to find a global minimum in situations where multiple local minima exist.</para>
        /// <para>Unlike <see cref="FindGlobalMinimum(Func{IReadOnlyList{double}, double}, IReadOnlyList{Interval}, MultiExtremumSettings)"/>,
        /// this method does not return an approximate Hessian matrix near the minimum.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="function"/> or <paramref name="volume"/> is <see langword="null"/>.</exception>
        /// <exception cref="NonconvergenceException">The minimum could not be found within the given evaluation budget.</exception>
        public static MultiExtremum FindGlobalMinimum (Func<IReadOnlyList<double>, double> function, IReadOnlyList<Interval> volume, MultiExtremumSettings settings) {
            return (FindGlobalExtremum(function, volume, settings, false));
        }

        /// <summary>
        /// Finds the maximum of a function within the given volume.
        /// </summary>
        /// <param name="function">The function.</param>
        /// <param name="volume">The volume to search.</param>
        /// <returns>The global maximum.</returns>
        public static MultiExtremum FindGlobalMaximum (Func<IReadOnlyList<double>, double> function, IReadOnlyList<Interval> volume) {
            return (FindGlobalExtremum(function, volume, null, true));
        }

        /// <summary>
        /// Finds the maximum of a function within the given volume, subject to the given evaluation constraints.
        /// </summary>
        /// <param name="function">The function.</param>
        /// <param name="volume">The volume to search.</param>
        /// <param name="settings">The evaluation constraints to apply.</param>
        /// <returns>The global maximum.</returns>
        public static MultiExtremum FindGlobalMaximum (Func<IReadOnlyList<double>, double> function, IReadOnlyList<Interval> volume, MultiExtremumSettings settings) {
            return (FindGlobalExtremum(function, volume, settings, true));
        }

        private static MultiExtremum FindGlobalExtremum (Func<IReadOnlyList<double>, double> function, IReadOnlyList<Interval> volume, MultiExtremumSettings settings, bool negate) {
            if (function == null) throw new ArgumentNullException(nameof(function));
            if (volume == null) throw new ArgumentNullException(nameof(volume));
            MultiFunctor f = new MultiFunctor(function, negate);
            DifferentialEvolutionSettings deSettings = GetDefaultSettings(settings, volume.Count);
            MultiExtremum extremum = FindGlobalExtremum(f, volume, deSettings);
            return (extremum);
        }

        // Differential evolution is a global optimization algorithm over continuous inputs that is adapted from genetic algorithms for finite inputs.
        // The idea is to maintain a population of input vectors ("agents") and to vary that population over cycles ("generations") according to rules that incorporate
        // random mutation but on average tend to bring them closer to optima ("fitter").

        private static MultiExtremum FindGlobalExtremum (MultiFunctor f, IReadOnlyList<Interval> volume, DifferentialEvolutionSettings settings) {

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


            while (f.EvaluationCount < settings.EvaluationBudget) {
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
                double tol = settings.ComputePrecision(minValue);
                if (range <= tol) {
                    MultiExtremum result = new MultiExtremum(f.EvaluationCount, settings, points[minIndex], f.IsNegated ? -values[minIndex] : values[minIndex], Math.Max(range, 0.75 * tol), null);
                    return (result);
                } else if (settings.Listener != null) {
                    MultiExtremum report = new MultiExtremum(f.EvaluationCount, settings, points[minIndex], f.IsNegated ? -values[minIndex] : values[minIndex], Math.Max(range, 0.75 * tol), null);
                    settings.Listener(report);
                }

            }

            throw new NonconvergenceException();

        }

    }


}
