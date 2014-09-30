using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Functions;
using Meta.Numerics.Statistics;

namespace Meta.Numerics.Analysis {

    // Interface for Monte Carlo point generators. Used to represent both pseudo-random and quasi-random point sequences.

    internal abstract class VectorGenerator {

        protected VectorGenerator (int dimension) {
            Debug.Assert(dimension > 0);
            this.dimension = dimension;
        }

        private int dimension;

        public int Dimension {
            get {
                return(dimension);
            }
        }

        public abstract double[] NextVector ();

    }

    internal class RandomVectorGenerator : VectorGenerator {

        public RandomVectorGenerator (int d, Random rng) : base(d) {
            this.rng = rng;
        }

        private Random rng;

        public override double[] NextVector () {
            double[] v = new double[Dimension];
            for (int i = 0; i < v.Length; i++) {
                v[i] = rng.NextDouble();
            }
            return (v);
        }

    }

    internal class SobolVectorGenerator : VectorGenerator {

        public SobolVectorGenerator (int d) : base(d) {
            sequences = new SobolSequence[d];
            for (int i = 0; i < sequences.Length; i++) {
                SobolSequenceParameters p = FunctionMath.sobolParameters[i];
                sequences[i] = new SobolSequence(p.Dimension, p.Coefficients, p.Seeds);
            }
        }

        private SobolSequence[] sequences;

        public override double[] NextVector () {
            double[] x = new double[sequences.Length];
            for (int i = 0; i < sequences.Length; i++) {
                sequences[i].MoveNext();
                x[i] = sequences[i].Current;
            }
            return (x);
        }
    }

    // LePage's method of importance sampling is described in "A New Algorithm for Adaptive Multidimensional Integration" 1977

    // The basic idea of importance sampling is to write
    //   \int dV f(x) = \int dV p(x) \frac{f(x)}{p(x)}
    // Intrepret p(x) as a probabily density to choose point x and sample f(x)/p(x) instead of f(x). If p(x) \propto f(x) then we are
    // just sampling a constant.

    // Approximating p(x) via a grid of n arbitrary bins in each dimension would require n^d bins, which grows very rapidly with d
    // to become impractical (In 8 dimensions, with 4 bins along each axis, the produces 65536 bins.) LePage proposed instead to
    // write a seperable approximation p(x) = p_1(x_1) \cdots p_d(x_d) and do a one-dimensional binned approximation in each dimension.
    // This requires only n d bins (In 8 dimensions, with 4 bins along each axis, this produces 32 bins.)

    internal class LePageGrid {

        internal LePageGrid (double[][] grid) {
            this.dimension = grid.Length;
            this.grid = grid;
            binSum = new double[dimension][];
            v0 = 1.0;
            for (int i = 0; i < binSum.Length; i++) {
                int binCount = grid[i].Length - 1;
                binSum[i] = new double[binCount];
                v0 *= binCount;
            }
        }

        public LePageGrid (int dimension, int binCount) {
            this.dimension = dimension;
            this.grid = new double[dimension][];
            this.binSum = new double[dimension][];
            for (int i = 0; i < grid.Length; i++) {
                grid[i] = new double[binCount + 1];
                binSum[i] = new double[binCount];
                for (int j = 0; j < grid[i].Length; j++) {
                    grid[i][j] = ((double) j) / binCount;
                }
            }
            v0 = MoreMath.Pow(binCount, dimension);
        }

        public LePageGrid (IList<Interval> box, int binCount) {
            this.dimension = box.Count;
            this.grid = new double[dimension][];
            this.binSum = new double[dimension][];
            for (int i = 0; i < grid.Length; i++) {
                grid[i] = new double[binCount + 1];
                binSum[i] = new double[binCount];
                for (int j = 0; j < grid[i].Length; j++) {
                    double z = ((double) j) / binCount;
                    grid[i][j] = (1.0 - z) * box[i].LeftEndpoint + z * box[i].RightEndpoint;
                }
            }
            v0 = MoreMath.Pow(binCount, dimension);
        }

        // The dimension of the grid.
        private readonly int dimension;

        // The bin borders for each dimension. This is a D X (N + 1) array.
        private readonly double[][] grid;

        // The values we are trying to smooth for each bin. This is a D X N array.
        private readonly double[][] binSum;

        // The probability density p along each dimension is the probability of the cell (1/N)
        // divided by the width of the cell w. The jacobian is the inverse of the probability
        // density, so the product of N w along each dimension. The total product will thus
        // always contain a factor N^D, which we pre-compute to avoid unnecessary multiplies.
        private readonly double v0;

        int count = 0;

        // Evaluate at x. x is assumed to be in [0,1]^d and is mapped to a bin.
        // Then f(y) / p(y) is evaluated and the result recorded for future refinements
        // before being returned. Note x is changed upon return.

        public double Evaluate (MultiFunctor f, double[] x) {
            return (Evaluate(f, null, x));
        }

        public double Evaluate (MultiFunctor f, CoordinateTransform[] map, double[] x) {

            Debug.Assert(x.Length == dimension);

            // Increase the evaluation count.
            count++;

            // We will need to record the bin number into which each coordinate falls
            // in order to acrue the result to the proper bin.
            int[] binIndexes = new int[dimension];

            // Map incoming x into a grid cell and value based on grid.
            double v = v0;
            for (int i = 0; i < x.Length; i++) {
                Debug.Assert((0.0 <= x[i]) && (x[i] < 1.0));
                double z = (grid[i].Length - 1) * x[i];
                int j = (int) Math.Floor(z);
                z = z - j;
                double w = grid[i][j + 1] - grid[i][j];
                x[i] = (1.0 - z) * grid[i][j] + z * grid[i][j + 1];
                v *= w;
                binIndexes[i] = j;
            }

            // Use the map to further transform that value.
            if (map != null) {
                Debug.Assert(map.Length == dimension);
                for (int i = 0; i < x.Length; i++) {
                    map[i].TransformInPlace(ref x[i], ref v);
                }
            }

            double y = f.Evaluate(x);

            // Record the value in the appropriate bins.
            for (int i = 0; i < binIndexes.Length; i++) {
                int j = binIndexes[i];
                //binSum[i][j] += y * y * v / ((grid[i][j + 1] - grid[i][j]) * (grid[i].Length - 1));
                //binSum[i][j] += Math.Abs(y) * (grid[i][j + 1] - grid[i][j]);
                //binSum[i][j] += Math.Abs(v * y) * (grid[i][j + 1] - grid[i][j]);
                binSum[i][j] += Math.Abs(v * y);
            }

            return (v * y);

        }

        public int EvaluationCount {
            get {
                return (count);
            }
        }

        public int BinCount {
            get {
                return (grid[0].Length - 1);
            }
        }

        public LePageGrid ComputeNewGrid (int m) {

            int d = grid.Length;
            double[][] newGrid = new double[d][];

            for (int i = 0; i < d; i++) {
                newGrid[i] = ComputeNewGrid(grid[i], binSum[i], m);
            }

            return (new LePageGrid(newGrid));
        }

        private static double[] ComputeNewGrid (double[] oldGrid, double[] distribution, int m) {

            int n = oldGrid.Length - 1;
            Debug.Assert(distribution.Length == n);

            // compute normalization factor
            double norm = 0.0;
            for (int i = 0; i < distribution.Length; i++) {
                Debug.Assert(distribution[i] >= 0.0);
                norm += distribution[i];
            }

            // allocate storage for the new grid points
            double[] newGrid = new double[m + 1];

            // the left endpoints of the new and old grids agree
            newGrid[0] = oldGrid[0];

            // j tracks the new grid bin we are currently forming
            int j = 0;

            // since an old bin may contribute an non-integer number of new bins,
            // when we enter an old bin, there is likely some incomplete new bin
            // at its right edge; we need to track this remainer.
            // For the first bin, though, this is zero.
            double remainder = 0.0;

            // iterate over the old bins
            for (int i = 0; i < n; i++) {

                // the remainder should always be a fraction of a bin
                Debug.Assert((0.0 <= remainder) && (remainder < 1.0));

                // compute the number of new bins contributed by the ith old bin
                double newBins = distribution[i] / norm * m;

                // we get more bins because of the remainder
                double netNewBins = newBins + remainder;
                int count = (int) Math.Floor(netNewBins);

                // if we still don't have enough to complete a new bin, add
                // what we have to the remainder and continue to the next old bin
                if (count < 1) {
                    remainder = netNewBins;
                    continue;
                }

                // compute the width of each new bin
                double newBinWidth = (oldGrid[i + 1] - oldGrid[i]) / newBins;
                // this division is safe because we would have continued already if newBins were zero

                // add the new grid points
                for (int k = 1; k <= count; k++) {
                    j++;
                    newGrid[j] = oldGrid[i] + (k - remainder) * newBinWidth;
                    Debug.Assert(newGrid[j - 1] < newGrid[j]);
                    //Debug.Assert((oldGrid[i] <= newGrid[j]) && (newGrid[j] <= oldGrid[i + 1]));
                }

                // pass the remainder to the next iteration
                remainder = netNewBins - count;

            }

            // the right endpoints of the new and old grids agree
            newGrid[m] = oldGrid[n];

            return (newGrid);

        }

    }


    public static partial class MultiFunctionMath {

        private static UncertainValue Integrate_MonteCarlo (MultiFunctor f, CoordinateTransform[] map, IList<Interval> box, EvaluationSettings settings) {

            int d = box.Count;

            // Use a Sobol quasi-random sequence. This give us 1/N accuracy instead of 1/\sqrt{N} accuracy.
            //VectorGenerator g = new RandomVectorGenerator(d, new Random(314159265));
            VectorGenerator g = new SobolVectorGenerator(d);

            // Start with a trivial Lepage grid.
            // We will increase the grid size every few cycles.
            // My tests indicate that trying to increase every cycle or even every other cycle is too often.
            // This makes sense, because we have no reason to believe our new grid will be better until we
            // have substantially more evaluations per grid cell than we did for the previous grid.
            LePageGrid grid = new LePageGrid(box, 1);
            int refineCount = 0;

            // Start with a reasonable number of evaluations per cycle that increases with the dimension.
            int cycleCount = 8 * d;

            //double lastValue = Integrate_MonteCarlo_Cycle(f, map, g, grid, cycleCount);

            // Each cycle consists of three sets of evaluations.
            // At first I did this with just two set and used the difference between the two sets as an error estimate.
            // I found that it was pretty common for that difference to be low just by chance, causing error underestimatation.
            double value1 = Integrate_MonteCarlo_Cycle(f, map, g, grid, cycleCount);
            double value2 = Integrate_MonteCarlo_Cycle(f, map, g, grid, cycleCount);
            double value3 = Integrate_MonteCarlo_Cycle(f, map, g, grid, cycleCount);

            while (f.EvaluationCount < settings.EvaluationBudget) {

                // Take the largest deviation as the error.
                double value = (value1 + value2 + value3) / 3.0;
                double error = Math.Max(Math.Abs(value1 - value3), Math.Max(Math.Abs(value1 - value2), Math.Abs(value2 - value3)));
                Debug.WriteLine("{0} {1} {2}", f.EvaluationCount, value, error);

                // Check for convergence.
                if ((error <= settings.AbsolutePrecision) || (error <= Math.Abs(value) * settings.RelativePrecision)) {
                    return (new UncertainValue(value, error));
                }

                // Do more cycles. In order for new sets to be equal-sized, one of those must be at the current count and the next at twice that.
                double smallValue = Integrate_MonteCarlo_Cycle(f, map, g, grid, cycleCount);
                cycleCount *= 2;
                double bigValue = Integrate_MonteCarlo_Cycle(f, map, g, grid, cycleCount);

                // Combine all the cycles into new ones with twice the number of evaluations each.
                value1 = (value1 + value2) / 2.0;
                value2 = (value3 + smallValue) / 2.0;
                value3 = bigValue;

                //double currentValue = Integrate_MonteCarlo_Cycle(f, map, g, grid, cycleCount);
                //double error = Math.Abs(currentValue - lastValue);
                //double value = (currentValue + lastValue) / 2.0;



                //lastValue = value;

                // Increase the number of evaluations for the next cycle.
                //cycleCount *= 2;

                // Refine the grid for the next cycle.
                refineCount++;
                if (refineCount == 2) {
                    Debug.WriteLine("Replacing grid with {0} bins after {1} evaluations", grid.BinCount, grid.EvaluationCount);
                    grid = grid.ComputeNewGrid(grid.BinCount * 2);
                    refineCount = 0;
                }


            }

            throw new NonconvergenceException();

        }

        // Sample a pre-determined number of points using the given generator and grid and return the average function value.

        private static double Integrate_MonteCarlo_Cycle (MultiFunctor f, CoordinateTransform[] map, VectorGenerator g, LePageGrid grid, int n) {

            double sum = 0.0;

            for (int i = 0; i < n; i++) {
                double[] x = g.NextVector();
                //sum += f.Evaluate(x);
                sum += grid.Evaluate(f, map, x);
            }

            return (sum / n);

        }

    }

}
