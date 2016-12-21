using System;
using System.Collections.Generic;
using System.Diagnostics;


namespace Meta.Numerics.Analysis {


    internal class BulrischStoerStrategy : OdeStrategy {

        public BulrischStoerStrategy (IBulrischStoerEngine engine) : base(engine) {
            Debug.Assert(engine != null);

            this.engine = engine;
            this.ratios = new double[engine.Sizes.Count];
            this.kMin = 2;
            this.kMax = engine.Sizes.Count - 2;
            Debug.Assert(this.kMax > this.kMin);
        }

        private int kMin;

        private int kMax;

        private double[] ratios;

        public override void Step () {

            Debug.Assert(engine.DeltaX != 0.0);

            engine.Clear();

            // Keep track of the total work (in RHS evaluations) to get to the kth column.
            int work = 1;

            bool accept = false;

            // We will use the total work and the computed expansion (or contraction) factor
            // that would target convergence in a given column to compute the most efficient (in terms
            // of distance stepped per evaluation) column in which to target convergence
            // on the next step. We need some variables to keep track of this.
            double bestEfficiency = 0.0;
            double bestFactor = 1.0;
            int bestK = -1;

            for (int k = 0; k <= kMax; k++) {

                engine.AddTrialStep();
                work += engine.Sizes[k];
                ratios[k] = engine.Ratio;

                // We need at least two trial steps in order for our estimate to have an associated
                // uncertainty. For the first step just record the value as the error and move on.
                if (k == 0) continue;

                // Compute the expansion factor required to target convergence in this column on the next step.
                // Compute the corresponding efficiency and determine if this is the best column to target.
                double factor = Math.Pow(1.0 / engine.Ratio, 1.0 / (2 * k + 1));
                double efficiency = factor / work;

                // To give ourselves a margin of safety, we never want to target the final column.
                if (((k + 1) < engine.Sizes.Count) && (efficiency > bestEfficiency)) {
                    bestEfficiency = efficiency;
                    bestFactor = factor;
                    bestK = k;
                }

                // Before the minimum k, don't even check for convergence. It might be spurious.
                // (Any example of this ever happening?)
                if (k < kMin) continue;

                // Check for convergence.
                if (engine.Ratio < 1.0) {

                    accept = true;

                    // We need to figure out whether we should target convergence in a higher column
                    // in the next step. (If we don't we'll never increase the step-size, only
                    // decrease it.) To do that, we need to estimate what the efficiency would be
                    // if we targeted convergence in the next column.

                    // One way to do this, which I did try, is to just do another evaluation
                    // That certainly works, but it is very costly. (It also messes up our
                    // efficiency calculations, since we are really doing more evaluations
                    // that it assumes.) In any case, I've tried it, and the methods below
                    // work better on my test set.

                    // Another way is to use the the estimate that NR quotes from Hairer,
                    // Noersett, and Wanner:
                    //   {\rm err}_{k+1} \approx \left \frac{n_0}{n_{k+1}} \right)^2 {\rm err}_{k}
                    // This estimate is order-of-magnitude as best. It also has a deeper problem:
                    // since it just multiplies the last error by a constant factor, the ratio
                    // of efficiencies it produces will always be the same, so it will always
                    // make the same decision about whether to try for the next column.

                    // Another way is to try to extrapolate from the pattern of previous errors.
                    // Like we do for the function value, we could do extrapolation from a full
                    // polynomial fit (to logs of errors, since they decrease by multiplicative
                    // factors). I tried this, and it works well, but it's a lot of extra
                    // code and a fair amount of increased computational effort.

                    // I also tried extrapolation from just the last three errors and from just
                    // the last two. Surprisingly, extrapolation from just the last two (i.e.
                    // assuming the next step will reduce the error by the same factor as the
                    // last one) resulted in the fewest function evaluations for our test set.

                    // To give ourselves a margin of safety, we never want to target the last
                    // column, so we'll only try this if we are converging at least two
                    // columns from the last.
                    if ((k + 2) < engine.Sizes.Count) {
                        int extrapolatedWork = work + engine.Sizes[k + 1];
                        double extrapolatedRatio = (ratios[k] / ratios[k - 1]) * ratios[k];
                        double extrapolatedFactor = Math.Pow(1.0 / extrapolatedRatio, 1.0 / (2 * k + 3));
                        double extrapolatedEfficiency = extrapolatedFactor / extrapolatedWork;

                        if (extrapolatedEfficiency > bestEfficiency) {
                            bestEfficiency = extrapolatedEfficiency;
                            bestFactor = extrapolatedFactor;
                            bestK = k + 1;
                        }
                    }

                    break;

                }

                // If we are nearly at the maximum, and extrapolation indicates we
                // won't make it, just quit already. This should reduce the number of
                // evaluations, and I've verified that it does on our test set.
                if ((k + 1) == kMax) {
                    double extrapolatedRatio = (ratios[k] / ratios[k - 1]) * ratios[k];
                    if (extrapolatedRatio > 1.0) break;
                }

            }

            engine.CompleteStep(accept);

            // Adjust the step size to target convergence in the optimum column in the next step.
            kMin = Math.Max(bestK - 2, 2);
            kMax = Math.Min(bestK + 2, engine.Sizes.Count - 1);
            if (bestFactor < 0.2) bestFactor = 0.2;
            if (bestFactor > 5.0) bestFactor = 5.0;
            engine.DeltaX *= 0.9375 * bestFactor;

            // Ensure minimum decrease on failed step?

        }

        private readonly IBulrischStoerEngine engine;

    }


    internal interface IBulrischStoerEngine : IOdeEngine {

        IList<int> Sizes { get; }

        void Clear ();

        void AddTrialStep ();

        double Ratio { get; }

    }

}
