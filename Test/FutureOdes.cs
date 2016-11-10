using System;


namespace Test {

    internal abstract class OdeStrategy {


        public abstract bool Step ();

    }

    internal abstract class BulrischStoerState {

        public abstract void ClearExtrapolation ();

        public abstract void AddExtrapolationPoint ();

        public abstract void AcceptExtrapolation ();

    }
    /*
    internal class BulrischStoerOdeStrategy : OdeStrategy {

        private BulrischStoerState state;

        private int kMin, kMax;

        private int[] N;

        public override bool Step () {

            state.ClearExtrapolation();

            int work = 1;

            double bestEfficiency = 0.0;
            double bestFactor = 1.0;
            int bestK = -1;

            for (int k = 0; k < kMax; k++) {

                double[] y1, yp1;
                //TrialStep(N[k], out y1, out yp1);
                //AddExtrapolationPoint(MoreMath.Sqr(1.0 / N[k]), y1, yp1);
                state.AddExtrapolationPoint();

                work += N[k];

                if (k == 0) {
                    double mAbs = 0.0;
                    for (int i = 0; i < y1.Length; i++) {
                        double iAbs = Math.Abs(y1[i]);
                        if (iAbs > mAbs) mAbs = iAbs;
                    }
                    errors[k] = mAbs;
                    continue;
                }

                double yNorm, yError, ypNorm, ypError;
                ComputeExtrapolation(ref y1, out yNorm, out yError, ref yp1, out ypNorm, out ypError);

                double yTol = Settings.ComputePrecision(yNorm);
                double ypTol = Settings.ComputePrecision(ypNorm);

                double yRatio = yTol / yError;
                double ypRatio = ypTol / ypError;
                double ratio = Math.Max(yRatio, ypRatio);
                double factor = Math.Pow(ratio, 1.0 / (2 * k + 1));
                double efficiency = factor / work;
                errors[k] = yError;

                if (((k + 1) < N.Length) && (efficiency > bestEfficiency)) {
                    bestEfficiency = efficiency;
                    bestFactor = factor;
                    bestK = k;
                }

                if (k < kMin) continue;

                if ((yError <= yTol) && (ypError <= ypTol)) {
                    //if (error <= tol) {

                    X += DeltaX;
                    Y = y1;
                    YPrime = yp1;
                    YPrimePrime = Evaluate(X, Y);


                    if ((k + 2) < N.Length) {
                        double extrapolatedError = (errors[k] / errors[k - 1]) * errors[k];
                        int extrapolatedWork = work + N[k + 1];
                        double extrapolatedFactor = Math.Pow(yTol / extrapolatedError, 1.0 / (2 * k + 3));
                        double extrapolatedEfficiency = extrapolatedFactor / extrapolatedWork;

                        if (extrapolatedEfficiency > bestEfficiency) {
                            bestEfficiency = extrapolatedEfficiency;
                            bestFactor = extrapolatedFactor;
                            bestK = k + 1;
                        }
                    }

                    break;



                }

                if (k == (kMax - 1)) {
                    double extrapolatedError = (errors[k] / errors[k - 1]) * errors[k];
                    if (extrapolatedError > yTol) break;
                }

            }

            kMin = Math.Max(bestK - 2, 2);
            kMax = Math.Min(bestK + 2, N.Length - 1);
            if (bestFactor < 0.2) bestFactor = 0.2;
            if (bestFactor > 5.0) bestFactor = 5.0;
            DeltaX *= 0.9375 * bestFactor;

        }

    } 
    */

}
