using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;


namespace Meta.Numerics.Analysis {

#if FUTURE
    public class OdeEvaulationSettings : EvaluationSettings {

        public Action<OdeResult> Listener { get; set; }

    }

    public class MultiOdeEvaluationSettings : EvaluationSettings {

        public Action<MultiOdeResult> Listener { get; set; }

    }


    public static class NewMethods {

        public static OdeResult IntegrateOde(Func<double, double, double> rhs, double x0, double y0, double x1) {

            if (rhs == null) throw new ArgumentNullException(nameof(rhs));

            SingleBulrischStoerEngine engine = new SingleBulrischStoerEngine(rhs, x0, y0);
            BulrischStoerStrategy<OdeResult> strategy = new BulrischStoerStrategy<OdeResult>(engine, null);
            strategy.IntegrateTo(x1);
            return (engine.GetResult());

        }

    }

    internal abstract class OdeStrategy<T> {

        protected OdeStrategy (OdeEngine<T> engine, Action<T> listener) {
            this.engine = engine;
            this.listener = listener;
        }

        private readonly OdeEngine<T> engine;

        private readonly Action<T> listener;

        public void IntegrateTo (double x1) {

            if (listener != null) {
                T result = engine.GetResult();
                listener(result);
            }

        }

    }

    internal class BulrischStoerStrategy<T> : OdeStrategy<T> {

        public BulrischStoerStrategy (OdeEngine<T> engine, Action<T> listener) : base(engine, listener) {

        }

        private int kMin;

        private int kMax;

        private double[] errors;

        private EvaluationSettings settings;

        public void Step () {

            engine.Clear();

            int work = 1;

            double bestEfficiency = 0.0;
            double bestFactor = 1.0;
            int bestK = -1;

            for (int k = 0; k < kMax; k++) {

                engine.AddTrialStep();
                work += engine.Sizes[k];

                if (k == 0) {
                    errors[k] = engine.Norm;
                    continue;
                }

                double norm = engine.Norm;
                double error = engine.Error;
                double tol = settings.ComputePrecision(norm);

                double factor = Math.Pow(tol / error, 1.0 / (2 * k + 1));
                double efficiency = factor / work;
                errors[k] = error;

                if (((k + 1) < engine.Sizes.Count) && (efficiency > bestEfficiency)) {
                    bestEfficiency = efficiency;
                    bestFactor = factor;
                    bestK = k;
                }

                if (k < kMin) continue;

                if (error <= tol) {

                    engine.Accept();

                    if ((k + 2) < engine.Sizes.Count) {
                        double extrapolatedError = (errors[k] / errors[k - 1]) * errors[k];
                        int extrapolatedWork = work + engine.Sizes[k + 1];
                        double extrapolatedFactor = Math.Pow(tol / extrapolatedError, 1.0 / (2 * k + 3));
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
                    if (extrapolatedError > tol) break;
                }

            }

            kMin = Math.Max(bestK - 2, 2);
            kMax = Math.Min(bestK + 2, engine.Sizes.Count - 1);
            if (bestFactor < 0.2) bestFactor = 0.2;
            if (bestFactor > 5.0) bestFactor = 5.0;
            engine.DeltaX *= 0.9375 * bestFactor;


        }

        private BulrischStoerEngine<T> engine;

    }


    internal class OdeEngine {

        public double X;

        public double DeltaX;


    }

    internal abstract class OdeEngine<T> : OdeEngine {

        public abstract T GetResult ();

    }


    internal interface IBulrischStoerEngine {

        IList<int> Sizes { get; }

        void Clear ();

        void AddTrialStep ();

        void Accept ();

    }

    internal abstract class BulrischStoerEngine<T> : OdeEngine<T> {

        public abstract IList<int> Sizes { get; }

        public abstract void Clear ();

        public abstract void Accept ();

        public abstract int AddTrialStep ();

        public double Error { get; set; }

        public double Norm { get; set; }

    }


    internal class OdeFunctor {

        public OdeFunctor (Func<double, double, double> rhs) {
            this.rhs = rhs;
        }

        private Func<double, double, double> rhs;

        public double Evaluate (double x, double y) {
            count++;
            return (rhs(x, y));
        }

        private int count;

        public int EvaluationCount {
            get {
                return (count);
            }
        }

    }

    internal class SingleStoermerEngine : OdeEngine<OdeResult>, IBulrischStoerEngine {

        private Func<double, double, double> rhs;
        private double Evaluate (double x, double y) {
            count++;
            return (rhs(x, y));
        }

        private int count;

        private double Y;

        private double YPrime;

        private double YPrimePrime;

        private static int[] N = new int[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

        public IList<int> Sizes {  get { return (N); } }

        private void TrialStep (int n, out double Y1, out double YP1) {

            Debug.Assert(n >= 1);

            double h = DeltaX / n;

            Y1 = Y;
            double D1 = h * (YPrime + h * YPrimePrime / 2.0);
            for (int k = 1; k < n; k++) {
                Y1 += D1;
                D1 += h * h * Evaluate(X + k * h, Y1);
            }
            Y1 += D1;

            YP1 = D1 / h + h * Evaluate(X + DeltaX, Y1) / 2.0;

        }

        private NevilleExtrapolator yExtrapolator = new NevilleExtrapolator(N.Length);
        private NevilleExtrapolator ypExtrapolator = new NevilleExtrapolator(N.Length);

        public void Clear () {
            yExtrapolator.Clear();
            ypExtrapolator.Clear();
        }

        public void AddTrialStep () {

            int k = yExtrapolator.Count;

            double y1, yp1;
            TrialStep(N[k], out y1, out yp1);

            yExtrapolator.Add(MoreMath.Sqr(1.0 / N[k]), y1);
            ypExtrapolator.Add(MoreMath.Sqr(1.0 / N[k]), yp1);

            Norm = Math.Max(
                Math.Abs(yExtrapolator.Estimate.Value),
                Math.Abs(ypExtrapolator.Estimate.Value)
            );
            Error = Math.Max(yExtrapolator.Estimate.Uncertainty, ypExtrapolator.Estimate.Uncertainty);

        }

        public double Norm { get; private set; }

        public double Error { get; private set; }

        public void Accept () {

            X += DeltaX;
            Y = yExtrapolator.Estimate.Value;
            YPrime = ypExtrapolator.Estimate.Value;
            YPrimePrime = Evaluate(X, Y);

        }
        
        public override OdeResult GetResult () {

            return (new OdeResult(
                this.X,
                this.Y,
                this.YPrime,
                this.count,
                null
            ));

        }

    }

    internal class SingleBulrischStoerEngine : BulrischStoerEngine<OdeResult> {

        public SingleBulrischStoerEngine(Func<double, double, double> rhs, double x, double y) {
            this.rhs = rhs;
            this.X = x;
            this.Y = y;
            this.YPrime = rhs(x, y);
        }

        private Func<double, double, double> rhs;

        private double Evaluate (double x, double y) {
            count++;
            return (rhs(x, y));
        }

        private int count = 0;

        private double Y;

        private double YPrime;

        private static readonly int[] N = new int[] { 2, 4, 6, 8, 10, 12, 14, 16, 18, 20 };

        public override IList<int> Sizes {  get { return (N); } }

        private double TrialStep (int n) {

            double h = DeltaX / n;
            double two_h = 2.0 * h;

            double nY0 = Y;
            double nY1 = Y + h * YPrime;

            for (int k = 1; k < n; k++) {
                double nY2 = nY0 + two_h * Evaluate(X + k * h, nY1);
                nY0 = nY1;
                nY1 = nY2;
            }

            nY1 = (nY0 + nY1 + h * Evaluate(X + DeltaX, nY1)) / 2.0;
            return (nY1);

        }

        private NevilleExtrapolator extrapolator;
        private double[] errors = new double[N.Length];


        public override void Clear () {
            extrapolator.Clear();
        }

        public override void Accept () {

            X += DeltaX;
            Y = extrapolator.Estimate.Value;
            YPrime = Evaluate(X, Y);

        }

        public override int AddTrialStep () {

            int k = extrapolator.Count;

            double y = TrialStep(N[k]);
            extrapolator.Add(MoreMath.Sqr(1.0 / N[k]), y);

            Norm = Math.Abs(extrapolator.Estimate.Value);
            Error = extrapolator.Estimate.Uncertainty;

            return (N[k]);
        }

        public int Work {
            get {
                throw new NotImplementedException();
            }
        }

        public override OdeResult GetResult () {

            return (new OdeResult(
                this.X,
                this.Y,
                this.YPrime,
                this.count,
                null
            ));

        }

    }


    internal class MultiBulrischStoerEngine : BulrischStoerEngine<MultiOdeResult> {

        private Func<double, IList<double>, IList<double>> rhs;

        private IList<double> Evaluate (double x, IList<double> y) {
            count++;
            return (rhs(x, new ReadOnlyCollection<double>(y)));
        }

        private int count;

        IList<double> Y;

        IList<double> YPrime;

        private static readonly int[] N = new int[] { 2, 4, 6, 8, 10, 12, 14, 16, 18 };

        public override IList<int> Sizes {  get { return (N); } }

        private double[] TrialStep (int n) {

            double h = DeltaX / n;

            double[] nY0 = new double[Y.Count];
            double[] nY1 = new double[Y.Count];
            for (int i = 0; i < nY1.Length; i++) { nY0[i] = Y[i]; nY1[i] = Y[i] + h * YPrime[i]; }

            IList<double> F;
            for (int k = 1; k < n; k++) {
                F = Evaluate(X + k * h, nY1);
                for (int i = 0; i < nY1.Length; i++) { double t = nY1[i]; nY1[i] = nY0[i] + 2.0 * h * F[i]; nY0[i] = t; }
            }

            F = Evaluate(X + DeltaX, nY1);
            for (int i = 0; i < nY1.Length; i++) nY1[i] = (nY0[i] + nY1[i] + h * F[i]) / 2.0;

            return (nY1);

        }

        private NevilleExtrapolator[] extrapolators;

        public override void Clear () {
            for (int i = 0; i < extrapolators.Length; i++) extrapolators[i].Clear();
        }

        public override int AddTrialStep () {

            int k = extrapolators[0].Count;
            double[] y1 = TrialStep(N[k]);
            double controlValue = MoreMath.Sqr(1.0 / N[k]);
            for (int i = 0; i < extrapolators.Length; i++) extrapolators[i].Add(controlValue, y1[i]);

            double[] value = new double[extrapolators.Length];
            double norm = 0.0;
            double error = 0.0;
            for (int i = 0; i < extrapolators.Length; i++) {
                UncertainValue estimate = extrapolators[i].Estimate;
                value[i] = estimate.Value;
                double aValue = Math.Abs(estimate.Value);
                if (aValue > norm) norm = aValue;
                if (estimate.Uncertainty > error) error = estimate.Uncertainty;
            }
            this.Norm = norm;
            this.Error = error;

            return (N[k]);
        }

        public override void Accept () {
            throw new NotImplementedException();
        }

        public override MultiOdeResult GetResult () {

            return (new MultiOdeResult(
                this.X, null, null, count, null
            ));
        }

    }
#endif
}
