using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using Meta.Numerics.Matrices;


namespace Meta.Numerics.Analysis {


    public static partial class MultiFunctionMath {

        /// <summary>
        /// Solves a set of coupled, conservative second order ordinary differential equation initial value problems.
        /// </summary>
        /// <param name="rhs">The right hand side function.</param>
        /// <param name="x0">The initial value of the independent variable.</param>
        /// <param name="y0">The initial values of the functions.</param>
        /// <param name="yPrime0">The intial values of the functions' derivatives.</param>
        /// <param name="x1">The final value of the independent variable.</param>
        /// <returns>The solution, including the final value of the functions and their derivatives.</returns>
        public static MultiOdeResult IntegrateConservativeOde (Func<double, IList<double>, IList<double>> rhs, double x0, IList<double> y0, IList<double> yPrime0, double x1) {
            return (IntegrateConservativeOde(rhs, x0, y0, yPrime0, x1, new MultiOdeEvaluationSettings()));
        }

        /// <summary>
        /// Solves a set of coupled, conservative second order ordinary differential equation initial value problems using the given settings.
        /// </summary>
        /// <param name="rhs">The right hand side function.</param>
        /// <param name="x0">The initial value of the independent variable.</param>
        /// <param name="y0">The initial values of the functions.</param>
        /// <param name="yPrime0">The intial values of the functions' derivatives.</param>
        /// <param name="x1">The final value of the independent variable.</param>
        /// <param name="settings">The settings to use when solving the problem.</param>
        /// <returns>The solution, including the final value of the functions and their derivatives.</returns>
        public static MultiOdeResult IntegrateConservativeOde (Func<double, IList<double>, IList<double>> rhs, double x0, IList<double> y0, IList<double> yPrime0, double x1, MultiOdeEvaluationSettings settings) {

            if (rhs == null) throw new ArgumentNullException(nameof(rhs));
            if (y0 == null) throw new ArgumentNullException(nameof(y0));
            if (yPrime0 == null) throw new ArgumentNullException(nameof(yPrime0));
            if (settings == null) throw new ArgumentNullException(nameof(settings));
            if (y0.Count != yPrime0.Count) throw new DimensionMismatchException();

            FunctionMath.SetOdeDefaults(settings);

            MultiStoermerEngine engine = new MultiStoermerEngine(rhs, x0, y0, yPrime0, settings);
            BulrischStoerStrategy strategy = new BulrischStoerStrategy(engine);
            strategy.IntegrateTo(x1);
            return (engine.GetResult());

        }

        /// <summary>
        /// Solves a set of coupled ordinary differential equation initial value problems.
        /// </summary>
        /// <param name="rhs">The right hand side function, which returns the value of the derivative given
        /// the values of the indepdent variable and the function.</param>
        /// <param name="x0">The initial value of the independent variable.</param>
        /// <param name="y0">The initial value of the function.</param>
        /// <param name="x1">The final value of the independent variable.</param>
        /// <returns>The solution, including the final value of the function and its derivative.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="rhs"/> or <paramref name="y0"/> is <see langword="null"/>.</exception>
        /// <exception cref="NonconvergenceException">The ODE could not be integrated to the required precision before exhausting
        /// the maximum allowed number of <paramref name="rhs"/> evaluations.</exception>
        /// <remarks>
        /// <para>The default settings for ODE integration are a relative precision of about 10<sup>-12</sup>, an absolute precision
        /// of about 10<sup>-24</sup>, and a maximum of about 8000 right-hand-side evaluations.</para>
        /// </remarks>
        public static MultiOdeResult IntegrateOde (Func<double, IList<double>, IList<double>> rhs, double x0, IList<double> y0, double x1) {
            return (IntegrateOde(rhs, x0, y0, x1, new MultiOdeEvaluationSettings()));
        }

        /// <summary>
        /// Solves a set of coupled ordinary differential equation initial value problems.
        /// </summary>
        /// <param name="rhs">The right hand side function, which returns the value of the derivative given
        /// the values of the indepdent variable and the function.</param>
        /// <param name="x0">The initial value of the independent variable.</param>
        /// <param name="y0">The initial value of the function.</param>
        /// <param name="x1">The final value of the independent variable.</param>
        /// <param name="settings">The settings to use when solving the problem.</param>
        /// <returns>The solution, including the final value of the function and its derivative.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="rhs"/>, <paramref name="y0"/>, or <paramref name="settings"/>
        /// is <see langword="null"/>.</exception>
        /// <exception cref="NonconvergenceException">The ODE could not be integrated to the required precision before exhausting
        /// the maximum allowed number of <paramref name="rhs"/>evaluations.</exception>
        /// <remarks>
        /// <para>This method integrates a set of coupled ordinary differential equations. The dependent variable y is a vector
        /// with an arbitrary number of components, and the right-hand-side is a vector-valued function that gives the derivative
        /// of each component, and may depend on the values of all the components as well as the independent variable.
        /// The independent variable x still takes only a single real value.</para>
        /// </remarks>
        public static MultiOdeResult IntegrateOde (Func<double, IList<double>, IList<double>> rhs, double x0, IList<double> y0, double x1, MultiOdeEvaluationSettings settings) {

            if (rhs == null) throw new ArgumentNullException(nameof(rhs));
            if (y0 == null) throw new ArgumentNullException(nameof(y0));
            if (settings == null) throw new ArgumentNullException(nameof(settings));

            FunctionMath.SetOdeDefaults(settings);

            MultiBulrischStoerEngine engine = new MultiBulrischStoerEngine(rhs, x0, y0, settings);
            BulrischStoerStrategy strategy = new BulrischStoerStrategy(engine);
            strategy.IntegrateTo(x1);
            return (engine.GetResult());

        }

    }

    internal abstract class MultiOdeEngine : OdeEngine {

        public MultiOdeEngine (Func<double, IList<double>, IList<double>> rhs, double x, MultiOdeEvaluationSettings settings) : base(x) {
            Debug.Assert(rhs != null);
            Debug.Assert(settings != null);
            this.rhs = rhs;
            this.settings = settings;
        }

        private readonly Func<double, IList<double>, IList<double>> rhs;

        protected readonly MultiOdeEvaluationSettings settings;

        public IList<double> Evaluate (double x, IList<double> y) {
            if (count >= settings.EvaluationBudget) throw new NonconvergenceException();
            count++;
            ReadOnlyCollection<double> yWrapper = new ReadOnlyCollection<double>(y);
            IList<double> z = rhs(x, yWrapper);
            if (z == null) throw new InvalidOperationException();
            return (z);
        }

        private int count = 0;

        public int EvaluationCount {
            get {
                return (count);
            }
        }

        public abstract MultiOdeResult GetResult ();

        protected override void AcceptStep () {
            base.AcceptStep();
            if (settings.Listener != null) settings.Listener(GetResult());
        }

    }

    internal class MultiBulrischStoerEngine : MultiOdeEngine, IBulrischStoerEngine {

        public MultiBulrischStoerEngine (Func<double, IList<double>, IList<double>> rhs, double x, IList<double> y, MultiOdeEvaluationSettings settings) : base(rhs, x, settings) {

            Debug.Assert(rhs != null);
            Debug.Assert(y != null);
            Debug.Assert(settings != null);

            Y = new double[y.Count];
            y.CopyTo(Y, 0);
            YPrime = new double[y.Count];
            Evaluate(X, Y).CopyTo(YPrime, 0);
            ComputeInitialStep();

            extrapolators = new NevilleExtrapolator[y.Count];
            for (int i = 0; i < extrapolators.Length; i++) extrapolators[i] = new NevilleExtrapolator(N.Length);
        }

        private void ComputeInitialStep () {
            double yNorm = 0.0;
            double yPrimeNorm = 0.0;
            for (int i = 0; i < YPrime.Length; i++) {
                double yAbs = Math.Abs(Y[i]);
                if (yAbs > yNorm) yNorm = yAbs;
                double yPrimeAbs = Math.Abs(YPrime[i]);
                if (yPrimeAbs > yPrimeNorm) yPrimeNorm = yPrimeAbs;
            }
            DeltaX = yNorm / yPrimeNorm;
            if (Double.IsNaN(DeltaX) || DeltaX == 0.0) DeltaX = 0.5;
        }

        double[] Y;

        double[] YPrime;

        private static readonly int[] N = new int[] { 2, 4, 6, 8, 10, 12, 14, 16, 18 };

        public virtual IList<int> Sizes { get { return (N); } }

        private double[] TrialStep (int n) {

            double h = DeltaX / n;

            double[] nY0 = new double[Y.Length];
            double[] nY1 = new double[Y.Length];
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

        public virtual void Clear () {
            for (int i = 0; i < extrapolators.Length; i++) extrapolators[i].Clear();
        }

        public virtual void AddTrialStep () {

            int k = extrapolators[0].Count;
            double[] y1 = TrialStep(N[k]);
            double controlValue = MoreMath.Sqr(1.0 / N[k]);
            for (int i = 0; i < extrapolators.Length; i++) extrapolators[i].Add(controlValue, y1[i]);

            double norm = 0.0;
            double error = 0.0;
            for (int i = 0; i < extrapolators.Length; i++) {
                double aValue = Math.Abs(extrapolators[i].Value);
                if (aValue > norm) norm = aValue;
                double uncertainty = extrapolators[i].Error;
                if (uncertainty > error) error = uncertainty;
            }
            double tol = settings.ComputePrecision(norm);
            this.Ratio = error / tol;

        }

        public virtual double Ratio {
            get; private set;
        }

        protected override void AcceptStep () {
            for (int i = 0; i < extrapolators.Length; i++) {
                Y[i] = extrapolators[i].Value;
            }
            Evaluate(X + DeltaX, Y).CopyTo(YPrime, 0);
            base.AcceptStep();
        }

        public override MultiOdeResult GetResult () {
            return (new MultiOdeResult(
                this.X, this.Y, this.YPrime, this.EvaluationCount, this.settings
            ));
        }

    }

    internal class MultiStoermerEngine : MultiOdeEngine, IBulrischStoerEngine {

        public MultiStoermerEngine (Func<double, IList<double>, IList<double>> rhs, double x, IList<double> y, IList<double> yPrime, MultiOdeEvaluationSettings settings) : base(rhs, x, settings) {

            Debug.Assert(rhs != null);
            Debug.Assert(y != null);
            Debug.Assert(yPrime != null);
            Debug.Assert(settings != null);
            Debug.Assert(y.Count == yPrime.Count);

            dimension = y.Count;

            Y = new double[dimension];
            y.CopyTo(Y, 0);
            YPrime = new double[dimension];
            yPrime.CopyTo(YPrime, 0);

            YPrimePrime = new double[dimension];
            Evaluate(x, y).CopyTo(YPrimePrime, 0);

            ComputeInitialStep();

            yExtrapolators = new NevilleExtrapolator[y.Count];
            yPrimeExtrapolators = new NevilleExtrapolator[yPrime.Count];
            for (int i = 0; i < y.Count; i++) {
                yExtrapolators[i] = new NevilleExtrapolator(N.Length);
                yPrimeExtrapolators[i] = new NevilleExtrapolator(N.Length);
            }
        }


        private void ComputeInitialStep () {
            this.DeltaX = 0.5;
        }

        private readonly int dimension;

        private double[] Y;

        private double[] YPrime;

        private double[] YPrimePrime;

        public override MultiOdeResult GetResult () {
            return (new MultiOdeResult(
                this.X, this.Y, this.YPrime, this.EvaluationCount, this.settings
            ));
        }

        public double Ratio {
            get; private set;
        }

        private static readonly int[] N = new int[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

        public IList<int> Sizes {
            get {
                return (N);
            }
        }

        protected override void AcceptStep () {
            for (int i = 0; i < dimension; i++) {
                Y[i] = yExtrapolators[i].Value;
                YPrime[i] = yPrimeExtrapolators[i].Value;
            }
            Evaluate(X + DeltaX, Y).CopyTo(YPrimePrime, 0);
            base.AcceptStep();
        }

        public void AddTrialStep () {

            int k = yExtrapolators[0].Count;
            double parameter = MoreMath.Sqr(1.0 / N[k]);

            double[] y1, yp1;
            TrialStep(N[k], out y1, out yp1);
            for (int i = 0; i < dimension; i++) {
                yExtrapolators[i].Add(parameter, y1[i]);
                yPrimeExtrapolators[i].Add(parameter, yp1[i]);
            }

            double yNorm, yError, ypNorm, ypError;
            ComputeNormAndError(out yNorm, out yError, out ypNorm, out ypError);

            double yTolerance = settings.ComputePrecision(yNorm);
            double yRatio = yError / yTolerance;

            double ypTolerance = settings.ComputePrecision(ypNorm);
            double ypRatio = ypError / ypTolerance;

            this.Ratio = Math.Max(yRatio, ypRatio);

        }

        private NevilleExtrapolator[] yExtrapolators;
        private NevilleExtrapolator[] yPrimeExtrapolators;

        private void ComputeNormAndError (out double yNorm, out double yError, out double ypNorm, out double ypError) {
            yNorm = 0.0;
            yError = 0.0;
            ypNorm = 0.0;
            ypError = 0.0;
            for (int i = 0; i < dimension; i++) {
                yNorm += MoreMath.Sqr(yExtrapolators[i].Value);
                yError += MoreMath.Sqr(yExtrapolators[i].Error);
                ypNorm += MoreMath.Sqr(yPrimeExtrapolators[i].Value);
                ypError += MoreMath.Sqr(yPrimeExtrapolators[i].Error);
            }
            yNorm = Math.Sqrt(yNorm);
            yError = Math.Sqrt(yError);
            ypNorm = Math.Sqrt(ypNorm);
            ypError = Math.Sqrt(ypError);
        }


        private void TrialStep (int n, out double[] Y1, out double[] YP1) {

            Debug.Assert(n >= 1);

            double h = DeltaX / n;

            int d = Y.Length;

            IList<double> F;
            Y1 = new double[d];
            double[] D1 = new double[d];

            for (int i = 0; i < d; i++) {
                Y1[i] = Y[i];
                D1[i] = h * (YPrime[i] + h * YPrimePrime[i] / 2.0);
            }

            for (int k = 1; k < n; k++) {
                for (int i = 0; i < d; i++) Y1[i] += D1[i];
                F = Evaluate(X + k * h, Y1);
                for (int i = 0; i < d; i++) D1[i] += h * h * F[i];
            }

            for (int i = 0; i < d; i++) Y1[i] += D1[i];
            F = Evaluate(X + DeltaX, Y1);
            for (int i = 0; i < d; i++) D1[i] = D1[i] / h + h * F[i] / 2.0;

            YP1 = D1;

        }

        public void Clear () {
            for (int i = 0; i < yExtrapolators.Length; i++) {
                yExtrapolators[i].Clear();
                yPrimeExtrapolators[i].Clear();
            }
        }
    }

}
