using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Text;

using Meta.Numerics;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Analysis {

    public class OdeResult<T> : EvaluationResult {

        public OdeResult (int count, EvaluationSettings settings, double x, T y, T yPrime) : base(count, settings) {
            this.X = x;
            this.Y = y;
            this.YPrime = yPrime;
        }

        public double X { get; private set; }

        public T Y { get; private set; }

        public T YPrime { get; private set; }

    }

    public static partial class FunctionMath {

        private static EvaluationSettings defaultOdeSettings = new EvaluationSettings() { RelativePrecision = 1.0E-12, AbsolutePrecision = 1.0E-24, EvaluationBudget = 5000 };

        public static double SolveConservativeOde (Func<double, double, double> rhs, double x0, double y0, double yp0, double x1) {
            return (SolveConservativeOde(rhs, x0, y0, yp0, x1, defaultOdeSettings));
        }


        public static double SolveConservativeOde (Func<double, double, double> rhs, double x0, double y0, double yp0, double x1, EvaluationSettings settings) {
            if (rhs == null) throw new ArgumentNullException("rhs");
            BaseOdeStepper<double> stepper = new StoermerStepper(rhs, x0, y0, yp0, settings);
            stepper.Integrate(x1);
            return (stepper.Y);

        }

        public static double SolveOde (Func<double, double, double> rhs, double x0, double y0, double x1) {
            return (SolveOde(rhs, x0, y0, x1, defaultOdeSettings));
        }

        public static double SolveOde (Func<double, double, double> rhs, double x0, double y0, double x1, EvaluationSettings settings) {
            if (rhs == null) throw new ArgumentNullException("rhs");

            MultiOdeStepper stepper = new BulrischStoerStepper((double x, IList<double> y) => new double[] { rhs(x, y[0]) }, x0, new double[] { y0 }, settings);
            //MultiOdeStepper stepper = new RungeKutta54Stepper((double x, IList<double> y) => new double[] { rhs(x, y[0]) }, range.LeftEndpoint, new double[] { start }, settings);
            stepper.Integrate(x1);
            return (stepper.Y[0]);

        }

    }

    public static partial class MultiFunctionMath {

        public static ColumnVector SolveConservativeOde (Func<double, IList<double>, IList<double>> rhs, double x0, IList<double> y0, IList<double> yp0, double x1) {
            EvaluationSettings settings = new EvaluationSettings() { RelativePrecision = 1.0E-12, AbsolutePrecision = 1.0E-24, EvaluationBudget = 4096 };
            return (SolveConservativeOde(rhs, x0, y0, yp0, x1, settings));
        }

        public static ColumnVector SolveConservativeOde (Func<double, IList<double>, IList<double>> rhs, double x0, IList<double> y0, IList<double> yp0, double x1, EvaluationSettings settings) {
            if (rhs == null) throw new ArgumentNullException("rhs");
            if (y0 == null) throw new ArgumentNullException("y0");

            MultiOdeStepper stepper = new MultiStoermerStepper(rhs, x0, y0, yp0, settings);
            stepper.Integrate(x1);
            return (new ColumnVector(stepper.Y));
        }

        public static OdeResult<IList<double>> SolveOde (Func<double, IList<double>, IList<double>> rhs, double x0, IList<double> y0, double x1) {
            EvaluationSettings settings = new EvaluationSettings() { RelativePrecision = 1.0E-12, AbsolutePrecision = 1.0E-24, EvaluationBudget = 4096 };
            return (SolveOde(rhs, x0, y0, x1, settings));
        }

        public static OdeResult<IList<double>> SolveOde (Func<double, IList<double>, IList<double>> rhs, double x0, IList<double> y0, double x1, EvaluationSettings settings) {
            if (rhs == null) throw new ArgumentNullException("rhs");
            if (y0 == null) throw new ArgumentNullException("y0");
            if (settings == null) throw new ArgumentException("settings");

            MultiOdeStepper stepper = new BulrischStoerStepper(rhs, x0, y0, settings);
            stepper.Integrate(x1);

            return (stepper.Result);
        }

    }

    // The base class for all ODE steppers.

    internal abstract class BaseOdeStepper<T> {

        public BaseOdeStepper(Func<double, T, T> rhs, double x0, T y0, EvaluationSettings settings) {
            if (rhs == null) throw new ArgumentNullException("rhs");
            if (settings == null) throw new ArgumentNullException("settings");
            this.rhs = rhs;
            this.X = x0;
            this.Y = y0;
            this.Settings = settings;
        }

        // The right-hand-side. This is kept private so it can only be called through Evaluate.

        private Func<double, T, T> rhs;

        /// <summary>
        /// The current value of the independent variable.
        /// </summary>
        public double X { get; protected set; }

        /// <summary>
        /// The current value of the dependent variable.
        /// </summary>
        public T Y { get; protected set; }

        /// <summary>
        /// The current value of the first derivative.
        /// </summary>
        public T YPrime { get; protected set; }

        protected virtual T Evaluate (double x, T y) {
            if (this.EvaluationCount >= Settings.EvaluationBudget) throw new NonconvergenceException();
            EvaluationCount++;
            return (rhs(x, y));
        }

        /// <summary>
        /// The number of right-hand-side evaluations.
        /// </summary>
        public int EvaluationCount { get; private set; }

        /// <summary>
        /// The current step size.
        /// </summary>
        public double DeltaX { get; protected set; }

        /// <summary>
        /// The number of RHS evaluations.
        /// </summary>
        public EvaluationSettings Settings { get; private set; }

        // This is the key method to be implemented by any stepper.
        // After it is called, if the step is sucessful, it should advance X (presumably by DeltaX but smaller is okay),
        // and update the value of Y to correspond to the new X. If the step is unsuccessful, it should leave X and Y
        // unchanged. If the stepper is adaptive, the method may also change DeltaX in preperation for the next step.

        public abstract void Step ();

        public virtual void Integrate (double X1) {

            double X0 = X;

            // reverse direction, if necessary
            if (Math.Sign(DeltaX) != Math.Sign(X1 - X0)) DeltaX = -DeltaX;

            // we can't just check (X < X1) because sometimes we integrate the other way
            // so instead check that "we are on the same side of X1 as X0"
            while (Math.Sign(X - X1) == Math.Sign(X0 - X1)) {

                // if we would overshoot in the next step, reduce it
                if (Math.Sign(X + DeltaX - X1) != Math.Sign(X0 - X1)) DeltaX = X1 - X;

                double X_old = X;

                Step();

                if (X != X_old) {
                    SuccessfulStepCount++;
                } else {
                    FailedStepCount++;
                }

                if (Settings.UpdateHandler != null) {
                    Settings.Update(this.Result);
                }

            }

        }

        /// <summary>
        /// The number of successful steps.
        /// </summary>
        public int SuccessfulStepCount { get; private set; }

        /// <summary>
        /// The number of failed steps.
        /// </summary>
        public int FailedStepCount { get; private set; }

        /// <summary>
        /// The total number of completed steps, both failed and successful.
        /// </summary>
        public int AttemptedStepCount {
            get {
                return (FailedStepCount + SuccessfulStepCount);
            }
        }

        public virtual OdeResult<T> Result {
            get {
                return (new OdeResult<T>(EvaluationCount, Settings, X, Y, YPrime));
            }
        }
    }

    // This is a stepper for 2nd order ODEs with where the 1st derivative does not appear, i.e.
    //   \frac{d^2 y}{dx^2} = f(x, y)
    // The evolution of such systems ensures a conserved quantity (e.g. energy).

    internal class StoermerStepper : BaseOdeStepper<double> {

        public StoermerStepper (Func<double, double, double> rhs, double x0, double y0, double yp0, EvaluationSettings settings) : base(rhs, x0, y0, settings) {
            YPrime = yp0;
            YPrimePrime = rhs(x0, y0);
            DeltaX = 1.0;
        }

        public double YPrimePrime { get; private set; }

        private NevilleExtrapolator yExtrapolator = new NevilleExtrapolator(N.Length);
        private NevilleExtrapolator ypExtrapolator = new NevilleExtrapolator(N.Length);

        public override void Step () {

            yExtrapolator.Clear();
            ypExtrapolator.Clear();

            int work = 1;

            double bestEfficiency = 0.0;
            double bestFactor = 1.0;
            int bestK = -1;

            for (int k = 0; k < N.Length; k++) {

                double y1, yp1;
                TrialStep(N[k], out y1, out yp1);
                yExtrapolator.Add(MoreMath.Sqr(1.0 / N[k]), y1);
                ypExtrapolator.Add(MoreMath.Sqr(1.0 / N[k]), yp1);

                work += N[k];

                if (k < 1) continue;

                UncertainValue yEstimate = yExtrapolator.Estimate;

                double error = yEstimate.Uncertainty;

                double tol = Settings.ComputePrecision(yEstimate.Value);

                double factor = Math.Pow(tol / error, 1.0 / (2 * k + 1));

                double efficiency = factor / work;

                if (((k+1) < N.Length) && (efficiency > bestEfficiency)) {
                    bestEfficiency = efficiency;
                    bestFactor = factor;
                    bestK = k;
                }

                
                if (error <= tol) {

                    X += DeltaX;
                    Y = yEstimate.Value;
                    YPrime = ypExtrapolator.Estimate.Value;
                    YPrimePrime = Evaluate(X, Y);

                    if ((k + 2) < N.Length) {
                        double extrapolatedError = MoreMath.Sqr(1.0 * N[0] / N[k + 1]) * error;
                        int extrapolatedWork = work + N[k + 1];
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
                
            }

            if (bestFactor < 0.2) bestFactor = 0.2;
            if (bestFactor > 5.0) bestFactor = 5.0;
            DeltaX *= 0.9375 * bestFactor;

        }

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

        private static int[] N = new int[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    }


    internal class MultiStoermerStepper : MultiOdeStepper {

        public MultiStoermerStepper (Func<double, IList<double>, IList<double>> rhs, double x0, IList<double> y0, IList<double> yp0, EvaluationSettings settings) : base(rhs, x0, y0, settings) {
            YPrime = yp0;
            YPrimePrime = rhs(x0, y0);
            DeltaX = 1.0;
            yExtrapolators = new NevilleExtrapolator[Dimension];
            ypExtrapolators = new NevilleExtrapolator[Dimension];
            for (int i = 0; i < Dimension; i++) {
                yExtrapolators[i] = new NevilleExtrapolator(N.Length);
                ypExtrapolators[i] = new NevilleExtrapolator(N.Length);
            }
        }

        public IList<double> YPrimePrime { get; private set; }

        private NevilleExtrapolator[] yExtrapolators;
        private NevilleExtrapolator[] ypExtrapolators;

        private void ClearExtrapolation () {
            for (int i = 0; i < Dimension; i++) {
                yExtrapolators[i].Clear();
                ypExtrapolators[i].Clear();
            }
        }

        private void AddExtrapolationPoint (double parameter, double[] y, double[] yp) {
            for (int i = 0; i < Dimension; i++) {
                yExtrapolators[i].Add(parameter, y[i]);
                ypExtrapolators[i].Add(parameter, yp[i]);
            }
        }

        private void ComputeExtrapolation (ref double[] y, out double yNorm, out double yError, ref double[] yp, out double ypNorm, out double ypError) {
            yNorm = 0.0;
            yError = 0.0;
            ypNorm = 0.0;
            ypError = 0.0;
            for (int i = 0; i < Dimension; i++) {
                UncertainValue yEstimate = yExtrapolators[i].Estimate;
                y[i] = yEstimate.Value;
                yNorm += MoreMath.Sqr(yEstimate.Value);
                yError += MoreMath.Sqr(yEstimate.Uncertainty);
                UncertainValue ypEstimate = ypExtrapolators[i].Estimate;
                yp[i] = ypEstimate.Value;
                ypNorm += MoreMath.Sqr(ypEstimate.Value);
                ypError += MoreMath.Sqr(ypEstimate.Uncertainty);
            }
            yNorm = Math.Sqrt(yNorm);
            yError += Math.Sqrt(yError);
            ypNorm = Math.Sqrt(ypNorm);
            ypError += Math.Sqrt(ypError);
        }

        public override void Step () {

            ClearExtrapolation();

            int work = 1;

            double bestEfficiency = 0.0;
            double bestFactor = 1.0;
            int bestK = -1;

            for (int k = 0; k < N.Length; k++) {

                double[] y1, yp1;
                TrialStep(N[k], out y1, out yp1);
                AddExtrapolationPoint(MoreMath.Sqr(1.0 / N[k]), y1, yp1);

                work += N[k];

                if (k < 1) continue;

                double yNorm, yError, ypNorm, ypError;
                ComputeExtrapolation(ref y1, out yNorm, out yError, ref yp1, out ypNorm, out ypError);

                double yTol = Settings.ComputePrecision(yNorm);
                double ypTol = Settings.ComputePrecision(ypNorm);

                double yRatio = yTol / yError;
                double ypRatio = ypTol / ypError;
                double ratio = Math.Max(yRatio, ypRatio);

                double factor = Math.Pow(ratio, 1.0 / (2 * k + 1));

                double efficiency = factor / work;

                if (((k + 1) < N.Length) && (efficiency > bestEfficiency)) {
                    bestEfficiency = efficiency;
                    bestFactor = factor;
                    bestK = k;
                }

                if ((yError <= yTol) && (ypError <= ypTol)) {
                //if (error <= tol) {

                    X += DeltaX;
                    Y = y1;
                    YPrime = yp1;
                    YPrimePrime = Evaluate(X, Y);

                    if ((k + 2) < N.Length) {
                        double yExtrapolatedError = MoreMath.Sqr(1.0 * N[0] / N[k + 1]) * yError;
                        double ypExtrapolatedError = MoreMath.Sqr(1.0 * N[0] / N[k + 1]) * ypError;
                        double yExtrapolatedRatio = yTol / yExtrapolatedError;
                        double ypExtrapolatedRatio = ypTol / ypExtrapolatedError;
                        double extrapolatedRatio = Math.Max(yExtrapolatedRatio, ypExtrapolatedRatio);
                        double extrapolatedFactor = Math.Pow(extrapolatedRatio, 1.0 / (2 * k + 3)); 
                        int extrapolatedWork = work + N[k + 1];
                        double extrapolatedEfficiency = extrapolatedFactor / extrapolatedWork;

                        if (extrapolatedEfficiency > bestEfficiency) {
                            bestEfficiency = extrapolatedEfficiency;
                            bestFactor = extrapolatedFactor;
                            bestK = k + 1;
                        }
                    }

                    break;

                }

            }

            if (bestFactor < 0.2) bestFactor = 0.2;
            if (bestFactor > 5.0) bestFactor = 5.0;
            DeltaX *= 0.9375 * bestFactor;

        }

        private void TrialStep (int n, out double[] Y1, out double[] YP1) {

            Debug.Assert(n >= 1);

            double h = DeltaX / n;

            IList<double> F;
            Y1 = new double[Dimension];
            double[] D1 = new double[Dimension];

            for (int i = 0; i < Dimension; i++) {
                Y1[i] = Y[i];
                D1[i] = h * (YPrime[i] + h * YPrimePrime[i] / 2.0);
            }

            for (int k = 1; k < n; k++) {
                for (int i = 0; i < Dimension; i++) Y1[i] += D1[i];
                F = Evaluate(X + k * h, Y1);
                for (int i = 0; i < Dimension; i++) D1[i] += h * h * F[i];
            }

            for (int i = 0; i < Dimension; i++) Y1[i] += D1[i];
            F = Evaluate(X + DeltaX, Y1);
            for (int i = 0; i < Dimension; i++) D1[i] = D1[i] / h + h * F[i] / 2.0;

            YP1 = D1;

        }

        private static int[] N = new int[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    }

    internal abstract class MultiOdeStepper : BaseOdeStepper<IList<double>> {

        public MultiOdeStepper (Func<double, IList<double>, IList<double>> rhs, double x0, IList<double> y0, EvaluationSettings settings) : base(rhs, x0, y0, settings) {
            this.Dimension = y0.Count;
        }

        /// <summary>
        /// The number of variables.
        /// </summary>
        public int Dimension { get; private set; }

        // Override the Evaluate function to ensure that the rhs function is given a read-only vector.

        protected override IList<double> Evaluate (double x, IList<double> y) {
            return (base.Evaluate(x, new ReadOnlyCollection<double>(y)));
        }

    }


    internal class RungeKutta54Stepper : MultiOdeStepper {

        public RungeKutta54Stepper (Func<double, IList<double>, IList<double>> rhs, double x0, IList<double> y0, EvaluationSettings settings) : base(rhs, x0, y0, settings) {
            YPrime = Evaluate(X, Y); 
            DeltaX = 1.0;
        }

        // These Dormand-Prince 5(4) parameters are shown to be particularly good in
        // Lawrence Shampine, "Some Practical Runge-Kutta Formulas", Mathematics of Computation 46 (1986) 135-150
        // They are quoted by NR 3rd edition

        private static double[][] a = new double[][] {
            new double[] { },
            new double[] { 1.0 / 5.0 },
            new double[] { 3.0 / 40.0, 9.0 / 40.0 },
            new double[] { 44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0 },
            new double[] { 19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0 },
            new double[] { 9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0 },
            new double[] { 35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0 }
        };

        /*
        private static double[][] a = new double[][] {
            new double[] { },
            new double[] { 1.0 / 20.0 },
            new double[] { -7161.0 / 1024000.0, 116281.0 / 1024000.0 },
            new double[] { 1023.0 / 25600.0, 0.0, 3069.0 / 25600.0 },
            new double[] { 4202367.0 / 11628100.0, 0.0, -3899844.0 / 2907025.0, 3982992.0 / 2907025.0 },
            new double[] { 5611.0 / 114400.0, 0.0, 0.0, 31744.0 / 135025.0, 923521.0 / 5106400.0 },
            new double[] { 21173.0 / 343200.0, 0.0, 0.0, 8602624.0 / 76559175.0, -26782109.0 / 689364000.0, 5611.0 / 283500.0 },
            new double[] { -1221101821869329.0 / 690812928000000.0, 0.0, 0.0, -125.0 / 2.0, -1024030607959889.0 / 168929280000000.0, 1501408353528689.0 / 265697280000000.0, 6070139212132283.0 / 92502016000000.0 }
            // more
        };
        */

        private static double[] c = new double[] { 0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0 };

        /*
        private static double[] c = new double[] {
            0.0, 1.0 / 20.0, 341.0 / 3200.0, 1023.0 / 6400.0, 39.0 / 100.0, 93.0 / 200.0,
            31.0 / 200.0, 943.0 / 1000.0, 7067558016280.0 / 7837150160667.0, 909.0 / 1000.0, 47.0 / 50.0, 1.0, 1.0
        };
        */

        /*
        private static double[] c = new double[] {
            0.0,
            0.020408163265306122448,
            0.088132939149981030086,
            0.13219940872497154512,
            0.42857142857142857142,
            0.53647553922432876813,
            0.22542922268043313662,
            0.63492063492063492063,
            0.47619047619047618947,
            1.05555555555555555555,
            0.77777777777777777777,
            0.14741696242609947604,
            0.9375, 0.975, 1.0, 1.0
        };
        */

        private static double[] b1 = new double[] { 35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0 };

        /*
        private static double[] b1 = new double[] {
            44901867737754616851973.0 / 1014046409980231013380680.0,
            0.0, 0.0, 0.0, 0.0,
            791638675191615279648100000.0 / 2235604725089973126411512319.0,
            3847749490868980348119500000.0 / 15517045062138271618141237517.0,
            -13734512432397741476562500000.0 / 875132892924995907746928783.0,
            12274765470313196878428812037740635050319234276006986398294443554969616342274215316330684448207141.0 / 489345147493715517650385834143510934888829280686609654482896526796523353052166757299452852166040.0,
            -9798363684577739445312500000.0 / 308722986341456031822630699.0,
            282035543183190840068750.0 / 12295407629873040425991.0,
            -306814272936976936753.0 / 1299331183183744997286.0,
            0.0
        };
        */

        /*
        private static double[] b1 = new double[] {
            0.041535560088059591688, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
           -0.42522808741698055893,
            0.49112696294176088418,
            0.45417824175884742543,
            1.00603264909442806518,
            0.23969807142877259184,
           -4.45549129773140746622,
            9.28897877570610182445,
           -6.41006164510035158839,
            10.0 / 13.0
        };
        */

        private static double[] b2 = new double[] { 5179.0 / 57600.0 , 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0 };

        /*
        private static double[] b2 = new double[] {
            10835401739407019406577.0 / 244521829356935137978320.0,
            0.0, 0.0, 0.0, 0.0,
            13908189778321895491375000.0 / 39221135527894265375640567.0,
            73487947527027243487625000.0 / 296504045773342769773399443.0,
            68293140641257649609375000.0 / 15353208647806945749946119.0,
            22060647948996678611017711379974578860522018208949721559448560203338437626022142776381.0 / 1111542009262325874512959185795727215759010577565736079641376621381577236680929558640.0,
            -547971229495642458203125000.0 / 23237214025700991642563601.0,
            0.0,
            0.0,
            -28735456870978964189.0 / 79783493704265043693.0
        };
        */

        // Tsitouras, "Optimized explicit Runge-Kutta pairs of orders 9(8)", Applied Numerical Mathematics 38 (2001) 123

        // Verner's 8(7) coefficients http://people.math.sfu.ca/~jverner/

        public override void Step () {

            Debug.Assert(a.Length == b1.Length);
            Debug.Assert(b1.Length == c.Length);

            // Compute the RK formula.
            //double[][] z = new double[b1.Length][];
            //for (int k = 0; k < Dimension; k++) z[0][k] = DeltaX * YPrime[k];
            double[] z = new double[b1.Length];
            z[0] = DeltaX * YPrime[0];

            for (int i = 1; i < a.Length; i++) {

                double x = X + c[i] * DeltaX;

                //double[] y = new double[Dimension];
                //for (int k = 0; k < Dimension; k++) {
                //    y[k] = Y[k];
                //    for (int j = 0; j < a[i].Length; j++) {
                //        y[k] += a[i][j] * z[j][k];
                //    }
                //}

                double y  = Y[0];
                for (int j = 0; j < a[i].Length; j++) {
                    y += a[i][j] * z[j];
                }

                //IList<double> FY = Evaluate(x, Y);
                //for (int k = 0; k < Dimension; k++) z[i][k] = DeltaX * FY[k];

                z[i] = DeltaX * Evaluate(x, new double[] { y })[0];

            }

            // Compute the new values.
            double dy1 = 0.0;
            double dy2 = 0.0;
            for (int i = 0; i < b1.Length; i++) {
                dy1 += b1[i] * z[i];
                dy2 += b2[i] * z[i]; 
            }
            double y1 = Y[0] + dy1;
            double err = Math.Abs(dy1 - dy2);

            // Check whether our error is within the required tolerance.
            double tol = Settings.ComputePrecision(y1);
            if (err <= tol) {
                X += DeltaX;
                Y[0] = y1;
                YPrime = Evaluate(X, Y);
            }

            // Adjust the step-size.
            double fac = 0.9375 * Math.Pow(tol / err, 1.0 / 5.0);
            if (fac < 0.2) fac = 0.2;
            if (fac > 5.0) fac = 5.0;
            DeltaX *= fac;
                
        }

    }

    public class NevilleExtrapolator {

        public NevilleExtrapolator (int capacity) {
            xValues = new double[capacity];
        }

        int count = 0;
        private double[] xValues;
        private double[] previousRow, row;

        public void Add (double x, double y) {

            xValues[count] = x;

            count++;

            previousRow = row;

            row = new double[count];

            row[0] = y;

            for (int j = 1; j < row.Length; j++) {
                row[j] = row[j - 1] + (row[j - 1] - previousRow[j - 1]) * x / (xValues[count - 1 - j] - x);
            }

        }

        public int Count {
            get {
                return (count);
            }
        }

        public UncertainValue Estimate {
            get {
                if (previousRow == null) throw new InvalidOperationException();
                double value = row[row.Length - 1];
                double uncertainty = Math.Max(Math.Abs(row[row.Length - 1] - row[row.Length - 2]), Math.Abs(row[row.Length - 1] - previousRow[previousRow.Length - 1]));
                return (new UncertainValue(value, uncertainty));
            }
        }

        public void Clear () {
            count = 0;
            previousRow = null;
            row = null;
        }

    }


    internal class BulrischStoerStepper : MultiOdeStepper {

        public BulrischStoerStepper (Func<double, IList<double>, IList<double>> rhs, double x0, IList<double> y0, EvaluationSettings settings) : base(rhs, x0, y0, settings) {
            YPrime = Evaluate(x0, y0);
            //DeltaX = 1.0;
            DeltaX = 0.1;
            //ComputeInitialSetp();

            extrapolators = new NevilleExtrapolator[Dimension];
            for (int i = 0; i < Dimension; i++) extrapolators[i] = new NevilleExtrapolator(N.Length);
            errorExtrapolator = new NevilleExtrapolator(N.Length);

        }


        private void ComputeInitialSetp () {

            double yNorm = 0.0;
            double yPrimeNorm = 0.0;
            for (int i = 0; i < YPrime.Count; i++) {
                yNorm += MoreMath.Sqr(Y[i]);
                yPrimeNorm += MoreMath.Sqr(YPrime[i]);
            }
            yNorm = Math.Sqrt(yNorm);
            yPrimeNorm = Math.Sqrt(yPrimeNorm);

            DeltaX = Settings.ComputePrecision(yNorm) / yPrimeNorm;
            for (int k = 0; k < N.Length / 2; k++) DeltaX *= MoreMath.Sqr(N[k]);

        }

        private NevilleExtrapolator[] extrapolators;

        private NevilleExtrapolator errorExtrapolator;

        private void ClearExtrapolation () {
            for (int i = 0; i < extrapolators.Length; i++) extrapolators[i].Clear();
            errorExtrapolator.Clear();
        }

        private void AddExtrapolationPoint (double controlValue, double[] y) {
            for (int i = 0; i < extrapolators.Length; i++) extrapolators[i].Add(controlValue, y[i]);
        }

        private void PerformExtrapolation (out double[] value, out double norm, out double error) {
            value = new double[extrapolators.Length];
            norm = 0.0;
            error = 0.0;
            for (int i = 0; i < extrapolators.Length; i++) {
                UncertainValue estimate = extrapolators[i].Estimate;
                value[i] = estimate.Value;
                //norm += MoreMath.Sqr(estimate.Value);
                //error += MoreMath.Sqr(estimate.Uncertainty);
                double aValue = Math.Abs(estimate.Value);
                if (aValue > norm) norm = aValue;
                if (estimate.Uncertainty > error) error = estimate.Uncertainty;
            }
            //norm = Math.Sqrt(norm);
            //error = Math.Sqrt(error);
        }

        public override void Step () {

            ClearExtrapolation();
            //NevilleExtrapolator[] extrapolators = new NevilleExtrapolator[Dimension];
            //for (int i = 0; i < Dimension; i++) extrapolators[i] = new NevilleExtrapolator(N.Length);
            //NevilleExtrapolator extrapolator = new NevilleExtrapolator(N.Length);
            //PolynomialInterpolator errorExtrapolator = new PolynomialInterpolator(N.Length);

            // Keep track of the total work (in RHS evaluations) to get to the kth column. 
            int work = 0;

            // We will use the total work and the computed expansion (or contraction) factor
            // to target convergence in a given column to compute the most efficient (in terms
            // of distance stepped per evaluation) column in which to target convergence
            // on the next step. We need some variables to keep track of this.
            double bestEfficiency = 0.0;
            double bestFactor = 1.0; 
            int bestK = -1;

            double previousError = 0.0;

            // Iterate over sub-steps.
            for (int k = 0; k < N.Length; k++) {

                double[] y1 = TrialStep(N[k]);
                AddExtrapolationPoint(MoreMath.Sqr(1.0 / N[k]), y1);
                //for (int i = 0; i < Dimension; i++) extrapolators[i].Add(MoreMath.Sqr(1.0 / N[k]), y1[i]);
                //double y1 = TrialStep(N[k])[0];
                //extrapolator.Add(MoreMath.Sqr(1.0 / N[k]), y1);

                work += N[k];

                // We need at least two trial steps in order for our estimate to have an associated uncertainty.
                if (k < 1) continue;
                //Debug.Assert(extrapolator.Count >= 2);

                double[] yExtrap;
                double norm, error;
                PerformExtrapolation(out yExtrap, out norm, out error);
                //double[] yExtrap = new double[Dimension];
                //double[] yExtrapError = new double[Dimension];
                //for (int i = 0; i < Dimension; i++) {
                //    UncertainValue estimate = extrapolators[i].Estimate;
                //    yExtrap[i] = estimate.Value;
                //    yExtrapError[i] = estimate.Uncertainty;
                //}
                //double norm = Blas1.dNrm2(yExtrap, 0, 1, Dimension);
                //double error = Blas1.dNrm2(yExtrapError, 0, 1, Dimension);
                double tol = Settings.ComputePrecision(norm);

                //UncertainValue estimate = extrapolator.Estimate;
                //double tol = Settings.ComputePrecision(estimate.Value);

                // Compute the expansion factor required to target convergence in this column on the next step.
                // Compute the corresponding efficiency and determine if this is the best column to target.
                double factor = Math.Pow(tol / error, 1.0 / (2 * k + 1));
                double efficiency = factor / work;
                //errorExtrapolator.Add(k, Math.Log(estimate.Uncertainty));
                // In give ourselves a margin of safety, we never want to target the final column.
                if (((k + 1) < N.Length) && (efficiency > bestEfficiency)) {
                    bestEfficiency = efficiency;
                    bestFactor = factor;
                    bestK = k;
                }

                // Check for convergence.
                if (error <= tol) {

                    X += DeltaX;
                    Y = yExtrap;
                    //Y[0] = estimate.Value;
                    YPrime = Evaluate(X, Y);

                    // It would be nice to try to figure out whether we should target convergence
                    // in a higher column on the next step. There are a few ways to go about this.
                    // (1) We could go ahead and do the next trial step; this is very expensive.
                    // (2) We could use polynomial extrapolation to predict the next error; this
                    // works well and isn't too expensive, but it's a good bit of extra code and
                    // memory dedicated to a rather marginal case. (3) We could use the estimate
                    // that NR quotes from Hairer, Noersett, and Wanner:
                    //   {\rm err}_{k+1} \approx \left \frac{n_0}{n_{k+1}} \right)^2 {\rm err}_{k}
                    // This estimate is order-of-magnitude as best, but it's very simple and in
                    // the end what I chose to use.

                    // To give ourselves a margin of safety, we never want to target the last
                    // column, so we'll only even try this if we are converging at least two
                    // columns from the last.
                    if ((k + 2) < N.Length) {

                        double extrapolatedError = error * (error / previousError);
                        //double extrapolatedError = MoreMath.Sqr(1.0 * N[0] / N[k + 1]) * error;
                        //double extrapolatedError = Math.Exp(errorExtrapolator.Evaluate(k + 1));
                        int extrapolatedWork = work + N[k + 1];
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

                previousError = error;
                    
            }

            // Adjust the step size to target convergence in the optimum column in the next step.
            if (bestFactor < 0.2) bestFactor = 0.2;
            if (bestFactor > 5.0) bestFactor = 5.0;
            DeltaX *= (0.9375 * bestFactor);

        }

        private static readonly int[] N = new int[] { 2, 4, 6, 8, 10, 12, 14, 16, 18 };

        // do a step consisting of n mini-steps

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
            
            return(nY1);
            

            /*
            double Y0 = Y[0];

            double Y1 = Y0 + h * YPrime[0];

            for (int k = 1; k < n; k++) {
                double Y2 = Y0 + 2.0 * h * Evaluate(X + k * h, new double[] { Y1 })[0];
                Y0 = Y1;
                Y1 = Y2;
            }

            Y1 = (Y1 + Y0 + h * Evaluate(X + DeltaX, new double[] { Y1 })[0]) / 2.0;

            return (Y1);
            */
        }

    }

    internal abstract class OdeStrategy<T> {

        public abstract double ComputeInitialStep (BaseOdeStepper<T> stepper);

        public abstract void Step (BaseOdeStepper<T> stepper);

    }

    internal abstract class BulrichStoerStrategy<T> : OdeStrategy<T> {

        protected abstract double Norm (IList<double> values);

        public override double ComputeInitialStep (BaseOdeStepper<T> stepper) {



            throw new NotImplementedException();
        }

        public override void Step (BaseOdeStepper<T> stepper) {
            throw new NotImplementedException();
        }

        protected abstract int StepCount (int k);

        protected abstract int KMax { get; }

    }

}
