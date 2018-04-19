using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Text;

using Meta.Numerics;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Analysis {

    public static partial class FunctionMath {

        /// <summary>
        /// Solves a conservative second order ordinary differential equation initial value problem.
        /// </summary>
        /// <param name="rhs">The right hand side function.</param>
        /// <param name="x0">The initial value of the independent variable.</param>
        /// <param name="y0">The initial value of the function variable.</param>
        /// <param name="yPrime0">The initial value of the function derivative.</param>
        /// <param name="x1">The final value of the independent variable.</param>
        /// <returns>The solution, including the final value of the function and its derivative.</returns>
        /// <exception cref="ArgumentNullException">The <paramref name="rhs"/> is null.</exception>
        /// <exception cref="NonconvergenceException">The ODE could not be integrated to the required precision
        /// before exhausting the maximum allowed number of <paramref name="rhs"/> evaluations.</exception>
        /// <remarks>
        /// <para>For information on integrating conservative ODEs, see
        /// <see cref="IntegrateConservativeOde(Func{double, double, double}, double, double, double, double, OdeSettings)"/>.</para>
        /// <para>This overload uses default values for precision and evaluation budget. It targets a relative precision of
        /// about 10<sup>-12</sup> and an absolute precision of about 10<sup>-24</sup>, with an evaluation budget of about 8000.
        /// </para>
        /// </remarks>
        public static OdeResult IntegrateConservativeOde (Func<double, double, double> rhs, double x0, double y0, double yPrime0, double x1) {
            return (IntegrateConservativeOde(rhs, x0, y0, yPrime0, x1, new OdeSettings()));
        }


        /// <summary>
        /// Solves a conservative second order ordinary differential equation initial value problem using the given settings.
        /// </summary>
        /// <param name="rhs">The right hand side function.</param>
        /// <param name="x0">The initial value of the independent variable.</param>
        /// <param name="y0">The initial value of the function variable.</param>
        /// <param name="yPrime0">The initial value of the function derivative.</param>
        /// <param name="x1">The final value of the independent variable.</param>
        /// <param name="settings">The settings to use when solving the problem.</param>
        /// <returns>The solution, including the final value of the function and its derivative.</returns>
        /// <exception cref="ArgumentNullException">The <paramref name="rhs"/> or <paramref name="settings"/> is null.</exception>
        /// <exception cref="NonconvergenceException">The ODE could not be integrated to the required precision before exhausting the maximum allowed number of <paramref name="rhs"/> evaluations.</exception>
        /// <remarks>
        /// <para>A conservative ODE is an ODE of the form</para>
        /// <img src="../images/ConservativeODE.png" />
        /// <para>where the right-hand-side depends only on x and y, not on the derivative y'. ODEs of this form are called conservative because
        /// they exhibit conserved quantities: combinations of y and y' that maintain the same value as the system evolves. Many forms of
        /// Newtonian equations of motion, for example, are conservative ODEs, with conserved quantities such as energy, momentum, and
        /// angular momentum. Our specialized conservative ODE integrator is not only more efficient for conservative ODEs, but does a
        /// better job of maintaining the conserved quantities.</para>
        /// </remarks>
        public static OdeResult IntegrateConservativeOde (Func<double, double, double> rhs, double x0, double y0, double yPrime0, double x1, OdeSettings settings) {
        
            if (rhs == null) throw new ArgumentNullException(nameof(rhs));
            if (settings == null) throw new ArgumentNullException(nameof(settings));

            SetOdeDefaults(settings);

            SingleStoermerEngine engine = new SingleStoermerEngine(rhs, x0, y0, yPrime0, settings);
            BulrischStoerStrategy strategy = new BulrischStoerStrategy(engine);
            strategy.IntegrateTo(x1);
            return (engine.GetResult());

        }

        /// <summary>
        /// Solves an ordinary differential equation initial value problem.
        /// </summary>
        /// <param name="rhs">The right hand side function, which returns the value of the derivative given
        /// the values of the independent variable and the function.</param>
        /// <param name="x0">The initial value of the independent variable.</param>
        /// <param name="y0">The initial value of the function.</param>
        /// <param name="x1">The final value of the independent variable.</param>
        /// <returns>The solution, including the final value of the function and its derivative.</returns>
        /// <exception cref="ArgumentNullException">The <paramref name="rhs"/> is null.</exception>
        /// <exception cref="NonconvergenceException">The ODE could not be integrated to the required precision
        /// before exhausting the maximum allowed number of <paramref name="rhs"/> evaluations.</exception>
        /// <remarks>
        /// <para>For information on integrating ODEs, see
        /// <see cref="IntegrateOde(Func{double, double, double}, double, double, double, OdeSettings)"/>.</para>
        /// <para>This overload uses default values for precision and evaluation budget. It targets a relative precision of
        /// about 10<sup>-12</sup> and an absolute precision of about 10<sup>-24</sup>,
        /// with an evaluation budget of about 8000.
        /// </para>
        /// </remarks>
        /// <seealso href="https://en.wikipedia.org/wiki/Ordinary_differential_equation"/>
        public static OdeResult IntegrateOde (Func<double, double, double> rhs, double x0, double y0, double x1) {
            return(IntegrateOde(rhs, x0, y0, x1, new OdeSettings()));
        }

        /// <summary>
        /// Solves an ordinary differential equation initial value problem.
        /// </summary>
        /// <param name="rhs">The right hand side function, which returns the value of the derivative given
        /// the values of the independent variable and the function.</param>
        /// <param name="x0">The initial value of the independent variable.</param>
        /// <param name="y0">The initial value of the function.</param>
        /// <param name="x1">The final value of the independent variable.</param>
        /// <param name="settings">The settings to use when solving the problem.</param>
        /// <returns>The solution, including the final value of the function and its derivative.</returns>
        /// <exception cref="ArgumentNullException">The <paramref name="rhs"/> or <paramref name="settings"/> is null.</exception>
        /// <exception cref="NonconvergenceException">The ODE could not be integrated to the required precision
        /// before exhausting the maximum allowed number of <paramref name="rhs"/> evaluations.</exception>
        /// <remarks>
        /// <para>An ordinary differential equation (ODE) has the form:</para>
        /// <img src="../images/ODE.png" />
        /// <para>The function specifying the derivative as a function of x and y is called the right-hand-side (RHS).</para>
        /// <para>The integration of an ODE consists of specifying the value of y at some initial x and computing its value
        /// at a different x in accordance with the differential equation. The terms "initial" and "final" are derived from
        /// the common case where the independent variable is time, but the algorithm works whether the independent variable
        /// resents a time, a location, or a completely non-physical quantity, as long as the problem has the form of an ODE.</para>
        /// <para>ODEs involving multiple, coupled dependent variables can be integrated using the
        /// <see cref="MultiFunctionMath.IntegrateOde(Func{double, IReadOnlyList{double}, IReadOnlyList{double}}, double, IReadOnlyList{double}, double, MultiOdeSettings)"/>
        /// method. Higher order ODEs can be integrated by representing them as coupled ODEs in which the zeroth component
        /// is the desired y, the first component is y', the second component is y'', etc. So-called conservative second order
        /// ODEs should be integrated using the
        /// <see cref="FunctionMath.IntegrateConservativeOde(Func{double, double, double}, double, double, double, double, OdeSettings)"/>
        /// method. If your ODE's RHS depends only on x, the problem reduces to a simple integral, which can be solved more rapidly and
        /// accurately using the <see cref="FunctionMath.Integrate(Func{double, double}, Interval, IntegrationSettings)"/> method.
        /// Analytic techniques can also be used to reduce several other kinds of ODEs to simple integrals or lower-order ODEs.</para>
        /// </remarks>
        /// <seealso href="https://en.wikipedia.org/wiki/Ordinary_differential_equation"/>
        public static OdeResult IntegrateOde (Func<double, double, double> rhs, double x0, double y0, double x1, OdeSettings settings) {

            if (rhs == null) throw new ArgumentNullException(nameof(rhs));
            if (settings == null) throw new ArgumentNullException(nameof(settings));

            SetOdeDefaults(settings);

            SingleBulrischStoerEngine engine = new SingleBulrischStoerEngine(rhs, x0, y0, settings);
            BulrischStoerStrategy strategy = new BulrischStoerStrategy(engine);
            strategy.IntegrateTo(x1);
            return (engine.GetResult());

        }

        internal static void SetOdeDefaults (EvaluationSettings settings) {
            if (settings.RelativePrecision < 0) settings.RelativePrecision = 1.0E-12;
            if (settings.AbsolutePrecision < 0) settings.AbsolutePrecision = 1.0E-24;
            if (settings.EvaluationBudget < 0) settings.EvaluationBudget = 8192;
        }

    }


#if FUTURE
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
#endif

    internal class NevilleExtrapolator {

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

        public double Value {
            get {
                return (row[row.Length - 1]);
            }
        }

        public double Error {
            get {
                if (previousRow == null) {
                    return (Math.Abs(Value));
                } else {
                    return (Math.Max(Math.Abs(row[row.Length - 1] - row[row.Length - 2]), Math.Abs(row[row.Length - 1] - previousRow[previousRow.Length - 1])));
                }
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


    internal abstract class SingleOdeEngine : OdeEngine {

        protected SingleOdeEngine (Func<double, double, double> rhs, double x, OdeSettings settings) : base(x) {
            Debug.Assert(rhs != null);
            Debug.Assert(settings != null);
            this.rhs = rhs;
            this.settings = settings;
        }

        private readonly Func<double, double, double> rhs;

        protected readonly OdeSettings settings;

        public double Evaluate (double x, double y) {
            if (count >= settings.EvaluationBudget) throw new NonconvergenceException();
            count++;
            return (rhs(x, y));
        }

        private int count = 0;

        public int EvaluationCount {
            get {
                return (count);
            }
        }

        public abstract OdeResult GetResult ();

        protected override void AcceptStep () {
            base.AcceptStep();
            if (settings.Listener != null) settings.Listener(GetResult());
        }

    }




    internal class SingleStoermerEngine : SingleOdeEngine, IBulrischStoerEngine {

        public SingleStoermerEngine (Func<double, double, double> rhs, double x, double y, double yPrime, OdeSettings settings) : base(rhs, x, settings) {
            this.Y = y;
            this.YPrime = yPrime;
            this.YPrimePrime = this.Evaluate(x, y);
            ComputeInitialStep();
        }

        private void ComputeInitialStep () {
            double x1 = Math.Abs(Y / YPrime);
            if (Double.IsNaN(x1)) x1 = 0.0;
            double x2 = Math.Abs(YPrime / YPrimePrime);
            if (Double.IsNaN(x2)) x2 = 0.0;
            // We could also do \sqrt{ y / y''}, but that's just \sqrt{x_1 x_2},
            // i.e. geometric mean, so it's guaranteed to lie between them.
            if (x1 == 0.0) {
                DeltaX = x2;
            } else if (x2 == 0.0) {
                DeltaX = x1;
            } else {
                DeltaX = Math.Min(x1, x2);
            }
            if (DeltaX == 0.0) DeltaX = 0.5;
        }

        private double Y;

        private double YPrime;

        private double YPrimePrime;

        private readonly static int[] N = new int[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

        public IList<int> Sizes { get { return (N); } }

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

            double yTol = settings.ComputePrecision(Math.Abs(yExtrapolator.Value));
            double yError = yExtrapolator.Error;
            double yRatio = yError / yTol;

            double ypTol = settings.ComputePrecision(Math.Abs(ypExtrapolator.Value));
            double ypError = ypExtrapolator.Error;
            double ypRatio = ypError / ypTol;

            this.Ratio = Math.Max(yRatio, ypRatio);

        }

        public double Ratio { get; private set; }

        protected override void AcceptStep () {
            Y = yExtrapolator.Estimate.Value;
            YPrime = ypExtrapolator.Estimate.Value;
            YPrimePrime = Evaluate(X + DeltaX, Y);
            base.AcceptStep();
        }

        public override OdeResult GetResult () {

            return (new OdeResult(
                this.X,
                this.Y,
                this.YPrime,
                this.EvaluationCount,
                this.settings
            ));

        }

    }

    internal class SingleBulrischStoerEngine : SingleOdeEngine, IBulrischStoerEngine {

        public SingleBulrischStoerEngine (Func<double, double, double> rhs, double x, double y, OdeSettings settings) : base(rhs, x, settings) {
            this.Y = y;
            this.YPrime = this.Evaluate(x, y);
            ComputeInitialStep();
        }

        private void ComputeInitialStep () {
            DeltaX = Math.Abs(Y / YPrime);
            // If DeltaX is too big (e.g. YPrime = 0), then Integrate will set it to integration interval.
            // If DeltaX is too small (e.g. Y = 0), we have no length scale
            if (Double.IsNaN(DeltaX) || (DeltaX == 0.0)) DeltaX = 0.5;
            // Can be NaN. If y = 0 and y' = 0, 0/0 = NaN

            // Surprisingly, the choice of initial step appears to have only a small effect
            // on evaluation count, presumably because we quickly find right step-size.
            // Tried 1.0, 0.1, |y/y'|, and 1.5 |y/y'|. 
        }

        private double Y;

        private double YPrime;

        private static readonly int[] N = new int[] { 2, 4, 6, 8, 10, 12, 14, 16, 18, 20 };

        public virtual IList<int> Sizes { get { return (N); } }

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

        private NevilleExtrapolator extrapolator = new NevilleExtrapolator(N.Length);

        public virtual void Clear () {
            extrapolator.Clear();
        }

        protected override void AcceptStep () {
            Y = extrapolator.Estimate.Value;
            YPrime = Evaluate(X + DeltaX, Y);
            base.AcceptStep();
        }

        public virtual void AddTrialStep () {

            int k = extrapolator.Count;

            double y = TrialStep(N[k]);
            extrapolator.Add(MoreMath.Sqr(1.0 / N[k]), y);

            double norm = Math.Abs(y);
            double tol = settings.ComputePrecision(norm);
            double error = extrapolator.Error;
            this.Ratio = error / tol;

        }

        public double Ratio {
            get; private set;
        }

        public override OdeResult GetResult () {

            return (new OdeResult(
                this.X,
                this.Y,
                this.YPrime,
                this.EvaluationCount,
                this.settings
            ));

        }

    }


}
