using System;
using System.Collections.Generic;

using Meta.Numerics.Matrices;


namespace Meta.Numerics.Analysis {


    public static partial class MultiFunctionMath {

        /// <summary>
        /// Solves a set of coupled, conservative second order ordinary differential equation initial value problems.
        /// </summary>
        /// <param name="rhs">The right hand side function.</param>
        /// <param name="x0">The initial value of the independent variable.</param>
        /// <param name="y0">The initial values of the functions.</param>
        /// <param name="yp0">The intial values of the functions' derivatives.</param>
        /// <param name="x1">The final value of the independent variable.</param>
        /// <returns>The solution, including the final value of the functions and their derivatives.</returns>
        public static MultiOdeResult SolveConservativeOde (Func<double, IList<double>, IList<double>> rhs, double x0, IList<double> y0, IList<double> yp0, double x1) {
            EvaluationSettings settings = new EvaluationSettings() { RelativePrecision = 1.0E-12, AbsolutePrecision = 1.0E-24, EvaluationBudget = 4096 };
            return (SolveConservativeOde(rhs, x0, y0, yp0, x1, settings));
        }

        /// <summary>
        /// Solves a set of coupled, conservative second order ordinary differential equation initial value problems using the given settings.
        /// </summary>
        /// <param name="rhs">The right hand side function.</param>
        /// <param name="x0">The initial value of the independent variable.</param>
        /// <param name="y0">The initial values of the functions.</param>
        /// <param name="yp0">The intial values of the functions' derivatives.</param>
        /// <param name="x1">The final value of the independent variable.</param>
        /// <param name="settings">The settings to use when solving the problem.</param>
        /// <returns>The solution, including the final value of the functions and their derivatives.</returns>
        public static MultiOdeResult SolveConservativeOde (Func<double, IList<double>, IList<double>> rhs, double x0, IList<double> y0, IList<double> yp0, double x1, EvaluationSettings settings) {

            if (rhs == null) throw new ArgumentNullException(nameof(rhs));
            if (y0 == null) throw new ArgumentNullException(nameof(y0));
            if (yp0 == null) throw new ArgumentNullException(nameof(yp0));
            if (y0.Count != yp0.Count) throw new DimensionMismatchException();
            if (settings == null) throw new ArgumentNullException("settings");

            MultiOdeStepper stepper = new MultiStoermerStepper(rhs, x0, y0, yp0, settings);
            stepper.Integrate(x1);
            return ((MultiOdeResult) stepper.GetResult());
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
        public static MultiOdeResult SolveOde (Func<double, IList<double>, IList<double>> rhs, double x0, IList<double> y0, double x1) {
            EvaluationSettings settings = new EvaluationSettings() { RelativePrecision = 1.0E-12, AbsolutePrecision = 1.0E-24, EvaluationBudget = 4096 };
            return (SolveOde(rhs, x0, y0, x1, settings));
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
        public static MultiOdeResult SolveOde (Func<double, IList<double>, IList<double>> rhs, double x0, IList<double> y0, double x1, EvaluationSettings settings) {

            if (rhs == null) throw new ArgumentNullException("rhs");
            if (y0 == null) throw new ArgumentNullException("y0");
            if (settings == null) throw new ArgumentNullException("settings");

            MultiOdeStepper stepper = new BulrischStoerStepper(rhs, x0, y0, settings);
            stepper.Integrate(x1);

            return ((MultiOdeResult) stepper.GetResult());
        }

    }

}
