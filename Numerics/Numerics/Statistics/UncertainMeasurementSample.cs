using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Diagnostics;

using Meta.Numerics.Analysis;
using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a set of measurements.
    /// </summary>
    /// <typeparam name="T">The type of independent variable associated with each measurement.</typeparam>
    public class UncertainMeasurementSample<T> : ICollection<UncertainMeasurement<T>>, IEnumerable<UncertainMeasurement<T>>, IEnumerable {

        /// <summary>
        /// Initializes a new, empty data set.
        /// </summary>
        public UncertainMeasurementSample () {
        }

        private List<UncertainMeasurement<T>> data = new List<UncertainMeasurement<T>>();

        /// <summary>
        /// Adds a new data point to the set.
        /// </summary>
        /// <param name="datum">The data point.</param>
        public void Add (UncertainMeasurement<T> datum) {
            if (datum == null) throw new ArgumentNullException("datum");
            data.Add(datum);
        }

        /// <summary>
        /// Adds a new data point to the set.
        /// </summary>
        /// <param name="x">The value of the ordinate (independent variable).</param>
        /// <param name="y">The value of the abcissa (dependent variable).</param>
        /// <param name="dy">The uncertainty of the abcissa (dependent variable).</param>
        public void Add (T x, double y, double dy) {
            data.Add(new UncertainMeasurement<T>(x, y, dy));
        }

        /// <summary>
        /// Adds a series of data points to the set.
        /// </summary>
        /// <param name="data">The data points.</param>
        public void Add (IEnumerable<UncertainMeasurement<T>> data) {
            if (data == null) throw new ArgumentNullException("data");
            foreach (UncertainMeasurement<T> datum in data) {
                this.data.Add(datum);
            }
        }


        /// <summary>
        /// Removes a data point from the set.
        /// </summary>
        /// <param name="datum">The data point to remove.</param>
        /// <returns>True if the data point was found and removed; otherwise false.</returns>
        public bool Remove (UncertainMeasurement<T> datum) {
            return (data.Remove(datum));
        }

        /// <summary>
        /// Determines whether the set contains the given data point.
        /// </summary>
        /// <param name="datum">The data point.</param>
        /// <returns>True if the set contains the given data point, otherwise false.</returns>
        public bool Contains (UncertainMeasurement<T> datum) {
            return (data.Contains(datum));
        }

        /// <summary>
        /// Removes all data points from the set.
        /// </summary>
        public void Clear () {
            data.Clear();
        }

        /// <summary>
        /// Gets the size of the data set.
        /// </summary>
        public int Count {
            get {
                return (data.Count);
            }
        }


        /// <summary>
        /// Fits the data to a linear combination of fit functions.
        /// </summary>
        /// <param name="functions">The component functions.</param>
        /// <returns>A fit result containing the best-fit coefficients of the component functions and a &#x3C7;<sup>2</sup> test
        /// of the quality of the fit.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="functions"/> is null.</exception>
        /// <exception cref="InsufficientDataException">There are fewer data points than fit parameters.</exception>
        public FitResult FitToLinearFunction (Func<T, double>[] functions) {

            if (functions == null) throw new ArgumentNullException("functions");
            if (functions.Length > data.Count) throw new InsufficientDataException();

            // construct the design matrix
            RectangularMatrix A = new RectangularMatrix(data.Count, functions.Length);
            for (int r = 0; r < data.Count; r++) {
                for (int c = 0; c < functions.Length; c++) {
                    A[r, c] = functions[c](data[r].X) / data[r].Y.Uncertainty;
                }
            }

            // construct the right-hand-side
            ColumnVector b = new ColumnVector(data.Count);
            for (int r = 0; r < data.Count; r++) {
                b[r] = data[r].Y.Value / data[r].Y.Uncertainty;
            }

            // Solve the system via QR
            ColumnVector a;
            SymmetricMatrix C;
            QRDecomposition.SolveLinearSystem(A, b, out a, out C);

            /*
            // construct the data matrix
            SymmetricMatrix A = new SymmetricMatrix(functions.Length);
            for (int i = 0; i < A.Dimension; i++) {
                for (int j = 0; j <= i; j++) {
                    double Aij = 0.0;
                    for (int k = 0; k < data.Count; k++) {
                        Aij += functions[i](data[k].X) * functions[j](data[k].X) / Math.Pow(data[k].Y.Uncertainty, 2);
                    }
                    A[i, j] = Aij;
                }
            }

            // construct the rhs
            double[] b = new double[functions.Length];
            for (int i = 0; i < b.Length; i++) {
                b[i] = 0.0;
                for (int j = 0; j < data.Count; j++) {
                    b[i] += data[j].Y.Value * functions[i](data[j].X) / Math.Pow(data[j].Y.Uncertainty, 2);
                }
            }

            // solve the system
            CholeskyDecomposition CD = A.CholeskyDecomposition();
            if (CD == null) throw new InvalidOperationException();
            Debug.Assert(CD != null);
            SymmetricMatrix C = CD.Inverse();
            ColumnVector a = CD.Solve(b);
            */

            // do a chi^2 test
            double chi2 = 0.0;
            for (int i = 0; i < data.Count; i++) {
                double f = 0.0;
                for (int j = 0; j < functions.Length; j++) {
                    f += functions[j](data[i].X) * a[j];
                }
                chi2 += Math.Pow((data[i].Y.Value - f) / data[i].Y.Uncertainty, 2);
            }
            TestResult test = new TestResult("ChiSquare", chi2, TestType.RightTailed, new ChiSquaredDistribution(data.Count - functions.Length));

            // return the results
            FitResult fit = new FitResult(a, C, test);
            return (fit);

        }

        /// <summary>
        /// Fits the data to an arbitrary parameterized function.
        /// </summary>
        /// <param name="function">The fit function.</param>
        /// <param name="start">An initial guess at the parameters.</param>
        /// <returns>A fit result containing the best-fitting function parameters
        /// and a &#x3C7;<sup>2</sup> test of the quality of the fit.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="function"/> or <paramref name="start"/> are null.</exception>
        /// <exception cref="InsufficientDataException">There are fewer data points than fit parameters.</exception>
        public FitResult FitToFunction (Func<double[], T, double> function, double[] start) {
            if (function == null) throw new ArgumentNullException("function");
            if (start == null) throw new ArgumentNullException("start");

            // you can't do a fit with less data than parameters
            if (this.Count < start.Length) throw new InsufficientDataException();

            /*
            Func<IList<double>, double> function0 = (IList<double> x0) => {
                double[] x = new double[x0.Count];
                x0.CopyTo(x, 0);
                return(function(x));
            };
            MultiExtremum minimum0 = MultiFunctionMath.FindMinimum(function0, start);
            */

            // create a chi^2 fit metric and minimize it 
            FitMetric<T> metric = new FitMetric<T>(this, function);
            SpaceExtremum minimum = FunctionMath.FindMinimum(new Func<double[], double>(metric.Evaluate), start);

            // compute the covariance (Hessian) matrix by inverting the curvature matrix
            SymmetricMatrix A = 0.5 * minimum.Curvature();
            CholeskyDecomposition CD = A.CholeskyDecomposition(); // should not return null if we were at a minimum
            SymmetricMatrix C = CD.Inverse();

            // package up the results and return them
            TestResult test = new TestResult("ChiSquare", minimum.Value, TestType.RightTailed, new ChiSquaredDistribution(this.Count - minimum.Dimension));
            FitResult fit = new FitResult(minimum.Location(), C, test);
            return (fit);

        }

        // can't expose this yet, because IEnumerable<DataPoint<double>> would conflict
        // with IEnumerable<DataPoint>; eliminating non-generic DataPoint would fix this

        /// <summary>
        /// Gets an enumerator over the measurements.
        /// </summary>
        /// <returns>An enumerator over the measurements.</returns>
        public IEnumerator<UncertainMeasurement<T>> GetEnumerator () {
            return (data.GetEnumerator());
        }

        bool ICollection<UncertainMeasurement<T>>.IsReadOnly {
            get {
                return (false);
            }
        }

        void ICollection<UncertainMeasurement<T>>.CopyTo (UncertainMeasurement<T>[] array, int offset) {
            if (array == null) throw new ArgumentNullException("array");
            data.CopyTo(array, offset);
        }

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator () {
            return (data.GetEnumerator());
        }

    }


    /// <summary>
    /// Represents a set of <see cref="UncertainMeasurement{Double}"/> measurements.
    /// </summary>
    /// <remarks>
    /// <para>This class adds functionality to the <see cref="UncertainMeasurementSample"/> class which applies when
    /// the independent variable (X variable) is a single real number. This includes fitting to a constant, line,
    /// or polynomial.</para>
    /// </remarks>
    public sealed class UncertainMeasurementSample : UncertainMeasurementSample<double> {

        /// <summary>
        /// Initializes a new, empty data set.
        /// </summary>
        public UncertainMeasurementSample ()
            : base() {
        }

        /// <summary>
        /// Initializes a new data set with the specified data.
        /// </summary>
        /// <param name="data">An enumerator over the <see cref="UncertainMeasurement{Double}" />s to place in the set.</param>
        public UncertainMeasurementSample (IEnumerable<UncertainMeasurement<double>> data)
            : base() {
            if (data == null) throw new ArgumentNullException("data");
            foreach (UncertainMeasurement<double> datum in data) {
                Add(datum);
            }
        }

        /// <summary>
        /// Fits the data to a constant value.
        /// </summary>
        /// <returns>A fit result containing the best combined value and a &#x3C7;<sup>2</sup> test of the quality of the fit.</returns>
        /// <remarks><para>This method provides a simple way to </para></remarks>
        public FitResult FitToConstant () {

            if (Count < 1) throw new InsufficientDataException();

            // do required sums
            double S = 0.0;
            double Sy = 0.0;
            foreach (UncertainMeasurement<double> datum in this) {
                double y = datum.Y.Value;
                double dy = datum.Y.Uncertainty;
                double w = 1.0 / (dy * dy);
                S += w;
                Sy += w * y;
            }

            // compute best value for m for model x2 = m
            double m = Sy / S;
            double dm = 1.0 / Math.Sqrt(S);

            // compute chi^2
            double chi2 = 0.0;
            foreach (UncertainMeasurement<double> datum in this) {
                double y = datum.Y.Value;
                double dy = datum.Y.Uncertainty;
                double z = (y - m) / dy;
                chi2 += z * z;
            }
            TestResult test = new TestResult("ChiSquared", chi2, TestType.RightTailed, new ChiSquaredDistribution(Count - 1));

            return (new FitResult(m, dm, test));
        }

        /// <summary>
        /// Fit the data to a proportionality relationship.
        /// </summary>
        /// <returns>A fit result containing the best-fit proportionality constant parameter and a &#x3C7;<sup>2</sup> test of the
        /// quality of the fit.</returns>
        public FitResult FitToProportionality () {

            if (Count < 1) throw new InsufficientDataException();

            // do required sums
            double Sxx = 0.0;
            double Sxy = 0.0;
            foreach (UncertainMeasurement<double> datum in this) {
                double x = datum.X;
                double y = datum.Y.Value;
                double dy = datum.Y.Uncertainty;
                double w = 1.0 / (dy * dy);
                Sxx += w * x * x;
                Sxy += w * x * y;
            }

            // fit to x2 = a x1
            double a = Sxy / Sxx;
            double da = 1.0 / Math.Sqrt(Sxx);

            // compute chi^2
            double chi2 = 0.0;
            foreach (UncertainMeasurement<double> datum in this) {
                double x = datum.X;
                double y = datum.Y.Value;
                double dy = datum.Y.Uncertainty;
                double z = (y - a * x) / dy;
                chi2 += z * z;
            }
            TestResult test = new TestResult("ChiSquared", chi2, TestType.RightTailed, new ChiSquaredDistribution(Count - 1));

            // return results
            return (new FitResult(a, da, test));
        }

        /// <summary>
        /// Fits the data to a line.
        /// </summary>
        /// <returns>A fit result containing the best-fit intercept and slope parameters and a &#x3C7;<sup>2</sup> test of
        /// the quality of the fit.</returns>
        public FitResult FitToLine () {

            if (Count < 2) throw new InsufficientDataException();

            // do required sums
            double S = 0.0;
            double Sx = 0.0;
            double Sy = 0.0;
            double Sxx = 0.0;
            double Sxy = 0.0;
            foreach (UncertainMeasurement<double> datum in this) {
                //for (int i = 0; i < data.Count; i++) {
                //DataPoint datum = data[i];
                double x = datum.X;
                double y = datum.Y.Value;
                double dy = datum.Y.Uncertainty;
                double w = 1.0 / (dy * dy);
                S += w;
                Sx += w * x;
                Sy += w * y;
                Sxx += w * x * x;
                Sxy += w * x * y;
            }
            double D = S * Sxx - Sx * Sx;

            // fit to x2 = m x1 + b
            double b = (Sy * Sxx - Sx * Sxy) / D;
            double db = Math.Sqrt(Sxx / D);
            double m = (S * Sxy - Sx * Sy) / D;
            double dm = Math.Sqrt(S / D);
            double cov = -Sx / D;

            // compute chi^2
            double chi2 = 0.0;
            foreach (UncertainMeasurement<double> datum in this) {
                //for (int i = 0; i < data.Count; i++) {
                //DataPoint datum = data[i];
                double x = datum.X;
                double y = datum.Y.Value;
                double dy = datum.Y.Uncertainty;
                double z = (y - (m * x + b)) / dy;
                chi2 += z * z;
            }
            TestResult test = new TestResult(chi2, new ChiSquaredDistribution(Count - 2));

            return (new FitResult(b, db, m, dm, cov, test));
        }

        /// <summary>
        /// Fits the data to a polynomial.
        /// </summary>
        /// <param name="order">The order of the polynomial to fit.</param>
        /// <returns>A fit result containg the best-fit polynomial coefficients, in order of ascending power from 0 to <paramref name="order"/>,
        /// and a &#x3C7;<sup>2</sup> test of the quality of the fit.</returns>
        /// <exception cref="InvalidOperationException">There are more polynomial coefficients than data points.</exception>
        public FitResult FitToPolynomial (int order) {

            if (order < 0) throw new ArgumentOutOfRangeException("order");
            if (Count < order) throw new InvalidOperationException();

            // create the functions
            Func<double, double>[] functions = new Func<double, double>[order + 1];
            for (int i = 0; i < functions.Length; i++) {
                int n = i;
                functions[i] = delegate(double x) { return (MoreMath.Pow(x, n)); };
                // the weird juggling of integer variables is necessary to ensure that each delegate gets is own variable
                // containing the power; without this trick, each points to the same variable, which contains the last value
            }

            // call the linear function fitter
            FitResult fit = FitToLinearFunction(functions);
            return (fit);

        }

    }


    internal class FitMetric<T> {

        public FitMetric (UncertainMeasurementSample<T> set, Func<double[], T, double> f) {
            this.set = set;
            this.f = f;
        }

        private UncertainMeasurementSample<T> set;

        private Func<double[], T, double> f;

        public double Evaluate (IList<double> p) {
            double chi2 = 0.0;
            foreach (UncertainMeasurement<T> point in set) {
                T x = point.X;
                // This is temporary, only works if p is double[]
                double fx = f((double[]) p, x);
                double y = point.Y.Value;
                double dy = point.Y.Uncertainty;
                double z = (y - fx) / dy;
                chi2 += z * z;
            }
            return (chi2);
        }

        public double Evaluate (double[] p) {
            double chi2 = 0.0;
            IEnumerator<UncertainMeasurement<T>> e = set.GetEnumerator();
            while (e.MoveNext()) {
                UncertainMeasurement<T> point = e.Current;
                T x = point.X;
                double fx = f(p, x);
                double y = point.Y.Value;
                double dy = point.Y.Uncertainty;
                double z = (y - fx) / dy;
                chi2 += z * z;
            }
            return (chi2);
        }

    }

}
