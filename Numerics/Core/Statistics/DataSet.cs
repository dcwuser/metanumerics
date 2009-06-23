using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Diagnostics;

using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    // There is a bit of ugliness here due to the fact that we origially published with the
    // non-generic DataPoint and DataSet. The generic versions were added later, and their
    // introducion is mostly smooth, but in order to be non-breaking DataSet<T> can't be
    // IEnumerable<DataPoint<T>>, because that would make DataSet IEnumerable<DataPoint<double>>,
    // but it's already IEnumerable<DataPoint> and it can't be both. In the future, we should
    // eliminate DataPoint and just use DataPoint<double> and re-jigger enumerablity

    /// <summary>
    /// Represents an experimental data point that is a function of an arbitrary variable.
    /// </summary>
    /// <typeparam name="T">The type of the ordinate (independent variable) characterizing the data point.</typeparam>
    public class DataPoint<T> {

        private T x;
        private UncertainValue y;

        /// <summary>
        /// Initializes a new data point with the given values for the ordinate and uncertain abcissa.
        /// </summary>
        /// <param name="x">The ordinate.</param>
        /// <param name="y">The abcissa.</param>
        public DataPoint (T x, UncertainValue y) {
            this.x = x;
            this.y = y;
        }

        /// <summary>
        /// Initializes a new data point with the given values for the ordinate, abcissa, and uncertainty.
        /// </summary>
        /// <param name="x">The ordinate.</param>
        /// <param name="y">The best estimate of the abcissa.</param>
        /// <param name="dy">The uncertainty in the abcissa.</param>
        public DataPoint (T x, double y, double dy) {
            this.x = x;
            this.y = new UncertainValue(y, dy);
        }

        /// <summary>
        /// Gets or sets the value of the ordinate (independent variable).
        /// </summary>
        public T X {
            get {
                return (x);
            }
        }

        /// <summary>
        /// Gets or sets the uncertain value of the abcissa (the depdent variable).
        /// </summary>
        public UncertainValue Y {
            get {
                return (y);
            }
        }

        // equality

        private static bool Equals (DataPoint<T> d1, DataPoint<T> d2) {
            if (Object.ReferenceEquals(d1, null)) {
                if (Object.ReferenceEquals(d2, null)) {
                    return (true);
                } else {
                    return (false);
                }
            } else {
                if (Object.ReferenceEquals(d2, null)) {
                    return (false);
                } else {
                    return ((d1.X.Equals(d2.X)) && (d1.Y == d2.Y));
                }
            }
        }

        /// <summary>
        /// Determines whether two data points are equal.
        /// </summary>
        /// <param name="d1">The first data point.</param>
        /// <param name="d2">The second data point.</param>
        /// <returns>True if the data points are equal, otherwise false.</returns>
        public static bool operator == (DataPoint<T> d1, DataPoint<T> d2) {
            return (Equals(d1, d2));
        }

        /// <summary>
        /// Determines whether two data points are not equal.
        /// </summary>
        /// <param name="d1">The first data point.</param>
        /// <param name="d2">The second data point.</param>
        /// <returns>True if the data points are not equal, otherwise false.</returns>
        public static bool operator != (DataPoint<T> d1, DataPoint<T> d2) {
            return (!Equals(d1, d2));
        }

        /// <summary>
        /// Determines whether the object represents the same data point.
        /// </summary>
        /// <param name="obj">The object.</param>
        /// <returns>True if the object represents the same data point, otherwise false.</returns>
        public override bool Equals (object obj) {
            DataPoint<T> d = obj as DataPoint<T>;
            if (Object.ReferenceEquals(d, null)) {
                return (false);
            } else {
                return (Equals(this, d));
            }
        }

        /// <summary>
        /// Gets a hash code for the data point.
        /// </summary>
        /// <returns>A hash code for the data point.</returns>
        public override int GetHashCode () {
            return (x.GetHashCode() ^ y.GetHashCode());
        }

    }

    /// <summary>
    /// Represents an experimental data point that is a function of a single real variable.
    /// </summary>
    /// <remarks><para></para></remarks>
    /// <seealso cref="DataSet" />
    public class DataPoint : DataPoint<double> {

        //private double x;
        //private UncertainValue y;

        /// <summary>
        /// Initializes a new data point with the given values for the ordinate and uncertain abcissa.
        /// </summary>
        /// <param name="x">The ordinate.</param>
        /// <param name="y">The abcissa.</param>
        public DataPoint (double x, UncertainValue y) : base(x,y) {
            //this.x = x;
            //this.y = y;
        }

        /// <summary>
        /// Initializes a new data point with the given values for the ordinate, abcissa, and uncertainty.
        /// </summary>
        /// <param name="x">The ordinate.</param>
        /// <param name="y">The best estimate of the abcissa.</param>
        /// <param name="dy">The uncertainty in the abcissa.</param>
        public DataPoint (double x, double y, double dy) : base (x, y, dy) {
            //this.x = x;
            //this.y = new UncertainValue(y, dy);
        }

    }

    /// <summary>
    /// Represents a set of measurements.
    /// </summary>
    /// <typeparam name="T">The type of independent variable associated with each measurement.</typeparam>
    public class DataSet<T> {

        /// <summary>
        /// Initializes a new, empty data set.
        /// </summary>
        public DataSet () {
        }

        private List<DataPoint<T>> data = new List<DataPoint<T>>();

        /// <summary>
        /// Adds a new data point to the set.
        /// </summary>
        /// <param name="datum">The data point.</param>
        public void Add (DataPoint<T> datum) {
            if (datum == null) throw new ArgumentNullException("datum");
            data.Add(datum);
        }

        /// <summary>
        /// Adds a series of data points to the set.
        /// </summary>
        /// <param name="data">The data points.</param>
        public void Add (IEnumerable<DataPoint<T>> data) {
            if (data == null) throw new ArgumentNullException("data");
            foreach (DataPoint<T> datum in data) {
                this.data.Add(datum);
            }
        }

        
        /// <summary>
        /// Removes a data point from the set.
        /// </summary>
        /// <param name="datum">The data point to remove.</param>
        /// <returns>True if the data point was found and removed; otherwise false.</returns>
        public bool Remove (DataPoint<T> datum) {
            return (data.Remove(datum));
        }

        /// <summary>
        /// Determines whether the set contains the given data point.
        /// </summary>
        /// <param name="datum">The data point.</param>
        /// <returns>True if the set contains the given data point, otherwise false.</returns>
        public bool Contains (DataPoint<T> datum) {
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
        public FitResult FitToLinearFunction (Function<T, double>[] functions) {

            if (functions == null) throw new ArgumentNullException("functions");
            if (functions.Length > data.Count) throw new InvalidOperationException();

            // construct the data matrix
            SymmetricMatrix A = new SymmetricMatrix(functions.Length);
            for (int i = 0; i < A.Dimension; i++) {
                //Debug.WriteLine(String.Format("F_{0}(2.0) = {1}", i, functions[i](2.0)));
                for (int j = 0; j <= i; j++) {
                    double Aij = 0.0;
                    for (int k = 0; k < data.Count; k++) {
                        Aij += functions[i](data[k].X) * functions[j](data[k].X) / Math.Pow(data[k].Y.Uncertainty, 2);
                    }
                    A[i, j] = Aij;
                    //Debug.WriteLine(String.Format("A[{0},{1}] = {2}", i, j, A[i, j]));
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

            // do a chi^2 test
            double chi2 = 0.0;
            for (int i = 0; i < data.Count; i++) {
                double f = 0.0;
                for (int j = 0; j < functions.Length; j++) {
                    f += functions[j](data[i].X) * a[j];
                }
                chi2 += Math.Pow((data[i].Y.Value - f) / data[i].Y.Uncertainty, 2);
            }
            TestResult test = new TestResult(chi2, new ChiSquaredDistribution(data.Count - functions.Length));

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
        /// <exception cref="InvalidOperationException">There are more fit function parameters than data points.</exception>
        public FitResult FitToFunction (Function<double[], T, double> function, double[] start) {
            if (function == null) throw new ArgumentNullException("function");
            if (start == null) throw new ArgumentNullException("start");

            // you can't do a fit with less data than parameters
            if (this.Count < start.Length) throw new InvalidOperationException();

            // create a chi^2 fit metric and minimize it 
            FitMetric<T> metric = new FitMetric<T>(this, function);
            SpaceExtremum minimum = FunctionMath.FindMinimum(new Function<double[], double>(metric.Evaluate), start);

            // compute the covariance (Hessian) matrix by inverting the curvature matrix
            SymmetricMatrix A = minimum.Curvature() / 2.0;
            CholeskyDecomposition CD = A.CholeskyDecomposition(); // should not return null if we were at a minimum
            SymmetricMatrix C = CD.Inverse();

            // package up the results and return them
            TestResult test = new TestResult(minimum.Value, new ChiSquaredDistribution(this.Count - minimum.Dimension));
            FitResult fit = new FitResult(minimum.Location(), C, test);
            return (fit);

        }

        // can't expose this yet, because IEnumerable<DataPoint<double>> would conflict
        // with IEnumerable<DataPoint>; eliminating non-generic DataPoint would fix this

        internal IEnumerator<DataPoint<T>> GetEnumerator () {
            return (data.GetEnumerator());
        }

    }


    /// <summary>
    /// Represents a set of <see cref="DataPoint"/> measurements.
    /// </summary>
    public sealed class DataSet : DataSet<double>, ICollection<DataPoint> {

        /// <summary>
        /// Initializes a new, empty data set.
        /// </summary>
        public DataSet () : base() {
        }

        /// <summary>
        /// Initializes a new data set with the specified data.
        /// </summary>
        /// <param name="data">An enumerator over the <see cref="DataPoint" />s to place in the set.</param>
        public DataSet (IEnumerable<DataPoint> data) : base() {
            if (data == null) throw new ArgumentNullException("data");
            foreach (DataPoint datum in data) {
                Add(datum);
            }
        }

        //private List<DataPoint> data = new List<DataPoint>();

        ///// <summary>
        ///// Adds a new data point to the set.
        ///// </summary>
        ///// <param name="datum">The data point.</param>
        //public void Add (DataPoint datum) {
        //    if (datum == null) throw new ArgumentNullException("datum");
        //    data.Add(datum);
        //}

        ///// <summary>
        ///// Adds a series of data points to the set.
        ///// </summary>
        ///// <param name="data">The data points.</param>
        //public void Add (IEnumerable<DataPoint> data) {
        //    if (data == null) throw new ArgumentNullException("data");
        //    foreach (DataPoint datum in data) {
        //        this.data.Add(datum);
        //    }
        //}

        ///// <summary>
        ///// Removes a data point from the set.
        ///// </summary>
        ///// <param name="datum">The data point to remove.</param>
        ///// <returns>True if the data point was found and removed; otherwise false.</returns>
        //public bool Remove (DataPoint datum) {
        //    return (data.Remove(datum));
        //}

        ///// <summary>
        ///// Determines whether the set contains the given data point.
        ///// </summary>
        ///// <param name="datum">The data point.</param>
        ///// <returns>True if the set contains the given data point, otherwise false.</returns>
        //public bool Contains (DataPoint datum) {
        //    return (data.Contains(datum));
        //}

        ///// <summary>
        ///// Removes all data points from the set.
        ///// </summary>
        //public void Clear () {
        //    data.Clear();
        //}

        ///// <summary>
        ///// Gets the size of the data set.
        ///// </summary>
        //public int Count {
        //    get {
        //        return (data.Count);
        //    }
        //}

        void ICollection<DataPoint>.Add (DataPoint datum) {
            Add(datum);
        }

        bool ICollection<DataPoint>.Contains (DataPoint datum) {
            return (Contains(datum));
        }

        bool ICollection<DataPoint>.Remove (DataPoint datum) {
            return (Remove(datum));
        }

        int ICollection<DataPoint>.Count {
            get {
                return (Count);
            }
        }

        bool ICollection<DataPoint>.IsReadOnly {
            get {
                return (false);
            }
        }

        void ICollection<DataPoint>.CopyTo (DataPoint[] array, int offset) {
            int i = offset;
            DataPointEnumerator e = new DataPointEnumerator(GetEnumerator());
            while (e.MoveNext()) {
                array[i] = e.Current;
                i++;
            }
        }

        IEnumerator<DataPoint> IEnumerable<DataPoint>.GetEnumerator () {
            return (new DataPointEnumerator(GetEnumerator()));
        }

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator () {
            return (new DataPointEnumerator(GetEnumerator()));
        }

        /*
        public FitResult FitToFunction (ParameterizedFunction function, double[] startParameters) {
            throw new NotImplementedException();
        }
        */

        /*
        /// <summary>
        /// Fits the data to a linear combination of fit functions.
        /// </summary>
        /// <param name="functions">The component functions.</param>
        /// <returns>A fit result containing the best-fit coefficients of the component functions and a &#x3C7;<sup>2</sup> test
        /// of the quality of the fit.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="functions"/> is null.</exception>
        public FitResult FitToLinearFunction (Function<double,double>[] functions) {

            if (functions == null) throw new ArgumentNullException("functions");
            if (functions.Length > data.Count) throw new InvalidOperationException();

            // construct the data matrix
            SymmetricMatrix A = new SymmetricMatrix(functions.Length);
            for (int i = 0; i < A.Dimension; i++) {
                //Debug.WriteLine(String.Format("F_{0}(2.0) = {1}", i, functions[i](2.0)));
                for (int j = 0; j <= i; j++) {
                    double Aij = 0.0;
                    for (int k = 0; k < data.Count; k++) {
                        Aij += functions[i](data[k].X) * functions[j](data[k].X) / Math.Pow(data[k].Y.Uncertainty, 2);
                    }
                    A[i, j] = Aij;
                    //Debug.WriteLine(String.Format("A[{0},{1}] = {2}", i, j, A[i, j]));
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

            // do a chi^2 test
            double chi2 = 0.0;
            for (int i = 0; i < data.Count; i++) {
                double f = 0.0;
                for (int j = 0; j < functions.Length; j++) {
                    f += functions[j](data[i].X) * a[j];
                }
                chi2 += Math.Pow((data[i].Y.Value - f) / data[i].Y.Uncertainty, 2);
            }
            TestResult test = new TestResult(chi2, new ChiSquaredDistribution(data.Count - functions.Length));

            // return the results
            FitResult fit = new FitResult(a, C, test);
            return (fit);

        }

        */
 
        /// <summary>
        /// Fits the data to a constant value.
        /// </summary>
        /// <returns>A fit result containing the best combined value and a &#x3C7;<sup>2</sup> test of the quality of the fit.</returns>
        /// <remarks><para>This method provides a simple way to </para></remarks>
        public FitResult FitToConstant () {

            // do required sums
            double S = 0.0;
            double Sy = 0.0;
            foreach (DataPoint datum in this) {
            //for (int i = 0; i < data.Count; i++) {
                //DataPoint datum = data[i];
                double y = datum.Y.Value;
                double dy = datum.Y.Uncertainty;
                double w = 1.0/(dy*dy);
                S += w;
                Sy += w * y; 
            }

            // compute best value for m for model x2 = m
            double m = Sy / S;
            double dm = 1.0 / Math.Sqrt(S);

            // compute chi^2
            double chi2 = 0.0;
            foreach (DataPoint datum in this) {
            //for (int i = 0; i < data.Count; i++) {
                //DataPoint datum = data[i];
                double x = datum.X;
                double y = datum.Y.Value;
                double dy = datum.Y.Uncertainty;
                double z = (y - m) / dy;
                chi2 += z * z;
            }
            TestResult test = new TestResult(chi2, new ChiSquaredDistribution(Count - 1));

            return (new FitResult(m, dm, test));
        }

        /// <summary>
        /// Fit the data to a proportionality relationship.
        /// </summary>
        /// <returns>A fit result containing the best-fit proportionality constant parameter and a &#x3C7;<sup>2</sup> test of the
        /// quality of the fit.</returns>
        public FitResult FitToProportionality () {

            // do required sums
            double Sxx = 0.0;
            double Sxy = 0.0;
            foreach (DataPoint datum in this) {
            //for (int i = 0; i < data.Count; i++) {
                //DataPoint datum = data[i];
                double x = datum.X;
                double y = datum.Y.Value;
                double dy = datum.Y.Uncertainty;
                double w = 1.0/(dy * dy);
                Sxx += w * x * x;
                Sxy += w * x * y;
            }

            // fit to x2 = a x1
            double a = Sxy / Sxx;
            double da = 1.0 / Math.Sqrt(Sxx);

            // compute chi^2
            double chi2 = 0.0;
            foreach (DataPoint datum in this) {
            //for (int i = 0; i < data.Count; i++) {
                //DataPoint datum = data[i];
                double x = datum.X;
                double y = datum.Y.Value;
                double dy = datum.Y.Uncertainty;
                double z = (y - a * x) / dy;
                chi2 += z * z;
            }
            TestResult test = new TestResult(chi2, new ChiSquaredDistribution(Count - 1));

            // return results
            return (new FitResult(a, da, test));
        }

        /// <summary>
        /// Fits the data to a line.
        /// </summary>
        /// <returns>A fit result containing the best-fit intercept and slope parameters and a &#x3C7;<sup>2</sup> test of
        /// the quality of the fit.</returns>
        public FitResult FitToLine () {

            if (Count < 2) throw new InvalidOperationException();

            // do required sums
            double S = 0.0;
            double Sx = 0.0;
            double Sy = 0.0;
            double Sxx = 0.0;
            double Sxy = 0.0;
            foreach (DataPoint datum in this) {
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
            foreach (DataPoint datum in this) {
            //for (int i = 0; i < data.Count; i++) {
                //DataPoint datum = data[i];
                double x = datum.X;
                double y = datum.Y.Value;
                double dy = datum.Y.Uncertainty;
                double z = (y - (m * x + b)) / dy;
                chi2 += z * z;
            }
            TestResult test = new TestResult(chi2, new ChiSquaredDistribution(Count-2));

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
            Function<double,double>[] functions = new Function<double,double>[order + 1];
            for (int i = 0; i < functions.Length; i++) {
                int n = i;
                functions[i] = delegate(double x) { return (Math.Pow(x, n)); };
                // the weird juggling of integer variables is necessary to ensure that each delegate gets is own variable
                // containing the power; without this trick, each points to the same variable, which contains the last value
            }

            // call the linear function fitter
            FitResult fit = FitToLinearFunction(functions);
            return (fit);

        }

        /*
        /// <summary>
        /// Fits the data to an arbitrary parameterized function.
        /// </summary>
        /// <param name="function">The fit function.</param>
        /// <param name="start">An initial guess at the parameters.</param>
        /// <returns>A fit result containing the best-fitting function parameters
        /// and a &#x3C7;<sup>2</sup> test of the quality of the fit.</returns>
        /// <exception cref="InvalidOperationException">There are more fit function parameters than data points.</exception>
        public FitResult FitToFunction (Function<double[], double, double> function, double[] start) {
            if (function == null) throw new ArgumentNullException("function");
            if (start == null) throw new ArgumentNullException("start");

            // you can't do a fit with less data than parameters
            if (this.Count < start.Length) throw new InvalidOperationException();

            // create a chi^2 fit metric and minimize it 
            FitMetric metric = new FitMetric(this, function);
            SpaceExtremum minimum = FunctionMath.FindMinimum(new Function<double[], double>(metric.Evaluate), start);

            // compute the covariance (Hessian) matrix by inverting the curvature matrix
            SymmetricMatrix A = minimum.Curvature() / 2.0;
            CholeskyDecomposition CD = A.CholeskyDecomposition(); // should not return null if we were at a minimum
            SymmetricMatrix C = CD.Inverse();

            // package up the results and return them
            TestResult test = new TestResult(minimum.Value, new ChiSquaredDistribution(this.Count - minimum.Dimension));
            FitResult fit = new FitResult(minimum.Location(), C, test);
            return (fit);

        }
        */

    }


    internal class FitMetric<T> {

        public FitMetric (DataSet<T> set, Function<double[], T, double> f) {
            this.set = set;
            this.f = f;
        }

        private DataSet<T> set;

        private Function<double[], T, double> f;

        public double Evaluate (double[] p) {
            //Console.WriteLine("Computing chi^2 with p[0]={0}", p[0]);
            double chi2 = 0.0;
            IEnumerator<DataPoint<T>> e = set.GetEnumerator();
            while (e.MoveNext()) {
                DataPoint<T> point = e.Current;
            //foreach (DataPoint<T> point in set) {
                T x = point.X;
                double fx = f(p, x);
                double y = point.Y.Value;
                double dy = point.Y.Uncertainty;
                double z = (y - fx) / dy;
                //Console.WriteLine("  x={0} f(x)={1} y={2} dy={3} z={4}", x, fx, y, dy, z);
                chi2 += z * z;
            }
            //Console.WriteLine("chi^2={0}", chi2);
            return (chi2);
        }

    }

    internal class DataPointEnumerator : IEnumerator<DataPoint>, IEnumerator {

        internal DataPointEnumerator (IEnumerator<DataPoint<double>> e) {
            this.e = e;
        }

        private IEnumerator<DataPoint<double>> e;

        public DataPoint Current {
            get {
                return (new DataPoint(e.Current.X, e.Current.Y));
            }
        }

        object IEnumerator.Current {
            get {
                return (new DataPoint(e.Current.X, e.Current.Y));
            }
        }

        public bool MoveNext () {
            return(e.MoveNext());
        }

        public void Reset () {
            e.Reset();
        }

        void IDisposable.Dispose () {
            e.Dispose();
        }

    }

    /*
    internal class FitMetric {

        public FitMetric (DataSet set, Function<double[], double, double> f) {
            this.set = set;
            this.f = f;
        }

        private DataSet set;

        private Function<double[], double, double> f;

        public double Evaluate (double[] p) {
            //Console.WriteLine("Computing chi^2 with p[0]={0}", p[0]);
            double chi2 = 0.0;
            foreach (DataPoint point in set) {
                double x = point.X;
                double fx = f(p, x);
                double y = point.Y.Value;
                double dy = point.Y.Uncertainty;
                double z = (y - fx) / dy;
                //Console.WriteLine("  x={0} f(x)={1} y={2} dy={3} z={4}", x, fx, y, dy, z);
                chi2 += z * z;
            }
            //Console.WriteLine("chi^2={0}", chi2);
            return (chi2);
        }

    }
    */

}
