using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Analysis {

    

    public static partial class MultiFunctionMath {

        /// <summary>
        /// Finds a local minimum of a multi-dimensional function in the vincinity of the given starting location.
        /// </summary>
        /// <param name="function">The multi-dimensional function to minimize.</param>
        /// <param name="start">The starting location for the search.</param>
        /// <returns>The local minimum.</returns>
        public static MultiExtremum FindMinimum (Func<IList<double>, double> function, IList<double> start) {
            if (function == null) throw new ArgumentNullException("function");
            if (start == null) throw new ArgumentNullException("start");

            EvaluationSettings settings = GetDefaultOptimizationSettings(start.Count);
            //EvaluationSettings settings = new EvaluationSettings() { RelativePrecision = 1.0E-12, AbsolutePrecision = 1.0E-14, EvaluationBudget = 32 * (d + 1) * (d + 2) };
            return (FindMinimum(function, start, settings));
        }

        private static EvaluationSettings GetDefaultOptimizationSettings (int d) {
            EvaluationSettings settings = new EvaluationSettings();
            settings.RelativePrecision = Math.Pow(10.0, -(8.0 + 8.0 / d));
            settings.AbsolutePrecision = settings.RelativePrecision;
            settings.EvaluationBudget = 16 * (d + 1) * (d + 2) * (d + 3);
            return (settings);
        }

        /// <summary>
        /// Finds a local minimum of a multi-dimensional function in the vincinity of the given starting location, subject to the given evaluation constraints.
        /// </summary>
        /// <param name="function">The multi-dimensional function to minimize.</param>
        /// <param name="start">The starting location for the search.</param>
        /// <param name="settings">The evaluation settings that govern the search for the minimum.</param>
        /// <returns>The local minimum.</returns>
        public static MultiExtremum FindMinimum (Func<IList<double>, double> function, IList<double> start, EvaluationSettings settings) {
            if (function == null) throw new ArgumentNullException("function");
            if (start == null) throw new ArgumentNullException("start");
            if (settings == null) throw new ArgumentNullException("settings");

            MultiFunctor f = new MultiFunctor(function);
            return (FindMinimum_ModelTrust(f, start, 0.25, settings));
        }

        /// <summary>
        /// Finds a local maximum of a multi-dimensional function in the vincinity of the given starting location, subject to the given evaluation constraints.
        /// </summary>
        /// <param name="function">The multi-dimensional function to maximize.</param>
        /// <param name="start">The starting location for the search.</param>
        /// <param name="settings">The evaluation settings that govern the search for the maximum.</param>
        /// <returns>The local maximum.</returns>
        public static MultiExtremum FindMaximum (Func<IList<double>, double> function, IList<double> start, EvaluationSettings settings) {
            MultiFunctor f = new MultiFunctor(function, true);
            throw new NotImplementedException();
        }

        private static MultiExtremum FindMinimum_ModelTrust (MultiFunctor f, IList<double> x, double s, EvaluationSettings settings) {

            // Construct an initial model.
            QuadraticInterpolationModel model = QuadraticInterpolationModel.Construct(f, x, s);
            double trustRadius = s;

            int terminationCount = 0;

            while (f.EvaluationCount < settings.EvaluationBudget) {

                // Find the minimum point of the model within the trust radius
                double[] z = model.FindMinimum(trustRadius);
                double expectedValue = model.Evaluate(z);

                // Evaluate the function at the suggested minimum
                double[] point = model.ConvertPoint(z);
                double value = f.Evaluate(point);

                // Check the termination criteria.
                double delta = model.MinimumValue - value;
                if ((Math.Abs(delta) < settings.AbsolutePrecision) || (Math.Abs(delta) <= Math.Abs(value) * settings.RelativePrecision)) {
                    terminationCount++;
                    if (terminationCount > 2) {
                        return (new MultiExtremum(f.EvaluationCount, settings, point, value, model.GetHessian()));
                    }
                } else {
                    terminationCount = 0;
                }

                // There are now three decisions to be made:
                //   1. How to change the trust radius
                //   2. Whether to accept the new point
                //   3. Which existing point to replace

                // Adjust the trust region radius based on the ratio of the actual reduction to the expected reduction.
                double r = (model.MinimumValue - value) / (model.MinimumValue - expectedValue);
                if (r < 0.25) {
                    // If we achieved less than 25% of the expected reduction, reduce the trust region.
                    trustRadius = trustRadius / 2.0;
                } else if (r > 0.75) {
                    // If we achieved at least 75% of the expected reduction, increase the trust region.
                    trustRadius = 2.0 * trustRadius;
                }

                // If our 
                // find maximum and minimum points
                int iMax = 0; double fMax = model.values[0];
                int iBad = 0; double fBad = model.ComputeBadness(0, z, point, value);
                for (int i = 1; i < model.values.Length; i++) {
                    if (model.values[i] > fMax) { iMax = i; fMax = model.values[i]; }
                    double bad = model.ComputeBadness(i, z, point, value);
                    if (bad > fBad) { iBad = i; fBad = bad; }
                }
                if (value < fMax) {
                    Debug.WriteLine("iMax={0}, iBad={1}", iMax, iBad);
                    //model.ReplacePoint(iMax, point, z, value);
                    model.ReplacePoint(iBad, point, z, value);
                }

            }

            throw new NonconvergenceException();
        }

    }

    // Represents a quadratic multinomial, i.e. all terms up to squares.
    // E.g. in two dimensions 1, x, y, x^2, xy, and y^2.

    internal class QuadraticModel {

        public QuadraticModel (int d) {
            g = new double[d];
            h = new double[d][];
            for (int i = 0; i < d; i++) h[i] = new double[i + 1];
        }


        // The constant coefficient.
        public double f;

        // Linear coefficients.
        public double[] g;

        // Quadratic coefficients.
        public double[][] h;

        public void ShiftOrigin (double[] z) {
            for (int i = 0; i < g.Length; i++) {
                f += g[i] * z[i];
                for (int j = 0; j <= i; j++) {
                    f += z[i] * h[i][j] * z[j];
                    g[i] += h[i][j] * z[j];
                }
                for (int j = i; j < g.Length; j++) {
                    g[i] += h[j][i] * z[j];
                }
            }
        }

        public double Evaluate (double[] z) {
            double f = this.f;
            for (int i = 0; i < g.Length; i++) {
                f += g[i] * z[i];
                for (int j = 0; j <= i; j++) {
                    f += h[i][j] * z[i] * z[j];
                }
            }
            return (f);
        }

        // Finds the minimum of the quadratic multinomial.
        public double[] FindMinimum (double trustRadius) {

            // Our model is y = f + g^T z + z^T H z / 2. Differentiating and setting equal to zero
            // gives g + H z = 0 or z = H^{-1} (-g). We could solve this linear system by direct
            // inversion, but that could take us to a maximum or H could be non-invertable.

            // Instead we use the truncated conjugate gradient method.

            // The location starts at the origin.
            double[] z = new double[g.Length];

            // The step direction starts as the negative gradient.
            double[] u = new double[g.Length];
            for (int i = 0; i < g.Length; i++) u[i] = -g[i];

            // We will need a vector to hold the deficit r = (-g) - H z.
            double[] r = new double[g.Length];

            // At first we need only its magnitude, and since we start with z = 0, this is just |g|.
            double old_r2 = 0.0;
            for (int i = 0; i < g.Length; i++) {
                old_r2 += MoreMath.Sqr(g[i]);
            }
            if (old_r2 == 0.0) return (z);

            // Since we have a perfect quadratic form, we need to loop at most d times
            // to obtain an exact solution.
            for (int k = 0; k < g.Length; k++) {

                // We will step z -> z + \alpha u.

                // First compute the maximum \alpha that will keep us within the trust radius,
                // i.e. \alpha s.t. |z + \alpha u| = \Delta.
                double z2 = 0.0, u2 = 0.0, zu = 0.0;
                for (int i = 0; i < z.Length; i++) {
                    z2 += MoreMath.Sqr(z[i]);
                    zu += z[i] * u[i];
                    u2 += MoreMath.Sqr(u[i]);
                }
                double t = trustRadius * trustRadius - z2;
                Debug.Assert(t >= 0.0);
                double alphaMax = t / (Math.Sqrt(zu * zu + u2 * t) + zu);

                // Now compute \alpha as prescribed by the conjugate gradient method.
                double uhu = 0.0;
                for (int i = 0; i < u.Length; i++) {
                    for (int j = 0; j <= i; j++) {
                        uhu += 2.0 * u[i] * h[i][j] * u[j];
                    }
                }

                // If we have negative curvature or alpha exceed the maximum, use the maximum.
                // Otherwise use the computed value.
                double alpha;
                if (uhu > 0.0) {
                    alpha = old_r2 / uhu;
                    if (alpha > alphaMax) alpha = alphaMax;
                } else {
                    alpha = alphaMax;
                }

                // Take the step.
                for (int i = 0; i < z.Length; i++) z[i] += alpha * u[i];

                // If we stepped to the boundary, terminate.
                if (alpha == alphaMax) return (z);

                // If the step was very small, terminate?

                // Compute the new deficit.
                double r2 = ComputeDeficit(z, ref r);

                // If we don't terminate here, in the next round
                // u = 0 and \delta z = 0/0 = NaN.
                // Can we terminate on some small value test?
                if (r2 == 0.0) return (z);

                // Compute the next direction.
                double beta = r2 / old_r2;
                for (int i = 0; i < u.Length; i++) {
                    u[i] = r[i] + beta * u[i];
                }

                old_r2 = r2;

            }

            return (z);

        }

        private double ComputeDeficit (double[] z, ref double[] r) {
            double r2 = 0.0;
            for (int i = 0; i < r.Length; i++) {
                r[i] = -g[i];
                for (int j = 0; j <= i; j++) {
                    r[i] -= h[i][j] * z[j];
                }
                for (int j = i; j < z.Length; j++) {
                    r[i] -= h[j][i] * z[j];
                }
                r2 += MoreMath.Sqr(r[i]);
            }
            return (r2);
        }

        public void DivideBy(double q) {
            f /= q;
            for (int i = 0; i < g.Length; i++) {
                g[i] /= q;
                for (int j = 0; j <= i; j++) {
                    h[i][j] /= q;
                }
            }
        }

        public void Subtract (double p, QuadraticModel Q) {
            f -= p * Q.f;
            for (int i = 0; i < g.Length; i++) {
                g[i] -= p * Q.g[i];
                for (int j = 0; j <= i; j++) {
                    h[i][j] -= p * Q.h[i][j];
                }
            }
        }

        public double[][] GetHessian () {
            double[][] H = new double[g.Length][];
            for (int i = 0; i < g.Length; i++) {
                H[i] = new double[i + 1];
                for (int j = 0; j < i; j++) {
                    H[i][j] = h[i][j];
                }
                H[i][i] = 2.0 * h[i][i];
            }
            return (H);
        }

    }

    // Represents a quadratic multinomial interpolation on a set
    // of given points.

    internal class QuadraticInterpolationModel {

        private int d;

        public double[] origin;

#if PAST
        // delete these when finished
        private double f0;
        private double[] g;
        private double[][] h;
#endif

        // The Lagrange polynomials for each interpolation point.

        private QuadraticModel[] polynomials;

        // The total interpolating polynomial.
        
        private QuadraticModel total;

        // The interpolation points and the function values at those points.

        private double[][] points;

        public double[] values;

        // We keep track of a few

        private int minValueIndex;

        //private int maxBadnessIndex;

        //public double[] badnesses;

        public double[][] GetHessian () {
            return (total.GetHessian());
        }

        public double MinimumValue {
            get {
                return (values[minValueIndex]);
            }
        }

        private double ComputeBadressAbsolute (double[] point, double f) {
            double s = 0.0;
            for (int i = 0; i < point.Length; i++) {
                s += MoreMath.Sqr(point[i] - origin[i]);
            }
            return (s * (f - MinimumValue));
        }

        private double ComputeBadness (double[] z, double f) {
            double s = 0.0;
            for (int i = 0; i < z.Length; i++) {
                s += MoreMath.Sqr(z[i]);
            }
            return (s * (f - MinimumValue));
        }

        private void Initialize (MultiFunctor f, IList<double> x, double s) {

            // Allocate storage
            d = x.Count;
            int m = (d + 1) * (d + 2) / 2;
            origin = new double[d];
            points = new double[m][];
            for (int i = 0; i < m; i++) points[i] = new double[d];
            values = new double[m];
            polynomials = new QuadraticModel[m];
            for (int i = 0; i < m; i++) polynomials[i] = new QuadraticModel(d);

            // Start with x as the origin.
            x.CopyTo(origin, 0);

            // The first interpolation point is the origin.
            x.CopyTo(points[0], 0);
            values[0] = f.Evaluate(points[0]);

            // Compute 2d more interpolation points one step along each axis. 
            int k = 0;
            for (int i = 0; i < d; i++) {
                k++;
                x.CopyTo(points[k], 0);
                points[k][i] += s;
                double plusValue = f.Evaluate(points[k]);
                values[k] = plusValue;
                k++;
                x.CopyTo(points[k], 0);
                points[k][i] -= s;
                double minusValue = f.Evaluate(points[k]);
                values[k] = minusValue;
            }

            // Compute d(d+1)/2 more interpolation points at the corners.
            for (int i = 0; i < d; i++) {
                for (int j = 0; j < i; j++) {
                    k++;
                    x.CopyTo(points[k], 0);
                    points[k][i] += s;
                    points[k][j] += s;
                    double cornerValue = f.Evaluate(points[k]);
                    values[k] = cornerValue;
                }
            }

            double s1 = 1.0 / s;
            double s2 = s1 * s1;

            // Compute the Lagrange polynomial for each point

            k = 0;

            for (int i = 0; i < d; i++) {
                k++;
                polynomials[2 * i + 1].g[i] = 0.5 * s1;
                polynomials[2 * i + 1].h[i][i] = 0.5 * s2;
                k++;
                polynomials[2 * i + 2].g[i] = -0.5 * s1;
                polynomials[2 * i + 2].h[i][i] = 0.5 * s2;
            }

            for (int i = 0; i < d; i++) {
                for (int j = 0; j < i; j++) {
                    k++;
                    polynomials[k].h[i][j] = s2;
                    polynomials[2 * i + 1].h[i][j] -= s2;
                    polynomials[2 * j + 1].h[i][j] -= s2;
                }
            }

            polynomials[0].f = 1.0;
            for (int l = 1; l < m; l++) {
                for (int i = 0; i < d; i++) {
                    polynomials[0].g[i] -= polynomials[l].g[i];
                    for (int j = 0; j <= i; j++) {
                        polynomials[0].h[i][j] -= polynomials[l].h[i][j];
                    }
                }
            }

            // Compute the total interpolating polynomial.

            total = new QuadraticModel(d);
            for (int l = 0; l < m; l++) {
                total.f += values[l] * polynomials[l].f;
                for (int i = 0; i < d; i++) {
                    total.g[i] += values[l] * polynomials[l].g[i];
                    for (int j = 0; j <= i; j++) {
                        total.h[i][j] += values[l] * polynomials[l].h[i][j];
                    }
                }
            }

            // Find the minimum point.

            minValueIndex = 0;
            for (int i = 1; i < m; i++) {
                if (values[i] < values[minValueIndex]) minValueIndex = i;
            }

            // move the origin to the minimum point
            double[] z = new double[d];
            for (int i = 0; i < z.Length; i++) z[i] = points[minValueIndex][i] - origin[i];
            ShiftOrigin(z);

            // compute badnesses
            //badnesses = new double[points.Length];
            //for (int i = 0; i < points.Length; i++) {
            //    badnesses[i] = ComputeBadressAbsolute(points[i], values[i]);
            //}

        }

#if PAST
        public static QuadraticInterpolationModel Construct (double[][] points, double[] values) {

            int m = points.Length;
            int d = points[0].Length;

            // find the minimum point, use it as the origin
            int iMin = 0; double fMin = values[0];
            for (int i = 1; i < values.Length; i++) {
                if (values[i] < fMin) { iMin = i; fMin = values[i]; }
            }

            SquareMatrix A = new SquareMatrix(m);
            int c = 0;
            for (int r = 0; r < m; r++) {
                A[r, 0] = 1.0;
            }
            for (int i = 0; i < d; i++) {
                c++;
                for (int r = 0; r < m; r++) {
                    A[r, c] = points[r][i] - points[iMin][i];
                }
            }
            for (int i = 0; i < d; i++) {
                for (int j = 0; j <= i; j++) {
                    c++;
                    for (int r = 0; r < m; r++) {
                        A[r, c] = (points[r][i] - points[iMin][i]) * (points[r][j] - points[iMin][j]);
                    }
                }
            }
            ColumnVector b = new ColumnVector(values);

            SquareQRDecomposition QR = A.QRDecomposition();
            ColumnVector a = QR.Solve(b);

            QuadraticInterpolationModel model = new QuadraticInterpolationModel();
            model.d = d;
            model.origin = points[iMin];
            model.f0 = a[0];
            model.g = new double[d];
            c = 0;
            for (int i = 0; i < d; i++) {
                c++;
                model.g[i] = a[c];
            }
            model.h = new double[d][];
            for (int i = 0; i < d; i++) {
                model.h[i] = new double[d];
            }
            for (int i = 0; i < d; i++) {
                for (int j = 0; j <= i; j++) {
                    c++;
                    if (i == j) {
                        model.h[i][j] = 2.0 * a[c];
                    } else {
                        model.h[i][j] = a[c];
                        model.h[j][i] = a[c];
                    }
                }
            }

            return (model);

        }

        public static QuadraticInterpolationModel Construct (Func<IList<double>, double> f, double[] x, double s) {
            MultiFunctor mf = new MultiFunctor(f);
            return (Construct(mf, x, s));
        }
#endif

        internal static QuadraticInterpolationModel Construct (MultiFunctor f, IList<double> x, double s) {
            QuadraticInterpolationModel model = new QuadraticInterpolationModel();
            model.Initialize(f, x, s);
            return (model);
        }

        // Shift to x0' = x0 + z, i.e. measure all x from point x0 + z instead of x0.
        // By multiplying out quadratic form f0 + g^T (x - x0) + (x - x0)^T h (x - x0) with and without primes and demanding equality, we find
        //   f0' = f0 + g^T z + 1/2 z^T h z
        //   g' = g + h z^T
        //   h' = h

        public void ShiftOrigin (double[] z) {
            for (int i = 0; i < origin.Length; i++) origin[i] += z[i];
            for (int i = 0; i < polynomials.Length; i++) polynomials[i].ShiftOrigin(z);
            total.ShiftOrigin(z);

        }

        public double[] FindMinimum (double trustRadius) {
            return (total.FindMinimum(trustRadius));
        }

        public void ReplacePoint (int index, double[] point, double[] z, double value) {

            // There are ~ m Lagrange polynomials, changing each requires ~ d^2 work, so
            // the total work to replace a point is ~ m * d^2.

            // For m ~ d^2, that is ~ d^4 work total. That is a lot, but still a lot
            // less than the direct inversion of a d^2 X d^2 design matrix, which would
            // be ~ d^6.

            // Rescale the polynomial being replaced so that it is 1 at the new point.
            // It is still zero at all the other points.
            double a = polynomials[index].Evaluate(z);
            polynomials[index].DivideBy(a);

            // Subtract from each of the other polynomials its value at the new point times the new polynomial.
            // The result is zero at the new point, and does not affect the values at the other points,
            // because the new polynomial is zero at those points.
            for (int i = 0; i < polynomials.Length; i++) {
                if (i == index) continue;
                double b = polynomials[i].Evaluate(z);
                polynomials[i].Subtract(b, polynomials[index]);
            }

            // Adjust the total polynomial to reflect the changes in all the Lagrange polynomials.
            double c = total.Evaluate(z) - value;
            total.Subtract(c, polynomials[index]);

            // Update the minimum.
            if (value < values[minValueIndex]) {
                minValueIndex = index;
                ShiftOrigin(z);
            }

            points[index] = point;
            values[index] = value;

            //badnesses[index] = ComputeBadressAbsolute(points[index], value);
        }


        public double Evaluate (double[] z) {
            return (total.Evaluate(z));
        }

        public double[] ConvertPoint (double[] z) {
            double[] x = new double[d];
            for (int i = 0; i < d; i++) x[i] = origin[i] + z[i];
            return (x);
        }

        public double ComputeBadness (int index, double[] z, double[] point, double value) {
            double s = 0.0;
            if (value < MinimumValue) {
                for (int i = 0; i < point.Length; i++) {
                    s += MoreMath.Sqr(points[index][i] - point[i]);
                }
            } else {
                if (index == minValueIndex) return (0.0);
                for (int i = 0; i < point.Length; i++) {
                    s += MoreMath.Sqr(points[index][i] - origin[i]);
                }
            }
            s = Math.Pow(s, 3.0 / 2.0) * Math.Abs(polynomials[index].Evaluate(z));
            return (s);
        }

    }


}
