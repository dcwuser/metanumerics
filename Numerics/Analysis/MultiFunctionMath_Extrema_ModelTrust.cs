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
        public static MultiExtremum FindLocalMinimum (Func<IList<double>, double> function, IReadOnlyList<double> start) {
            return (FindLocalMinimum(function, start, new EvaluationSettings()));
        }

        private static void SetDefaultOptimizationSettings (EvaluationSettings settings, int d) {
            if (settings.RelativePrecision < 0.0) settings.RelativePrecision = Math.Pow(10.0, -(10.0 + 4.0 / d));
            if (settings.AbsolutePrecision < 0.0) settings.AbsolutePrecision = Math.Pow(10.0, -(10.0 + 4.0 / d));
            if (settings.EvaluationBudget < 0) settings.EvaluationBudget = 16 * (d + 1) * (d + 2) * (d + 3);
        }

        /// <summary>
        /// Finds a local minimum of a multi-dimensional function in the vincinity of the given starting location, subject to the given evaluation constraints.
        /// </summary>
        /// <param name="function">The multi-dimensional function to minimize.</param>
        /// <param name="start">The starting location for the search.</param>
        /// <param name="settings">The evaluation settings that govern the search for the minimum.</param>
        /// <returns>The local minimum.</returns>
        /// <remarks>
        /// <para>The Hessian (matrix of second derivatives) returned with the minimum is an approximation that is constructed in the course of search. It should be
        /// considered a crude approximation, and may not even be that if the minimum is highly non-quadratic.</para>
        /// <para>If you have a constrained minimization problem, require a high-precision solution, and do not have a good initial guess, consider first feeding
        /// your constrained problem into <see cref="FindGlobalMinimum(Func{IList{double}, double}, IList{Interval}, EvaluationSettings)"/>, which supports constraints but gives relatively lower precision solutions, then
        /// feeding the result of that method into this method, which finds relatively precise solutions but does not support constraints.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="function"/>, <paramref name="start"/>, or <paramref name="settings"/> is <see langword="null"/>.</exception>
        /// <exception cref="NonconvergenceException">The number of function evaluations required exceeded the evaluation budget.</exception>
        public static MultiExtremum FindLocalMinimum (Func<IList<double>, double> function, IReadOnlyList<double> start, EvaluationSettings settings) {
            if (function == null) throw new ArgumentNullException("function");
            if (start == null) throw new ArgumentNullException("start");
            if (settings == null) throw new ArgumentNullException("settings");
            return (FindLocalExtremum(function, start, settings, false));
        }

        /// <summary>
        /// Finds a local maximum of a multi-dimensional function in the vincinity of the given starting location.
        /// </summary>
        /// <param name="function">The multi-dimensional function to maximize.</param>
        /// <param name="start">The starting location for the search.</param>
        /// <returns>The local maximum.</returns>
        public static MultiExtremum FindLocalMaximum (Func<IList<double>, double> function, IReadOnlyList<double> start) {
            return (FindLocalMaximum(function, start, new EvaluationSettings()));
        }

        /// <summary>
        /// Finds a local maximum of a multi-dimensional function in the vincinity of the given starting location, subject to the given evaluation constraints.
        /// </summary>
        /// <param name="function">The multi-dimensional function to maximize.</param>
        /// <param name="start">The starting location for the search.</param>
        /// <param name="settings">The evaluation settings that govern the search for the maximum.</param>
        /// <returns>The local maximum.</returns>
        public static MultiExtremum FindLocalMaximum (Func<IList<double>, double> function, IReadOnlyList<double> start, EvaluationSettings settings) {
            if (function == null) throw new ArgumentNullException("function");
            if (start == null) throw new ArgumentNullException("start");
            if (settings == null) throw new ArgumentNullException("settings");
            return (FindLocalExtremum(function, start, settings, true));
        }

        private static MultiExtremum FindLocalExtremum (Func<IList<double>, double> function, IReadOnlyList<double> start, EvaluationSettings settings, bool negate) {
            MultiFunctor f = new MultiFunctor(function, negate);

            // Pick an initial radius; we need to do this better.
            /*
            double s = Double.MaxValue;
            foreach (double x in start) s = Math.Min((Math.Abs(x) + 1.0 / 8.0) / 8.0, s);
            */

            double s = 0.0;
            foreach (double x in start) s += (Math.Abs(x) + 1.0 / 4.0) / 4.0;
            s = s / start.Count;

            //double s = 0.2;
            Debug.WriteLine("s={0}", s);

            SetDefaultOptimizationSettings(settings, start.Count);

            return (FindMinimum_ModelTrust(f, start, s, settings));
        }

        // This method is due to Powell (http://en.wikipedia.org/wiki/Michael_J._D._Powell), but it is not what
        // is usually called Powell's Method (http://en.wikipedia.org/wiki/Powell%27s_method); Powell
        // developed that method in the 1960s, it was included in Numerical Recipies and is very popular.
        // This is a model trust algorithm developed by Powell in the 2000s. It typically uses many
        // fewer function evaluations, but does more intensive calcuations between each evaluation.

        // This is basically the UOBYQA variant of Powell's new methods. It maintains a quadratic model
        // that interpolates between (d + 1) (d + 2) / 2 points. The model is trusted
        // within a given radius. At each step, it moves to the minimum of the model (or the boundary of
        // the trust region in that direction) and evaluates the function. The new value is incorporated
        // into the model and the trust region expanded or contracted depending on how accurate its
        // prediction of the function value was.

        // Papers on these methods are collected at http://mat.uc.pt/~zhang/software.html#powell_software.
        // The UOBYQA paper is here: http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.28.1756.
        // The NEWUOA paper is here: http://www.damtp.cam.ac.uk/user/na/NA_papers/NA2004_08.pdf.
        // The CONDOR system (http://www.applied-mathematics.net/optimization/CONDORdownload.html) is based on these same ideas.
        // The thesis of CONDOR's author (http://www.applied-mathematics.net/mythesis/index.html) was also helpful.

        // It should be very easy to extend this method to constrained optimization, either by incorporating the bounds into
        // the step limits or by mapping hyper-space into a hyper-cube.

        private static MultiExtremum FindMinimum_ModelTrust (MultiFunctor f, IReadOnlyList<double> x, double s, EvaluationSettings settings) {

            // Construct an initial model.
            QuadraticInterpolationModel model = QuadraticInterpolationModel.Construct(f, x, s);
            double trustRadius = s;

            while (f.EvaluationCount < settings.EvaluationBudget) {

                // Find the minimum point of the model within the trust radius
                double[] z = model.FindMinimum(trustRadius);
                double expectedValue = model.Evaluate(z);

                double deltaExpected = model.MinimumValue - expectedValue;

                // Evaluate the function at the suggested minimum
                double[] point = model.ConvertPoint(z);
                double value = f.Evaluate(point);

                double delta = model.MinimumValue - value;
                double tol = settings.ComputePrecision(value);

                // To terminate, we demand: a reduction, the reduction be small, the reduction be in line with its expected value, that we have run up against trust boundary,
                // and that the gradient is small.
                // I had wanted to demand delta > 0, but we run into some cases where delta keeps being very slightly negative, typically orders of magnitude less than tol,
                // causing the trust radius to shrink in and endless cycle that causes our approximation to ultimately go sour, even though terminating on the original
                // very slightly negative delta would have produced an accurate estimate. So we tolerate this case for now.
                if ((-tol / 4.0 <= delta) && (delta <= tol)) {
                    // We demand that the model be decent, i.e. that the expected delta was within tol of the measured delta.
                    if (Math.Abs(delta - deltaExpected) <= tol) {
                        // We demand that the step not just be small because it ran up against the trust radius. If it ran up against the trust radius,
                        // there is probably more to be hand by continuing.
                        double zm = Blas1.dNrm2(z, 0, 1, z.Length);
                        if (zm < trustRadius) {
                            // Finally, we demand that the gradient be small. You might think this was obvious since z was small, but if the Hessian is not positive definite
                            // the interplay of the Hessian and the gradient can produce a small z even if the model looks nothing like a quadratic minimum.
                            double gm = Blas1.dNrm2(model.GetGradient(), 0, 1, z.Length);
                            if (gm * zm <= tol) {
                                if (f.IsNegated) value = -value;
                                return (new MultiExtremum(f.EvaluationCount, settings, point, value, Math.Max(Math.Abs(delta), 0.75 * tol), model.GetHessian()));
                            }
                        }
                    }
                }


                // There are now three decisions to be made:
                //   1. How to change the trust radius
                //   2. Whether to accept the new point
                //   3. Which existing point to replace

                // If the actual change was very far from the expected change, reduce the trust radius.
                // If the expected change did a good job of predicting the actual change, increase the trust radius.
                if ((delta < 0.25 * deltaExpected) /*|| (8.0 * deltaExpected < delta)*/) {
                    trustRadius = trustRadius / 2.0;
                } else if ((0.75 * deltaExpected <= delta) /*&& (delta <= 2.0 * deltaExpected)*/) {
                    trustRadius = 2.0 * trustRadius;
                }
                // It appears that the limits on delta being too large don't help, and even hurt if made too stringent.

                // Replace an old point with the new point.
                int iMax = 0; double fMax = model.values[0];
                int iBad = 0; double fBad = model.ComputeBadness(0, z, point, value);
                for (int i = 1; i < model.values.Length; i++) {
                    if (model.values[i] > fMax) { iMax = i; fMax = model.values[i]; }
                    double bad = model.ComputeBadness(i, z, point, value);
                    if (bad > fBad) { iBad = i; fBad = bad; }
                }
                if (value < fMax) {
                    Debug.WriteLine("iMax={0}, iBad={1}", iMax, iBad);
                    model.ReplacePoint(iBad, point, z, value);
                }
                // There is some question about how best to choose which point to replace.
                // The largest value? The furthest away? The one closest to new min?

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

        public void DivideBy (double q) {
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

        public double[] GetGradient () {
            return (g);
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

        public double[] GetGradient () {
            return (total.GetGradient());
        }

        public double[][] GetHessian () {
            return (total.GetHessian());
        }

        public double MinimumValue {
            get {
                return (values[minValueIndex]);
            }
        }

        private void Initialize (MultiFunctor f, IReadOnlyList<double> x, double s) {

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

        internal static QuadraticInterpolationModel Construct (MultiFunctor f, IReadOnlyList<double> x, double s) {
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
                s = Math.Pow(s, 3.0 / 2.0) * (values[index] - value);
            } else {
                if (index == minValueIndex) return (0.0);
                for (int i = 0; i < point.Length; i++) {
                    s += MoreMath.Sqr(points[index][i] - origin[i]);
                }
                s = Math.Pow(s, 3.0 / 2.0) * (values[index] - values[minValueIndex]);
            }
            return (s);
        }

        public double FindMinimumSeperation () {
            double min = Double.MaxValue;
            for (int i = 0; i < points.Length; i++) {
                for (int j = 0; j < i; j++) {
                    double s = 0.0;
                    for (int k = 0; k < d; k++) {
                        s += MoreMath.Sqr(points[i][k] - points[j][k]);
                    }
                    s = Math.Sqrt(s);
                    if (s < min) min = s;
                }
            }
            return (min);
        }

    }

    // IList defines CopyTo. IReadOnlyList doesn't. So in order to switch from IList to IReadOnlyList,
    // we define an extension method to implement it.

    internal static class CopyHelper {

        public static void CopyTo<T> (this IReadOnlyList<T> source, T[] target, int startIndex) {
            Debug.Assert(source != null);
            Debug.Assert(target != null);
            Debug.Assert(startIndex >= 0);
            for (int i = 0; i < source.Count; i++) {
                target[startIndex + i] = source[i];
            }
        }

    }
}
