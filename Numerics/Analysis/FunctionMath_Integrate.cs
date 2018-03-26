using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Analysis {

    public static partial class FunctionMath {

        // the public API

        /// <summary>
        /// Evaluates a definite integral.
        /// </summary>
        /// <param name="integrand">The function to be integrated.</param>
        /// <param name="start">The lower integration endpoint.</param>
        /// <param name="end">The upper integration endpoint.</param>
        /// <returns>The result of the integral.</returns>
        /// <exception cref="ArgumentNullException">The <paramref name="integrand"/> is <see langword="null"/>.</exception>
        /// <exception cref="NonconvergenceException">The maximum number of function evaluations was exceeded before the integral
        /// could be determined to the required precision.</exception>
        /// <remarks>
        /// <para>By default, integrals are evaluated to a relative precision of about 10<sup>-14</sup>, about two digits short of full
        /// precision, or an absolute precision of about 10<sup>-16</sup>, using a budget of about 5000 evaluations.
        /// To specify different evaluation settings use
        /// <see cref="Integrate(Func{double, double}, double, double, IntegrationSettings)"/>.</para>
        /// <para>See <see cref="Integrate(Func{double, double}, double, double, IntegrationSettings)"/> for detailed remarks on
        /// numerical integration.</para>
        /// </remarks>
        public static IntegrationResult Integrate(Func<double, double> integrand, double start, double end) {
            IntegrationSettings settings = new IntegrationSettings();
            return (Integrate(integrand, start, end, settings));
        }

        /// <summary>
        /// Evaluates a definite integral.
        /// </summary>
        /// <param name="integrand">The function to be integrated.</param>
        /// <param name="range">The range of integration.</param>
        /// <returns>The result of the integral.</returns>
        /// <exception cref="ArgumentNullException">The <paramref name="integrand"/> is <see langword="null"/>.</exception>
        /// <exception cref="NonconvergenceException">The maximum number of function evaluations was exceeded before the integral
        /// could be determined to the required precision.</exception>
        /// <remarks>
        /// <para>By default, integrals are evaluated to a relative precision of about 10<sup>-14</sup>, about two digits short of full
        /// precision, or an absolute precision of about 10<sup>-16</sup>, using a budget of about 5000 evaluations.
        /// To specify different evaluation settings use
        /// <see cref="Integrate(Func{double, double}, Interval, IntegrationSettings)"/>.</para>
        /// <para>See <see cref="Integrate(Func{double, double}, Interval, IntegrationSettings)"/> for detailed remarks on
        /// numerical integration.</para>
        /// </remarks>
        public static IntegrationResult Integrate (Func<double, double> integrand, Interval range) {
            IntegrationSettings settings = new IntegrationSettings();
            return (Integrate(integrand, range, settings));
        }

        internal static IntegrationSettings SetIntegrationDefaults (IntegrationSettings settings) {
            IntegrationSettings result = new IntegrationSettings();
            result.RelativePrecision = (settings.RelativePrecision < 0.0) ? 1.0E-14 : settings.RelativePrecision;
            result.AbsolutePrecision = (settings.AbsolutePrecision < 0.0) ? 1.0E-15 : settings.AbsolutePrecision;
            result.EvaluationBudget = (settings.EvaluationBudget < 0) ? 5000 : settings.EvaluationBudget;
            result.Listener = settings.Listener;
            return (result);
        }

        /// <summary>
        /// Evaluates a definite integral with the given evaluation settings.
        /// </summary>
        /// <param name="integrand">The function to be integrated.</param>
        /// <param name="range">The range of integration.</param>
        /// <param name="settings">The settings which control the evaluation of the integral.</param>
        /// <returns>The result of the integral.</returns>
        /// <exception cref="ArgumentNullException">The <paramref name="integrand"/> is <see langword="null"/>.</exception>
        /// <exception cref="NonconvergenceException">The maximum number of function evaluations was exceeded before the integral
        /// could be determined to the required precision.</exception>
        /// <remarks>
        /// <para>For information, see <see cref="Integrate(Func{double, double}, double, double, IntegrationSettings)"/>.</para>
        /// </remarks>
        public static IntegrationResult Integrate (Func<double, double> integrand, Interval range, IntegrationSettings settings) {
            return (Integrate(integrand, range.LeftEndpoint, range.RightEndpoint, settings));
        }

        /// <summary>
        /// Evaluates a definite integral with the given evaluation settings.
        /// </summary>
        /// <param name="integrand">The function to be integrated.</param>
        /// <param name="start">The left integration endpoint.</param>
        /// <param name="end">The right integration endpoint.</param>
        /// <param name="settings">The settings which control the evaluation of the integral.</param>
        /// <returns>The result of the integral.</returns>
        /// <remarks>
        /// <para>To do integrals over infinite regions, simply set <paramref name="start"/> or <paramref name="end"/>
        /// to <see cref="System.Double.NegativeInfinity"/> or <see cref="System.Double.PositiveInfinity"/>.</para>
        /// <para>Our integrator handles smooth functions extremely efficiently. It handles integrands with
        /// discontinuities or kinks at the price of slightly more evaluations of the integrand.
        /// It can handle oscillatory functions, as long as cancelation between positive and negative regions
        /// is not too severe. It can integrate logarithmic and mild power-law singularities.</para>
        /// <para>Strong power-law singularities will cause the algorithm to fail with a <see cref="NonconvergenceException"/>.
        /// This is unavoidable for essentially any double-precision numerical integrator. Consider, for example,
        /// the integrable singularity x<sup>-1/2</sup>. Since
        /// &#x3B5; = &#x222B;<sub>0</sub><sup>&#x3B4;</sup> x<sup>-1/2</sup> dx = 2 &#x3B4;<sup>1/2</sup>,
        /// points within &#x3B4; &#x223C; 10<sup>-16</sup> of the end-points, which as a close as you can get to
        /// a point in double precision without being on top of it, contribute at the &#x3B5; &#x223C; 10<sup>-8</sup>
        /// level to our integral, well beyond limit that nearly-full double precision requires. Said differently,
        /// to know the value of the integral to &#x3B5; &#x223C; 10<sup>-16</sup> precision, we would need to
        /// evaluate the contributions of points within &#x3B4; &#x223C; 10<sup>-32</sup> of the endpoints,
        /// which is far closer than we can get.</para>
        /// <para>If you need to evaluate an integral with such a strong singularity, try to make an analytic
        /// change of variable to absorb the singularity before attempting numerical integration. For example,
        /// to evaluate I = &#x222B;<sub>0</sub><sup>b</sup> f(x) x<sup>-1/2</sup> dx, substitute y = x<sup>1/2</sup>
        /// to obtain I = 2 &#x222B;<sub>0</sub><sup>&#x221A;b</sup> f(y<sup>2</sup>) dy.</para>
        /// <para>To do multi-dimensional integrals, use
        /// <see cref="MultiFunctionMath.Integrate(Func{IReadOnlyList{double}, double}, IReadOnlyList{Interval}, IntegrationSettings)"/>.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException">The <paramref name="integrand"/> is <see langword="null"/>.</exception>
        /// <exception cref="NonconvergenceException">The maximum number of function evaluations was exceeded before the integral
        /// could be determined to the required precision.</exception>
        public static IntegrationResult Integrate (Func<double,double> integrand, double start, double end, IntegrationSettings settings) {

            if (integrand == null) throw new ArgumentNullException(nameof(integrand));
            if (settings == null) throw new ArgumentNullException(nameof(settings));

            // Deal with right-to-left integrals
            if (end < start) {
                IntegrationResult r = Integrate(integrand, end, start, settings);
                return (new IntegrationResult(-r.Estimate, r.EvaluationCount, r.Settings));
            }

            // Re-map infinite integrals to finite integrals
            if (Double.IsNegativeInfinity(start) && Double.IsPositiveInfinity(end)) {
                // -\infty to +\infty
                // remap to (-\pi/2,\pi/2)
                Func<double, double> f1 = delegate (double t) {
                    double x = Math.Tan(t);
                    return (integrand(x) * (1.0 + x * x));
                };
                return (Integrate(f1, -Global.HalfPI, +Global.HalfPI, settings));
            } else if (Double.IsPositiveInfinity(end)) {
                // finite to +\infty
                // remap to interval (-1,1)
                Func<double, double> f1 = delegate (double t) {
                    double q = 1.0 / (1.0 - t);
                    double x = start + (1.0 + t) * q;
                    return (integrand(x) * 2.0 * q * q);
                };
                return (Integrate(f1, -1.0, +1.0, settings));
            } else if (Double.IsNegativeInfinity(start)) {
                // -\infty to finite
                // remap to interval (-1,1)
                Func<double, double> f1 = delegate (double t) {
                    double q = t + 1.0;
                    double x = end + (t - 1.0) / q;
                    return(integrand(x) * (2.0 / q / q));
                };
                return(Integrate(f1, -1.0, +1.0, settings));
            }

            // Fix settings.
            settings = SetIntegrationDefaults(settings);

            // normal integral over a finite range
            Debug.Assert(end >= start);
            IAdaptiveIntegrator integrator = new GaussKronrodIntegrator(integrand, Interval.FromEndpoints(start, end));
            IntegrationResult result = Integrate_Adaptive(integrator, settings);
            return (result);
        }

        // the drivers

        private static IntegrationResult Integrate_Adaptive (IAdaptiveIntegrator integrator, IntegrationSettings settings) {

            Debug.Assert(integrator != null);
            Debug.Assert(settings != null);

            LinkedList<IAdaptiveIntegrator> list = new LinkedList<IAdaptiveIntegrator>();
            list.AddFirst(integrator);

            int n = integrator.EvaluationCount;

            while (true) {

                // go through the intervals, adding estimates (and errors)
                // and noting which contributes the most error
                // keep track of the total value and uncertainty
                UncertainValue vTotal = new UncertainValue();

                // keep track of which node contributes the most error
                LinkedListNode<IAdaptiveIntegrator> maxNode = null;
                double maxError = 0.0;

                LinkedListNode<IAdaptiveIntegrator> node = list.First;
                while (node != null) {

                    IAdaptiveIntegrator i = node.Value;

                    UncertainValue v = i.Estimate;
                    vTotal += v;

                    if (v.Uncertainty > maxError) {
                        maxNode = node;
                        maxError = v.Uncertainty;
                    }

                    node = node.Next;

                }

                // Inform listeners of our latest result.
                if (settings.Listener != null) {
                    settings.Listener(new IntegrationResult(vTotal, n, settings));
                }

                // if our error is small enough, return
                double tol = settings.ComputePrecision(vTotal.Value);
                if (vTotal.Uncertainty <= tol) {
                    // Don't claim uncertainty significantly less than tol.
                    if (vTotal.Uncertainty < tol / 2.0) vTotal = new UncertainValue(vTotal.Value, tol / 2.0);
                    return (new IntegrationResult(vTotal, n, settings));
                }

                // if our evaluation count is too big, throw
                if (n > settings.EvaluationBudget) throw new NonconvergenceException();

                // Subdivide the interval with the largest error
                IEnumerable<IAdaptiveIntegrator> divisions = maxNode.Value.Divide();
                foreach (IAdaptiveIntegrator division in divisions) {
                    list.AddBefore(maxNode, division);
                    n += division.EvaluationCount;
                }
                list.Remove(maxNode);

            }

        }

    }

    // integrator contracts

    internal interface IIntegrator {

        Func<double, double> Integrand { get; }

        Interval Range { get; }

        UncertainValue Estimate { get; }

        int EvaluationCount { get; }

    }

    internal interface IAdaptiveIntegrator : IIntegrator {

        IList<IAdaptiveIntegrator> Divide ();

        int Generation { get; }

    }

    // Integrator base class
    // All integrators inherit from this, whether adaptive or iterative.

    internal abstract class Integrator : IIntegrator {

        protected Integrator (Func<double, double> integrand, Interval range) {
            Debug.Assert(integrand != null);
            this.f = integrand;
            this.r = range;
        }

        private Func<double, double> f;
        private Interval r;

        public Func<double, double> Integrand {
            get {
                return (f);
            }
        }

        public Interval Range {
            get {
                return (r);
            }
        }

        public int EvaluationCount {
            get {
                return (n);
            }
        }

        private int n = 0;

        protected double Evaluate (double x) {
            n++;
            double y = f(x);
            // Ignore singularities
            if (Double.IsInfinity(y) || Double.IsNaN(y)) y = 0.0;
            return (y);
        }

        public abstract UncertainValue Estimate { get; }

    }

    // Integrators

    internal class GaussKronrodIntegrator : Integrator, IAdaptiveIntegrator {

        public GaussKronrodIntegrator (Func<double, double> integrand, Interval range) : base(integrand, range) {
            Generation = 1;
            Compute();
        }

        // This is a 15-point Gauss-Kronrod role. It is exact for up to 30th order polynomials.
        // There is an embedded 7-point rule (embedded meaning it uses the same points with different weights).
        // Combining them gives a high-order estimate along with an error estimate without requiring additional evaluations.

        // These values appear at https://en.wikipedia.org/wiki/Gauss%E2%80%93Kronrod_quadrature_formula.
        // Other possibilities appear at http://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights/

        private static readonly double[] x = new double[] {
            -0.9914553711208126,
            -0.9491079123427585,
            -0.8648644233597691,
            -0.7415311855993944,
            -0.5860872354676911,
            -0.4058451513773972,
            -0.2077849550078985,
            0.0,
            0.2077849550078985,
            0.4058451513773972,
            0.5860872354676911,
            0.7415311855993944,
            0.8648644233597691,
            0.9491079123427585,
            0.9914553711208126
        };

        private static readonly double[] c1 = new double[] {
            0.0,
            0.1294849661688697,
            0.0,
            0.2797053914892767,
            0.0,
            0.3818300505051189,
            0.0,
            0.4179591836734694,
            0.0,
            0.3818300505051189,
            0.0,
            0.2797053914892767,
            0.0,
            0.1294849661688697,
            0.0
        };

        private static readonly double[] c2 = new double[] {
            0.02293532201052922,
            0.06309209262997855,
            0.1047900103222502,
            0.1406532597155259,
            0.1690047266392679,
            0.1903505780647854,
            0.2044329400752989,
            0.2094821410847278,
            0.2044329400752989,
            0.1903505780647854,
            0.1690047266392679,
            0.1406532597155259,
            0.1047900103222502,
            0.06309209262997855,
            0.02293532201052922
        };

        private void Compute () {

            // we will use these repeatedly
            double mx = Range.Midpoint;
            double w2 = Range.Width / 2.0;

            // evaluate at the points
            double s1 = 0.0;
            double s2 = 0.0;
            for (int i = 0; i < x.Length; i++) {
                double xx = mx + w2 * x[i];
                double fx = Evaluate(xx);
                s1 += fx * c1[i];
                s2 += fx * c2[i];
            }

            double I1 = w2 * s1;
            double I2 = w2 * s2;

            estimate = new UncertainValue(I2, Math.Abs(I2 - I1));

        }

        private UncertainValue estimate;

        public override UncertainValue Estimate {
            get {
                return (estimate);
            }
        }

        public IList<IAdaptiveIntegrator> Divide () {
            GaussKronrodIntegrator i1 = new GaussKronrodIntegrator(Integrand, Interval.FromEndpoints(Range.LeftEndpoint, Range.Midpoint));
            i1.Generation = this.Generation + 1;
            GaussKronrodIntegrator i2 = new GaussKronrodIntegrator(Integrand, Interval.FromEndpoints(Range.Midpoint, Range.RightEndpoint));
            i2.Generation = this.Generation + 1;
            return (new GaussKronrodIntegrator[] { i1, i2 });
        }

        public int Generation { get; private set; }
    }

}