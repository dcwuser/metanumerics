using System;
using System.Collections.Generic;

namespace Meta.Numerics.Analysis {

    public static partial class FunctionMath {

        // the public API

        /// <summary>
        /// Evaluates a definite integral.
        /// </summary>
        /// <param name="integrand">The function to be integrated.</param>
        /// <param name="range">The range of integration.</param>
        /// <returns>A numerical estimate of the given integral.</returns>
        /// <remarks>
        /// <para>Integral values are accurate to within about a digit of full double precision.</para>
        /// <para>To do integrals over infinite regions, simply set the lower bound of the <paramref name="range"/>
        /// to <see cref="System.Double.NegativeInfinity"/> or the upper bound to <see cref="System.Double.PositiveInfinity"/>.</para>
        /// <para>Our numerical integrator uses a Gauss-Kronrod rule that can integrate efficiently,
        /// combined with an adaptive strategy that limits function
        /// evaluations to those regions required to achieve the desired accuracy.</para>
        /// <para>Our integrator handles smooth functions extremely efficiently. It handles integrands with
        /// discontinuities, or discontinuities of derivatives, at the price of slightly more evaluations
        /// of the integrand. It handles oscilatory functions, so long as not too many periods contribute
        /// significantly to the integral. It can integrate logarithmic and mild power-law singularities.
        /// </para>
        /// <para>Strong power-law singularities will cause the alrorighm to fail with a NonconvergenceException.
        /// This is unavoidable for essentially any double-precision numerical integrator. Consider, for example,
        /// the integrable singularity 1/&#x221A;x. Since
        /// &#x3B5; = &#x222B;<sub>0</sub><sup>&#x3B4;</sup> x<sup>-1/2</sup> dx = 2 &#x3B4;<sup>1/2</sup>,
        /// points within &#x3B4; &#x223C; 10<sup>-16</sup> of the end-points, which as a close as you can get to
        /// a point in double precision without being on top of it, contribute at the &#x3B5; &#x223C; 10<sup>-8</sup>
        /// level to our integral, well beyond limit that nearly-full double precision requires. Said differently,
        /// to know the value of the integral to &#x3B5; &#x223C; 10<sup>-16</sup> prescision, we would need to
        /// evaluate the contributions of points within &#x3B4; &#x223C; 10<sup>-32</sup> of the endpoints,
        /// far closer than we can get.</para>
        /// <para>If you need to evaluate an integral with such a strong singularity, make an analytic
        /// change of variable to absorb the singularity before attempting numerical integration. For example,
        /// to evaluate I = &#x222B;<sub>0</sub><sup>b</sup> f(x) x<sup>-1/2</sup> dx, substitute y = x<sup>1/2</sup>
        /// to obtain I = 2 &#x222B;<sub>0</sub><sup>&#x221A;b</sup> f(y<sup>2</sup>) dy.</para>
        /// </remarks>
        public static double Integrate (Func<double, double> integrand, Interval range) {
            EvaluationSettings settings = new EvaluationSettings();
            return (Integrate(integrand, range, settings).Value);
        }


        /// <summary>
        /// Evaluates a definite integral with the given evaluation settings.
        /// </summary>
        /// <param name="integrand">The function to be integrated.</param>
        /// <param name="range">The range of integration.</param>
        /// <param name="settings">The settings which control the evaulation of the integal.</param>
        /// <returns>The result of the integral, which includes an estimated value and estimated uncertainty of that value.</returns>
        public static IntegrationResult Integrate (Func<double,double> integrand, Interval range, EvaluationSettings settings) {

            if (integrand == null) throw new ArgumentNullException("integrand");

            // remap infinite integrals to finite integrals

            if (Double.IsNegativeInfinity(range.LeftEndpoint) && Double.IsPositiveInfinity(range.RightEndpoint)) {

                // -infinity to +infinity

                // remap to (-pi/2,pi/2)
                Func<double, double> f0 = integrand;
                Func<double, double> f1 = delegate (double t) {
                    double x = Math.Tan(t);
                    return (f0(x) * (1.0 + x * x));
                };
                Interval r1 = Interval.FromEndpoints(-Global.HalfPI, Global.HalfPI);

                return (Integrate(f1, r1, settings));

            } else if (Double.IsPositiveInfinity(range.RightEndpoint)) {

                // finite to +infinity

                // remap to interval (-1,1)
                double a0 = range.LeftEndpoint;
                Func<double, double> f0 = integrand;
                Func<double, double> f1 = delegate (double t) {
                    double q = 1.0 - t;
                    double x = a0 + (1 + t) / q;
                    return (f0(x) * (2.0 / q / q));
                };
                Interval r1 = Interval.FromEndpoints(-1.0, 1.0);

                return (Integrate(f1, r1, settings));

            } else if (Double.IsNegativeInfinity(range.LeftEndpoint)) {

                // -infinity to finite

                // remap to interval (-1,1)
                double b0 = range.RightEndpoint;
                Func<double, double> f0 = integrand;
                Func<double, double> f1 = delegate (double t) {
                    double q = t + 1.0;
                    double x = b0 + (t - 1.0) / q;
                    return(f0(x) * (2.0 / q / q));
                };
                Interval r1 = Interval.FromEndpoints(-1.0, 1.0);

                return(Integrate(f1, r1, settings));

            }

            // normal integral over a finitite range

            IAdaptiveIntegrator integrator = new GaussKronrodIntegrator(integrand, range);
            IntegrationResult result = Integrate_Adaptive(integrator, settings);
            return (result);

        }

        // the drivers

        private static IntegrationResult Integrate_Adaptive (IAdaptiveIntegrator integrator, EvaluationSettings s) {

            LinkedList<IAdaptiveIntegrator> list = new LinkedList<IAdaptiveIntegrator>();
            list.AddFirst(integrator);

            int n = integrator.EvaluationCount;

            while (true) {

                // go through the intervals, adding estimates (and errors)
                // and noting which contributes the most error
                // keep track of the total value and uncertainty
                UncertainValue vTotal = new UncertainValue();
                //double v = 0.0;
                //double u = 0.0;

                // keep track of which node contributes the most error
                LinkedListNode<IAdaptiveIntegrator> maxNode = null;
                double maxError = 0.0;

                LinkedListNode<IAdaptiveIntegrator> node = list.First;
                while (node != null) {

                    IAdaptiveIntegrator i = node.Value;

                    UncertainValue v = i.Estimate;
                    vTotal += v;
                    //UncertainValue uv = i.Estimate;
                    //v += uv.Value;
                    //u += uv.Uncertainty;

                    if (v.Uncertainty > maxError) {
                        maxNode = node;
                        maxError = v.Uncertainty;
                    }

                    node = node.Next;

                }

                // if our error is small enough, return
                if ((vTotal.Uncertainty <= Math.Abs(vTotal.Value) * s.RelativePrecision) || (vTotal.Uncertainty <= s.AbsolutePrecision)) {
                    return (new IntegrationResult(vTotal, n, s));
                }
                //if ((vTotal.Uncertainty <= Math.Abs(vTotal.Value) * s.RelativePrecision) || (vTotal.Uncertainty <= s.AbsolutePrecision)) {
                //    return (new IntegrationResult(vTotal.Value, n));
                //}

                // if our evaluation count is too big, throw
                if (n > s.EvaluationBudget) throw new NonconvergenceException();

                // subdivide the interval with the largest error
                IEnumerable<IAdaptiveIntegrator> divisions = maxNode.Value.Divide();
                foreach (IAdaptiveIntegrator division in divisions) {
                    list.AddBefore(maxNode, division);
                    n += division.EvaluationCount;
                    //v2 += division.Estimate;
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

    internal interface IIterativeIntegrator : IIntegrator {

        void Iterate ();

    }

    // integrator base class

    internal abstract class Integrator : IIntegrator {

        protected Integrator (Func<double, double> integrand, Interval range) {
            this.f = integrand;
            this.r = range;
        }

        private Func<double, double> f;
        private Interval r;

        public virtual Func<double, double> Integrand {
            get {
                return (f);
            }
        }

        public virtual Interval Range {
            get {
                return (r);
            }
        }

        public virtual int EvaluationCount {
            get {
                return (n);
            }
        }

        private int n = 0;

        protected double Evaluate (double x) {
            n++;
            return (Integrand(x));
        }

        public abstract UncertainValue Estimate { get; }

    }

    // integrators

    internal class GaussKronrodIntegrator : Integrator, IAdaptiveIntegrator {

        public GaussKronrodIntegrator (Func<double, double> integrand, Interval range) : base(integrand, range) {
            Generation = 1;
            Compute();
        }

        // abcissas and weights 

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


        public virtual IList<IAdaptiveIntegrator> Divide () {
            GaussKronrodIntegrator i1 = new GaussKronrodIntegrator(Integrand, Interval.FromEndpoints(Range.LeftEndpoint, Range.Midpoint));
            i1.Generation = this.Generation + 1;
            GaussKronrodIntegrator i2 = new GaussKronrodIntegrator(Integrand, Interval.FromEndpoints(Range.Midpoint, Range.RightEndpoint));
            i2.Generation = this.Generation + 1;
            return (new GaussKronrodIntegrator[] { i1, i2 });
        }

        public virtual int Generation { get; private set; }
    }

}