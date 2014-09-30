using System;
using System.Collections.Generic;
using System.Threading;
using System.Collections.ObjectModel;

using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Meta.Numerics.Analysis {

    public static partial class MultiFunctionMath {

        /// <summary>Estimates a multi-dimensional integral.</summary>
        /// <param name="function">The function to integrate.</param>
        /// <param name="volume">The volume over which to integrate.</param>
        /// <remarks>
        /// <para>By default, our multidimensional integration system targets a relative accuracy of about 10<sup>-7</sup> (close to full single precision) for d=2, falling gradually
        /// to about 10<sup>-2</sup> (1%) for d=12. To achieve that accuracy, it allows up to about 10<sup>5</sup> evaluations of the integrand for d=2, rising
        /// up to about 10<sup>8</sup> evaluations for d=12.</para>
        /// <para>You can change the accuracy demands and evaluation budget by passing an <see cref="EvaluationSettings"/> object to the integration method. By decreasing
        /// the accuracy you require or increasing the evaluation budget, you may be able to successfully complete integrals that would fail for the default settings.</para>
        /// </remarks>
        public static IntegrationResult Integrate (Func<IList<double>, double> function, IList<Interval> volume) {

            if (function == null) throw new ArgumentNullException("function");
            if (volume == null) throw new ArgumentNullException("volume");

            // Compute a target accuracy that decreases with dimension from full at d = 1 to ~1/8 at d = 8.
            //   d = 2, eps = 3E-8, max = 65K
            //   d = 3, eps = 1E-5, max = 330K
            //   d = 4, eps = 2E-4, max = 1M
            //   d = 6, eps =, max = 
            //   d = 8, eps =, max =
            //   d = 10, eps =, max =

            // Choose default evaluation settings based on dimension.
            int d = volume.Count;
            double precision = Math.Pow(10.0, -(1.0 + 12.0 / d)); // change to 10^{-(1+12/d)} or 2^{-3(1 + 12/d)}?
            int count = (d * d * d * d) * 4096;

            EvaluationSettings settings = new EvaluationSettings() {
                AbsolutePrecision = precision * precision,
                RelativePrecision = precision,
                EvaluationBudget = d * d * d * d * 4096
            };

            return (Integrate(function, volume, settings));

        }

        /// <summary>
        /// Estimates a multi-dimensional integral using the given evaluation settings.
        /// </summary>
        /// <param name="function">The function to be integrated.</param>
        /// <param name="volume">The volume over which to integrate.</param>
        /// <param name="settings">The integration settings.</param>
        /// <returns>A numerical estimate of the multi-dimensional integral.</returns>
        /// <remarks>
        /// <para>Note that the integration function must not attempt to modify the argument passed to it.</para>
        /// <para>Note that the integration volume must be a hyper-rectangle. You can integrate over regions with more complex boundaries by specifying the integration
        /// volume as a bounding hyper-rectangle that encloses your desired integration region, and returing the value 0 for the integrand outside of the desired integration
        /// region. For example, to find the volume of a unit d-sphere, you can integrate a function that is 1 inside the unit d-sphere and 0 outside it over the volume
        /// [-1,1]<sup>d</sup>. You can integrate over infinite volumes by specifing volume endpoints of <see cref="Double.PositiveInfinity"/> and/or
        /// <see cref="Double.NegativeInfinity"/>. Volumes with dimension greater than 12 are not currently supported.</para>
        /// <para>Integrals with hard boundaries (like our hyper-sphere volume problem) typically require more evaluations than integrals of smooth functions to achieve the same accuracy.
        /// Integrals with canceling positive and negative contributions also typically require more evaluations than integtrals of purely positive functions.</para>
        /// <para>Numerical multi-dimensional integration is computationally expensive. To make problems more tractable, keep in mind some rules of thumb:</para>
        /// <ul>
        ///   <li>Reduce the required accuracy to the minimum required. Alternatively, if you are willing to wait longer, increase the evaluation budget.</li>
        ///   <li>Exploit symmetries of the problem to reduce the integration volume. For example, to compute the volume of the unit d-sphere, it is better to
        ///   integrate of [0,1]<sup>d</sup> and multiply the result by 2<sup>d</sup> than to simply integrate over [-1,1]<sup>d</sup>.</li>
        ///   <li>Apply analytic techniques to reduce the dimension of the integral. For example, when computing the volume of the unit d-sphere, it is better
        ///   to do a (d-1)-dimesional integral over the function that is the height of the sphere in the dth dimension than to do a d-dimensional integral over the function
        ///   that is 1 inside the sphere and 0 outside it.</li>
        /// </ul>
        /// </remarks>
        /// <exception cref="ArgumentException"><paramref name="function"/>, <paramref name="volume"/>, or <paramref name="settings"/> are null, or
        /// the dimension of <paramref name="volume"/> is larger than 12.</exception>
        /// <exception cref="NonconvergenceException">The prescribed accuracy could not be achieved with the given evaluation budget.</exception>
        public static IntegrationResult Integrate (Func<IList<double>, double> function, IList<Interval> volume, EvaluationSettings settings) {

            if (function == null) throw new ArgumentNullException("function");
            if (volume == null) throw new ArgumentNullException("volume");
            if (settings == null) throw new ArgumentNullException("settings");

            // Get the dimension from the box
            int d = volume.Count;

            // Translate the integration volume, which may be infinite, into a bounding box plus a coordinate transform.
            Interval[] box = new Interval[d];
            CoordinateTransform[] map = new CoordinateTransform[d];
            for (int i = 0; i < d; i++) {
                Interval limits = volume[i];
                if (Double.IsInfinity(limits.RightEndpoint)) {
                    if (Double.IsInfinity(limits.LeftEndpoint)) {
                        // -\infinity to +\infinity
                        box[i] = Interval.FromEndpoints(-1.0, +1.0);
                        map[i] = new TangentCoordinateTransform(0.0);
                    } else {
                        // a to +\infinity
                        box[i] = Interval.FromEndpoints(0.0, 1.0);
                        map[i] = new TangentCoordinateTransform(limits.LeftEndpoint);
                    }
                } else {
                    if (Double.IsInfinity(limits.LeftEndpoint)) {
                        // -\infinity to b
                        throw new NotImplementedException();
                    } else {
                        // a to b
                        box[i] = limits;
                        map[i] = new IdentityCoordinateTransform();
                    }
                }
            }

            // Use adaptive cubature for small dimensions, Monte-Carlo for large dimensions.

            MultiFunctor f;
            UncertainValue estimate;

            if (d < 1) {
                throw new InvalidOperationException();
            } else if (d < 4) {
                f = new MultiFunctor(function);
                IntegrationRegion r = new IntegrationRegion(box);
                estimate = Integrate_Adaptave(f, map, r, settings);
            } else if (d < 12) {
                f = new MultiFunctor(function);
                estimate = Integrate_MonteCarlo(f, map, box, settings);
            } else {
                throw new ArgumentException();
            }

            // Sometimes the estimated uncertainty drops precipitiously.
            double minUncertainty = Math.Max(0.75 * settings.AbsolutePrecision, Math.Abs(estimate.Value) * 0.75 * settings.RelativePrecision);
            if (estimate.Uncertainty < minUncertainty) estimate = new UncertainValue(estimate.Value, minUncertainty);

            return (new IntegrationResult(estimate, f.EvaluationCount, settings));
        }

    }


    internal class MultiFunctor {

        private readonly Func<IList<double>, double> function;

        private readonly CoordinateTransform[] map;

        private readonly bool negate = false;

        private int count = 0;

        public MultiFunctor (Func<IList<double>, double> function, bool negate) : this(function) {
            this.negate = negate;
        }

        public MultiFunctor (Func<IList<double>, double> function) : this(function, null) { }

        public MultiFunctor (Func<IList<double>, double> function, CoordinateTransform[] map) {
            this.function = function;
            this.map = map;
        }

        public double Evaluate (double[] x) {
            Interlocked.Increment(ref count);
            if (map == null) {
                double z = function(new ReadOnlyCollection<double>(x));
                //if (Double.IsInfinity(z) || Double.IsNaN(z)) z = 0.0;
                if (negate) z = -z;
                return (z);
            } else {
                double j = 1.0;
                for (int i = 0; i < x.Length; i++) {
                    CoordinateTransform transform = map[i];
                    if (transform != null) transform.TransformInPlace(ref x[i], ref j);
                }
                double z = function(new ReadOnlyCollection<double>(x));
                if (Double.IsInfinity(z) || Double.IsNaN(z)) z = 0.0;
                return (j * z);
            }

        }

        public int EvaluationCount {
            get {
                return (count);
            }
        }

        public bool IsNegated {
            get {
                return (negate);
            }
        }

    }


    internal abstract class CoordinateTransform {

        public virtual void TransformInPlace (ref double x, ref double j) {
            double j1;
            x = Transform(x, out j1);
            j *= j1;
        }

        public abstract double Transform (double x, out double j);

    }

    internal class IdentityCoordinateTransform : CoordinateTransform {

        public override void TransformInPlace (ref double x, ref double j) {
            // no change
        }

        public override double Transform (double x, out double j) {
            j = 1.0;
            return (x);
        }
    }


    internal class TangentCoordinateTransform : CoordinateTransform {

        public TangentCoordinateTransform (double a) {
            this.a = a;
        }

        private double a;

        public override void TransformInPlace (ref double x, ref double j) {
            double u = Math.PI / 2.0 * x;
            x = a + Math.Tan(u);
            j *= Math.PI / 2.0 / MoreMath.Sqr(Math.Cos(u));

        }

        public override double Transform (double x, out double j) {
            double u = Math.PI / 2.0 * x;
            j = Math.PI / 2.0 / MoreMath.Sqr(Math.Cos(u));
            return (a + Math.Tan(u));
        }
    }

}
