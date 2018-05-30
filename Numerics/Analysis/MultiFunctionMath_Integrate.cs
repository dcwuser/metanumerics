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
        /// <returns>A numerical estimate of the multi-dimensional integral.</returns>
        /// <remarks>
        /// <para>For multi-dimensional integration, the default accurary targets and evaluation budget varies with the dimension of the integral.
        /// The default relative accuracy target falls from about 10<sup>-7</sup> (close to full single precision) for d=2
        /// to about 10<sup>-1</sup> (10%) for d=15. The default absolute accuracy target falls from about 10<sup>-8</sup> to about 10<sup>-2</sup> over the same range.
        /// To achieve those targets, it allows up to about 100,000 evaluations of the integrand for d=2, rising up to about 1 billion evaluations for d=15.</para>
        /// <para>You can change the accuracy targets and evaluation budget by passing an <see cref="EvaluationSettings"/> object to
        /// <see cref="Integrate(Func{IReadOnlyList{double}, double}, IReadOnlyList{Interval}, IntegrationSettings)"/> overload. By decreasing
        /// the accuracy you demand or increasing the evaluation budget, you may be able to successfully complete integrals that would fail for the default settings.</para>
        /// <para>For more information on multi-dimensional integration, see <see cref="Integrate(Func{IReadOnlyList{double}, double}, IReadOnlyList{Interval}, IntegrationSettings)"/>.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="function"/> or <paramref name="volume"/> is <see langword="null"/>.</exception>
        /// <exception cref="NonconvergenceException">The required accuracy could not be achieved using the evaulation budget.</exception>
        public static IntegrationResult Integrate (Func<IReadOnlyList<double>, double> function, IReadOnlyList<Interval> volume) {
            return (Integrate(function, volume, new IntegrationSettings()));
        }

        internal static IntegrationSettings SetMultiIntegrationDefaults (IntegrationSettings original, int d) {
            IntegrationSettings settings = new IntegrationSettings();
            settings.RelativePrecision = (original.RelativePrecision < 0.0) ? Math.Pow(10.0, -(0.0 + 14.0 / d)) : original.RelativePrecision;
            settings.AbsolutePrecision = (original.AbsolutePrecision < 0.0) ? Math.Pow(10.0, -(1.0 + 14.0 / d)) : original.AbsolutePrecision;
            settings.EvaluationBudget = (original.EvaluationBudget < 0.0) ? (int) Math.Round(Math.Pow(10.0, 9.0 - 8.0 / d)) : original.EvaluationBudget;
            settings.Listener = original.Listener;
            return (settings);
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
        /// <see cref="Double.NegativeInfinity"/>. Volumes with dimension greater than 15 are not currently supported.</para>
        /// <para>Integrals with hard boundaries (like our hyper-sphere volume problem) typically require more evaluations than integrals of smooth functions to achieve the same accuracy.
        /// Integrals with canceling positive and negative contributions also typically require more evaluations than integtrals of purely positive functions.</para>
        /// <para>Numerical multi-dimensional integration is computationally expensive. To make problems more tractable, keep in mind some rules of thumb:</para>
        /// <ul>
        ///   <li>Reduce the required accuracy to the minimum required. Alternatively, if you are willing to wait longer, increase the evaluation budget.</li>
        ///   <li>Exploit symmetries of the problem to reduce the integration volume. For example, to compute the volume of the unit d-sphere, it is better to
        ///   integrate over [0,1]<sup>d</sup> and multiply the result by 2<sup>d</sup> than to simply integrate over [-1,1]<sup>d</sup>.</li>
        ///   <li>Apply analytic techniques to reduce the dimension of the integral. For example, when computing the volume of the unit d-sphere, it is better
        ///   to do a (d-1)-dimesional integral over the function that is the height of the sphere in the dth dimension than to do a d-dimensional integral over the
        ///   indicator function that is 1 inside the sphere and 0 outside it.</li>
        /// </ul>
        /// </remarks>
        /// <exception cref="ArgumentException"><paramref name="function"/>, <paramref name="volume"/>, or <paramref name="settings"/> are null, or
        /// the dimension of <paramref name="volume"/> is larger than 15.</exception>
        /// <exception cref="NonconvergenceException">The prescribed accuracy could not be achieved with the given evaluation budget.</exception>
        public static IntegrationResult Integrate (Func<IReadOnlyList<double>, double> function, IReadOnlyList<Interval> volume, IntegrationSettings settings) {

            if (function == null) throw new ArgumentNullException(nameof(function));
            if (volume == null) throw new ArgumentNullException(nameof(volume));
            if (settings == null) throw new ArgumentNullException(nameof(settings));

            // Get the dimension from the box
            int d = volume.Count;

            settings = SetMultiIntegrationDefaults(settings, d);

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
                        box[i] = Interval.FromEndpoints(-1.0, 0.0);
                        map[i] = new TangentCoordinateTransform(limits.RightEndpoint);
                    } else {
                        // a to b
                        box[i] = limits;
                        map[i] = new IdentityCoordinateTransform();
                    }
                }
            }

            // Use adaptive cubature for small dimensions, Monte-Carlo for large dimensions.

            MultiFunctor f = new MultiFunctor(function);
            f.IgnoreInfinity = true;
            f.IgnoreNaN = true;

            UncertainValue estimate;

            if (d < 1) {
                throw new ArgumentException("The dimension of the integration volume must be at least 1.", nameof(volume));
            } else if (d < 4) {
                IntegrationRegion r = new IntegrationRegion(box);
                estimate = Integrate_Adaptive(f, map, r, settings);
            } else if (d < 16) {
                estimate = Integrate_MonteCarlo(f, map, box, settings);
            } else {
                throw new ArgumentException("The dimension of the integrtion volume must be less than 16.", nameof(volume));
            }

            // Sometimes the estimated uncertainty drops precipitiously. We will not report an uncertainty less than 3/4 of that demanded.
            double minUncertainty = Math.Max(0.75 * settings.AbsolutePrecision, Math.Abs(estimate.Value) * 0.75 * settings.RelativePrecision);
            if (estimate.Uncertainty < minUncertainty) estimate = new UncertainValue(estimate.Value, minUncertainty);

            return (new IntegrationResult(estimate, f.EvaluationCount, settings));
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
