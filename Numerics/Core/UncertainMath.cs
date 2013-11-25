using System;

namespace Meta.Numerics {


    /// <summary>
    /// Contains methods for computing basic mathematical functions of uncertain values.
    /// </summary>
    /// <remarks><para>The methods in this static class perform the same basic mathematical operations as the methods of
    /// the <see cref="System.Math"/> class, accounting for the uncertainty in the inputs to produce a corresponding
    /// uncertainty in the output.</para>
    /// <para>As with operations on uncertain values, the methods assume that the uncertainty in input parameters represents the
    /// standard deviation of a distribution of measurements, and produce a value for the uncertainty in the output which
    /// represent a corresponding standard deviation, under the assumption that the standard deviations are small relative to
    /// the best values.</para></remarks>
    public static class UncertainMath {

        // basic functions of uncertain values

        /// <summary>
        /// Computes the square root of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The square root of the argument.</returns>
        public static UncertainValue Sqrt (UncertainValue x) {
            double v = Math.Sqrt(x.Value);
            double u = 0.5 * x.Uncertainty / v;
            return (new UncertainValue(v, u));
        }

        /// <summary>
        /// Computes the sine of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The sine of the argument.</returns>
        public static UncertainValue Sin (UncertainValue x) {
            UncertainValue y = new UncertainValue();
            y.Value = Math.Sin(x.Value);
            y.Uncertainty = Math.Abs(Math.Cos(x.Value)) * x.Uncertainty;
            return (y);
        }

        /// <summary>
        /// Computes the cosine of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The cosine of the argument.</returns>
        public static UncertainValue Cos (UncertainValue x) {
            double v = Math.Cos(x.Value);
            double u = Math.Abs(Math.Sin(x.Value)) * x.Uncertainty;
            return (new UncertainValue(v, u));
            //UncertainValue y = new UncertainValue();
            //y.Value = Math.Cos(x.Value);
            //y.Uncertainty = Math.Abs(Math.Sin(x.Value)) * x.Uncertainty;
            //return (y);
        }

        /// <summary>
        /// Computes the tangent of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The tanget of the argument.</returns>
        public static UncertainValue Tan (UncertainValue x) {
            double v = Math.Tan(x.Value);
            double u = (1.0 + MoreMath.Pow2(v)) * x.Uncertainty;
            return (new UncertainValue(v, u));
            //UncertainValue y = new UncertainValue();
            //y.Value = Math.Tan(x.Value);
            //y.Uncertainty = (1.0 + Math.Pow(y.Value, 2)) * x.Uncertainty;
            //return (y);
        }

        /// <summary>
        /// Computes the arcsine of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The arcsine of the argument.</returns>
        public static UncertainValue Asin (UncertainValue x) {
            double v = Math.Asin(x.Value);
            double u = x.Uncertainty / Math.Abs(Math.Cos(v));
            return (new UncertainValue(v, u));
            //UncertainValue y = new UncertainValue();
            //y.Value = Math.Asin(x.Value);
            //y.Uncertainty = x.Uncertainty / Math.Abs(Math.Cos(y.Value));
            //return (y);
        }

        /// <summary>
        /// Computes the arccosine of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The arccosine of the argument.</returns>
        public static UncertainValue Acos (UncertainValue x) {
            UncertainValue y = new UncertainValue();
            y.Value = Math.Acos(x.Value);
            y.Uncertainty = x.Uncertainty / Math.Abs(Math.Sin(y.Value));
            return (y);
        }

        /// <summary>
        /// Computes the arctangent of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The arctanget of the argument.</returns>
        public static UncertainValue Atan (UncertainValue x) {
            double xv = x.Value;
            double v = Math.Atan(xv);
            double u = x.Uncertainty / (1 + xv * xv);
            return (new UncertainValue(v, u));
        }

        /// <summary>
        /// Computes the arctangent of the ratio of two uncertain values.
        /// </summary>
        /// <param name="x">The argument of the numerator.</param>
        /// <param name="y">The argument of the denominator.</param>
        /// <returns>The arctangent of the quotient.</returns>
        public static UncertainValue Atan2 (UncertainValue x, UncertainValue y) {
            double xv = x.Value;
            double yv = y.Value;
            double v = Math.Atan2(xv, yv);
            double u1 = x.Uncertainty * yv;
            double u2 = y.Uncertainty * xv;
            double u = Math.Sqrt(u1 * u1 + u2 * u2) / (xv * xv + yv * yv);
            return (new UncertainValue(v, u));
        }

        /// <summary>
        /// Computes e to the power of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of e^<sup>x1</sup>.</returns>
        public static UncertainValue Exp (UncertainValue x) {
            UncertainValue y = new UncertainValue();
            y.Value = Math.Exp(x.Value);
            y.Uncertainty = y.Value * x.Uncertainty;
            return (y);
        }

        /// <summary>
        /// Computes the natural logarithm of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of ln(x1).</returns>
        public static UncertainValue Log (UncertainValue x) {
            UncertainValue y = new UncertainValue();
            y.Value = Math.Log(x.Value);
            y.Uncertainty = x.RelativeUncertainty;
            return (y);
        }

        /// <summary>
        /// Computes an uncertain value raised to an arbitrary power.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <param name="p">The power.</param>
        /// <returns>The argument raised to the specified power.</returns>
        public static UncertainValue Pow (UncertainValue x, double p) {
            UncertainValue y = new UncertainValue();
            y.Value = Math.Pow(x.Value, p);
            y.Uncertainty = Math.Abs(p * y.Value * x.RelativeUncertainty);
            return (y);
        }

        /// <summary>
        /// Computes the hyperbolic sine of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The hyperbolic sine of the argument.</returns>
        public static UncertainValue Sinh (UncertainValue x) {
            UncertainValue y = new UncertainValue();
            y.Value = Math.Sinh(x.Value);
            y.Uncertainty = Math.Cosh(x.Value) * x.Uncertainty;
            return (y);
        }

        /// <summary>
        /// Computes the hyperbolic cosine of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The hyperbolic cosine of the argument.</returns>
        public static UncertainValue Cosh (UncertainValue x) {
            UncertainValue y = new UncertainValue();
            y.Value = Math.Cosh(x.Value);
            y.Uncertainty = Math.Abs(Math.Sinh(x.Value)) * x.Uncertainty;
            return (y);
        }

        /// <summary>
        /// Computes the hyperbolic tangent of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The hyperbolic tanget of the argument.</returns>
        public static UncertainValue Tanh (UncertainValue x) {
            UncertainValue y = new UncertainValue();
            y.Value = Math.Tanh(x.Value);
            y.Uncertainty = y.Uncertainty / Math.Pow(Math.Cosh(x.Value), 2);
            return (y);
        }

    }

}