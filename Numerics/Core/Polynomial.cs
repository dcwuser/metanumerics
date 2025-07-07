using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Text;

using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics {

#if FUTURE

    internal class HermiteHePolynomial : Polynomial {

        int n;

        public override int Degree {
            get {
                return n;
            }
        }

        public override double Evaluate(double x) {
            return OrthogonalPolynomials.HermiteHe(n, x);
        }

        public override double Coefficient(int n) {
            throw new NotImplementedException();
        }

    }

#endif

    internal class CoefficientPolynomial : Polynomial {

        internal CoefficientPolynomial(double[] coefficients) {
            Debug.Assert(coefficients is object);
            this.coefficients = coefficients;
            degree = -1;
            for (int i = coefficients.Length - 1; i >= 0; i--) {
                if (coefficients[i] != 0.0) {
                    degree = i;
                    break;
                }
            }
        }

        private readonly double[] coefficients;
        private readonly int degree;

        public override int Degree {
            get {
                return degree;
            }
        }

        public override double Coefficient(int n) {
            if (n < 0 || n > degree) {
                return 0.0;
            } else {
                Debug.Assert(0 <= n && n < coefficients.Length);
                return coefficients[n];
            }
        }

        public override double Evaluate(double x) {

            // This is Horner's scheme, which is well-known but non-obvious, at least to me.
            //   a_0 + a_1 x + a_2 x^2 + a_3 x^3 = ((a_3 x + a_2) x + a_1) x + a_0
            // It requires 2d flops, i.e. only one multiply and one add for each degree added,
            // so it is more efficient than naïve evaluation.
            // It is also ostensibly less subject to rounding errors than naïve evaluation.

            double t = 0.0;
            for (int i = degree; i >= 0; i--) {
                t = t * x + coefficients[i];
            }
            return t;
        }

    }


    /// <summary>
    /// Represents a polynomial with real coefficients.
    /// </summary>
    public abstract class Polynomial {

        // Make constructor internal so others cannot inherit
        internal Polynomial () { }

        /// <summary>
        /// Initializes a new polynomial with the given coefficients.
        /// </summary>
        /// <param name="coefficients">The coefficients of the polynomial.</param>
        /// <returns>The specified polynomial.</returns>
        /// <remarks>
        /// <para>Coefficients should be arranged from low to high order, so that the kth entry is the coefficient of x<sup>k</sup>. For example,
        /// to specify the polynomial 5 - 6 x + 7 x<sup>3</sup>, give the values 5, -6, 0, 7.</para>
        /// </remarks>
        public static Polynomial FromCoefficients(params double[] coefficients) {
            if (coefficients is null) throw new ArgumentNullException(nameof(coefficients));
            if (coefficients.Length == 0) throw new InvalidOperationException();
            return new CoefficientPolynomial(coefficients);
        }


        /// <summary>
        /// Initializes a new polynomial that passes through the given points.
        /// </summary>
        /// <param name="points">An N X 2 array whose first column contains the x values of points and whose second column contains the corresponding y values.</param>
        /// <returns>A polynomial of degree N-1 that passes through all the given points.</returns>
        public static Polynomial FromPoints(double[,] points) {
            if (points == null) throw new ArgumentNullException(nameof(points));
            if (points.GetLength(0) == 0) throw new ArgumentException("The first dimension of the points array must have length at least one.", nameof(points));
            if (points.GetLength(1) != 2) throw new ArgumentException("The second dimension of the points array must have length two.", nameof(points));
            double[] x = new double[points.GetLength(0)];
            double[] y = new double[points.GetLength(0)];
            for (int i = 0; i < points.GetLength(0); i++) {
                x[i] = points[i, 0];
                y[i] = points[i, 1];
            }
            PolynomialInterpolator interpolator = new PolynomialInterpolator(x, y);
            return (new InterpolatingPolynomial(interpolator));
        }

        /// <summary>
        /// Initializes a new polynomial that passes through the given points.
        /// </summary>
        /// <param name="points">A collection of points.</param>
        /// <returns>A polynomial that passes through all the given points.</returns>
        public static Polynomial FromPoints(ICollection<XY> points) {
            if (points == null) throw new ArgumentNullException(nameof(points));
            if (points.Count == 0) throw new ArgumentException("There must be at least one point in the points collection.", nameof(points));
            double[] x = new double[points.Count];
            double[] y = new double[points.Count];
            int i = 0;
            foreach (XY point in points) {
                x[i] = point.X;
                y[i] = point.Y;
                i++;
            }
            PolynomialInterpolator interpolator = new PolynomialInterpolator(x, y);
            return (new InterpolatingPolynomial(interpolator));
        }

        /// <summary>
        /// Gets the degree of the polynomial.
        /// </summary>
        /// <remarks>
        /// <para>The degree of a polynomial is the highest power of the variable that appears. For example, the degree of 5 + 6 x + 7 x<sup>2</sup> is 2.</para>
        /// </remarks>
        public abstract int Degree { get; }

        /// <summary>
        /// Gets the specified coefficient.
        /// </summary>
        /// <param name="n">The power of the variable for which the coefficient is desired.</param>
        /// <returns>The coefficient of x<sup>n</sup>.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative.</exception>
        public abstract double Coefficient(int n);

        /// <summary>
        /// Evaluates the polynomial for the given input value.
        /// </summary>
        /// <param name="x">The value of the variable.</param>
        /// <returns>The value of the polynomial.</returns>
        public abstract double Evaluate(double x);

        /// <summary>
        /// Differentiates the polynomial.
        /// </summary>
        /// <returns>The derivative of the polynomial.</returns>
        public virtual Polynomial Differentiate() {
            double[] newCoefficients = new double[this.Degree];
            for (int i = 0; i < newCoefficients.Length; i++) newCoefficients[i] = (i + 1) * this.Coefficient(i + 1);
            return new CoefficientPolynomial(newCoefficients);
        }

        /// <summary>
        /// Integrates the polynomial.
        /// </summary>
        /// <param name="C">The integration constant.</param>
        /// <returns>The integral of the polynomial.</returns>
        public virtual Polynomial Integrate(double C) {
            double[] newCoefficients = new double[this.Degree + 2];
            newCoefficients[0] = C;
            for (int i = 1; i < newCoefficients.Length; i++) newCoefficients[i] = this.Coefficient(i - 1) / i;
            return new CoefficientPolynomial(newCoefficients);
        }

#if FUTURE

        /// <summary>
        /// Computes the roots of the polynomial.
        /// </summary>
        /// <returns>An array containing complex values for which the polynomial vanishes.</returns>
        public virtual Complex[] FindRoots () {

            // We should really implement Box-Jenkins, but since we already have a well-developed eigenvalue finder,
            // just form the companion matrix and find its eigenvalues. 
            // https://en.wikipedia.org/wiki/Companion_matrix
            SquareMatrix C = new SquareMatrix(this.Degree);
            double q = this.Coefficient(this.Degree);
            Debug.Assert(q != 0.0);
            for (int r = 1; r < C.Dimension; r++) C[r, r - 1] = 1.0;
            for (int r = 0; r < C.Dimension; r++) C[r, C.Dimension - 1] = -this.Coefficient(r) / q;
            Complex[] roots = C.Eigenvalues();
            return roots;

        }

        /// <summary>
        /// Computes a composed polynomial.
        /// </summary>
        /// <param name="q">The polynomial argument.</param>
        /// <returns>The polynomial p(q(x)).</returns>
        public virtual Polynomial ComposeWith (Polynomial q) {
            throw new NotImplementedException();
        }

#endif

        /// <summary>
        /// Generates a string representation of the polynomial.
        /// </summary>
        /// <returns>A string representation of the polynomial.</returns>
        public override string ToString() {
            return ToString(CultureInfo.CurrentCulture);
        }

        private string ToString(IFormatProvider provider) {
            NumberFormatInfo numberFormat = (NumberFormatInfo) provider.GetFormat(typeof(NumberFormatInfo));
            StringBuilder text = new StringBuilder();
            bool leading = true;
            for (int i = 0; i <= this.Degree; i++) {
                double coefficient = Coefficient(i);
                // Skip zero terms
                if (coefficient == 0.0) continue;
                if (leading) {
                    // Leading terms get formatted a little differently since there is no + or - in front of them.
                    leading = false;
                } else {
                    // For non-leading coefficients, use binary + or - instead signed coefficient.
                    if (coefficient < 0.0) {
                        text.Append($" {numberFormat.NegativeSign} ");
                        coefficient = -coefficient;
                    } else {
                        text.Append($" {numberFormat.PositiveSign} ");
                    }       
                }
                if (i == 0) {
                    // For terms with no power of x, write the coefficient in all cases.
                    // Note that if this term is non-zero, then it is the leading term.
                    text.Append(coefficient);
                } else {
                    // For terms that get a power of x, writing a 1 is unnecessary.
                    if (coefficient != 1.0) text.AppendFormat(provider,"{0} ", coefficient);
                    if (i == 1) {
                        text.Append('x');
                    } else {
                        text.AppendFormat(provider, "x^{0}", i);
                    }
                }
            }
            return text.ToString();
        }

        /// <summary>
        /// Negates a polynomial.
        /// </summary>
        /// <param name="p">The polynomial.</param>
        /// <returns>The additive inverse of the polynomial.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="p"/> is <see langword="null"/>.</exception>
        public static Polynomial operator - (Polynomial p) {
            if (p == null) throw new ArgumentNullException(nameof(p));
            double[] coefficients = new double[p.Degree + 1];
            for (int i = 0; i < coefficients.Length; i++) {
                coefficients[i] = -p.Coefficient(i);
            }
            return new CoefficientPolynomial(coefficients);
        }

        /// <summary>
        /// Computes the sum of two polynomials.
        /// </summary>
        /// <param name="p1">The first polynomial.</param>
        /// <param name="p2">The second polynomial.</param>
        /// <returns>The sum polynomial.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="p1"/> or <paramref name="p2"/> is <see langword="null"/>.</exception>
        public static Polynomial operator + (Polynomial p1, Polynomial p2) {
            if (p1 == null) throw new ArgumentNullException(nameof(p1));
            if (p2 == null) throw new ArgumentNullException(nameof(p2));
            double[] coefficients = new double[Math.Max(p1.Degree, p2.Degree) + 1];
            for (int i = 0; i < coefficients.Length; i++) coefficients[i] = p1.Coefficient(i) + p2.Coefficient(i);
            return new CoefficientPolynomial(coefficients);
        }

        /// <summary>
        /// Computes the difference of two polynomials.
        /// </summary>
        /// <param name="p1">The first polynomial.</param>
        /// <param name="p2">The second polynomial.</param>
        /// <returns>The difference polynomial.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="p1"/> or <paramref name="p2"/> is <see langword="null"/>.</exception>
        public static Polynomial operator - (Polynomial p1, Polynomial p2) {
            if (p1 == null) throw new ArgumentNullException(nameof(p1));
            if (p2 == null) throw new ArgumentNullException(nameof(p2));
            double[] coefficients = new double[Math.Max(p1.Degree, p2.Degree) + 1];
            for (int i = 0; i < coefficients.Length; i++) coefficients[i] = p1.Coefficient(i) - p2.Coefficient(i);
            return new CoefficientPolynomial(coefficients);
        }

        /// <summary>
        /// Computes the product of two polynomials.
        /// </summary>
        /// <param name="p1">The first polynomial.</param>
        /// <param name="p2">The second polynomial.</param>
        /// <returns>The product polynomial.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="p1"/> or <paramref name="p2"/> is null.</exception>
        public static Polynomial operator * (Polynomial p1, Polynomial p2) {
            if (p1 == null) throw new ArgumentNullException(nameof(p1));
            if (p2 == null) throw new ArgumentNullException(nameof(p2));
            double[] coefficients = new double[p1.Degree + p2.Degree + 1];
            for (int i1 = 0; i1 <= p1.Degree; i1++) {
                for (int i2 = 0; i2 <= p2.Degree; i2++) {
                    coefficients[i1 + i2] += p1.Coefficient(i1) * p2.Coefficient(i2);
                }
            }
            return new CoefficientPolynomial(coefficients);
        }

        /// <summary>
        /// Computes the quotient of two polynomials.
        /// </summary>
        /// <param name="p1">The dividend polynomial.</param>
        /// <param name="p2">The divisor polynomial.</param>
        /// <param name="remainder">The remainder polynomial.</param>
        /// <returns>The quotient polynomial.</returns>
        /// <remarks>
        /// <para>p<sub>1</sub> = q p<sub>2</sub> + r</para>
        /// </remarks>
        public static Polynomial Divide (Polynomial p1, Polynomial p2, out Polynomial remainder) {

            if (p1 == null) throw new ArgumentNullException(nameof(p1));
            if (p2 == null) throw new ArgumentNullException(nameof(p2));

            if (p2.Degree >= p1.Degree) {
                throw new InvalidOperationException();
            }

            // compute and store some value we will use repeatedly
            int d1 = p1.Degree;
            int d2 = p2.Degree;
            double a = p2.Coefficient(d2);
            // a is the leading coefficient of p2; it is the only number we ever divide by
            // (so if p2 is monic, polynomial division does not involve division at all!)

            // copy the coefficients of p1 into an working array; this will become our remainder when we are done
            double[] r = new double[d1 + 1];
            for (int i = 0; i < r.Length; i++) r[i] = p1.Coefficient(i);

            // create space for the coefficients of the quotient polynomial, which has order p1.Order - p2.Order
            double[] q = new double[d1 - d2 + 1];

            // do long division
            for (int k = q.Length - 1; k >= 0; k--) {
                q[k] = r[d2 + k] / a;
                r[d2 + k] = 0.0;
                for (int j = d2 + k - 1; j >= k; j--) {
                    r[j] = r[j] - q[k] * p2.Coefficient(j - k);
                }
            }

            // form the remainder and quotient polynomials from the arrays
            remainder = new CoefficientPolynomial(r);
            return new CoefficientPolynomial(q);

        }

    }

}