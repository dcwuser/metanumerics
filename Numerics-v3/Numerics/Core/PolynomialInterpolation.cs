using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Meta.Numerics {

    internal class InterpolatingPolynomial : Polynomial {

        internal InterpolatingPolynomial (PolynomialInterpolator interpolator) : base(interpolator.GetCoefficients()) {
            this.interpolator = interpolator;
        }

        private PolynomialInterpolator interpolator;

        public override double Evaluate (double x) {
            return (interpolator.Evaluate(x));
        }

    }


    internal class PolynomialInterpolator {

        public PolynomialInterpolator (int capacity) {
            if (capacity < 1) throw new ArgumentOutOfRangeException("capacity");
            this.n = 0;
            this.x = new double[capacity];
            this.y = new double[capacity][];
            for (int i = 0; i < capacity; i++) {
                this.y[i] = new double[capacity - i];
            }
        }

        internal PolynomialInterpolator (double[] x, double[] y) {
            if (x == null) throw new ArgumentNullException("x");
            if (y == null) throw new ArgumentNullException("y");
            if (x.Length != y.Length) throw new InvalidOperationException();

            n = x.Length;

            this.x = x;
            this.y = new double[n][];
            this.y[0] = y;
            for (int i = 1; i < n; i++) {
                this.y[i] = new double[n - i];
            }

        }

        private double[] x;
        private double[][] y;
        int n;

        public void Add (double x, double y) {
            if (n >= this.x.Length) throw new InvalidOperationException();
            this.x[n] = x;
            this.y[0][n] = y;
            n++;
        }

        public int Order {
            get {
                return (n - 1);
            }
        }

        public double Evaluate (double xm) {
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < (n - i); j++) {
                    y[i][j] = ((xm - x[i + j]) * y[i - 1][j] + (x[j] - xm) * y[i - 1][j + 1]) / (x[j] - x[i + j]);
                }
            }
            return (y[n - 1][0]);
        }

        public double[] GetCoefficients () {

            // The algorithm for a Vandermonde system is adapted from Golub & Van Loan, "Matrix Computations", p. 185

            // Copy the values into scratch space. At the end of the algorithm this space will hold the coefficients.
            double[] c = new double[n];
            Array.Copy(y[0], c, n);

            int m = n - 1; // degree of polynomial

            for (int i = 0; i < m; i++) {
                for (int j = m; j >= i + 1; j--) {
                    c[j] = (c[j] - c[j - 1]) / (x[j] - x[j - i - 1]);
                }
            }
            for (int k = m - 1; k >= 0; k--) {
                for (int i = k; i < m; i++) {
                    c[i] = c[i] - c[i + 1] * x[k];
                }
            }

            return (c);


        }

    }

}


