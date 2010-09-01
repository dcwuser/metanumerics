using System;
using System.Collections.Generic;

namespace Meta.Numerics {

#if FUTURE

    public class Polynomial {

        public Polynomial (int n) {
            if (n < 0) throw new ArgumentOutOfRangeException("n");
            c = new double[n + 1];
        }

        public Polynomial (IList<double> coeffcients) {
            if (coeffcients == null) throw new ArgumentNullException("coefficients");
            c = new double[coeffcients.Count];
            for (int i = 0; i < coeffcients.Count; i++) {
                c[i] = coeffcients[i];
            }
        }

        private double[] c;

        public int Degree {
            get {
                return (c.Length - 1);
            }
        }

        public double Evaluate (double x) {
            throw new NotImplementedException();
        }

    }

    public static class PolynomialMath {

        public static Polynomial Add (Polynomial p1, Polynomial p2) {
            throw new NotImplementedException();
        }

        public static Polynomial Subtract (Polynomial p1, Polynomial p2) {
            throw new NotImplementedException();
        }

        public static Polynomial Multiply (Polynomial p1, Polynomial p2) {
            throw new NotImplementedException();
        }

        public static Polynomial Divide (Polynomial p1, Polynomial p2, out Polynomial remainder) {
            throw new NotImplementedException();
        }

        public static Polynomial Differentiate (Polynomial p) {
            throw new NotImplementedException();
        }

        public static Polynomial Integrate (Polynomial p) {
            throw new NotImplementedException();
        }

        public static Complex[] FindRoots (Polynomial p) {
            throw new NotImplementedException();
        }

    }

#endif

}