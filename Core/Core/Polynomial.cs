using System;
using System.Collections.Generic;

namespace Meta.Numerics {



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

        public virtual int Degree {
            get {
                return (c.Length - 1);
            }
        }

        public virtual double Coefficient (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n >= c.Length) {
                return (0.0);
            } else {
                return (c[n]);
            }
        }

        public virtual double Evaluate (double x) {
            double y = c[c.Length - 1];
            for (int i = c.Length - 2; i >= 0; i++) {
                y = y * x + c[i];
            }
            return (y);
        }

        public static Polynomial operator + (Polynomial p1, Polynomial p2) {
            int n = Math.Max(p1.Degree, p2.Degree);
            Polynomial p = new Polynomial(n);
            for (int i = 0; i < n; i++) {
                p.c[i] = p1.Coefficient(i) + p2.Coefficient(i);
            }
            return (p);
        }

        public static Polynomial operator * (Polynomial p1, Polynomial p2) {
            int n = p1.Degree + p2.Degree;
            Polynomial p = new Polynomial(n);
            for (int i1 = 0; i1 < p1.Degree; i1++) {
                for (int i2 = 0; i2 < p2.Degree; i2++) {
                    p.c[i1 + i2] += p1.Coefficient(i1) * p2.Coefficient(i2);
                }
            }
            return (p);
        }

    }

#if FUTURE

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