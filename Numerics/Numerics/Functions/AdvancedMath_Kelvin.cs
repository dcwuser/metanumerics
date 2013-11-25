using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Meta.Numerics.Functions {

#if FUTURE

    public static partial class AdvancedMath {

        public static double KelvinBei (double nu, double x) {
            throw new NotImplementedException();
        }

        public static double KelvinBer (double nu, double x) {
            if (nu < 0.0) throw new ArgumentOutOfRangeException("nu");
            if (x < 0.0) throw new ArgumentOutOfRangeException("x");

            if (x < 4.0 + 2.0 * Math.Sqrt(nu)) {
                return (BerSeries(nu, x));
            } else {
                throw new NotImplementedException();
            }
        }

        public static double KelvenKei (double nu, double x) {
            throw new NotImplementedException();
        }

        public static double KelvinKer (double nu, double x) {
            throw new NotImplementedException();
        }

        private static double BerSeries (double nu, double x) {

            double c = Math.Cos(3.0 * nu * Math.PI / 4.0);
            double s = Math.Sin(3.0 * nu * Math.PI / 4.0);
            double xh = x / 2.0; double xh2 = xh * xh;

            double df = Math.Pow(xh, nu) / AdvancedMath.Gamma(nu + 1.0);
            double f_old = 0.0;
            double f = c * df;

            for (int k = 1; k < Global.SeriesMax; k++) {
                // we look two values back because for some values of nu (e.g. 0.0) only alternating
                // powers contribute; for these values we would always have f_old == f for non-contributing
                // powers and the series would terminate early
                double f_old_old = f_old; f_old = f;
                df *= xh2 / k / (nu + k);
                switch (k % 4) {
                    case 0:
                        f += c * df;
                        break;
                    case 1:
                        f -= s * df;
                        break;
                    case 2:
                        f -= c * df;
                        break;
                    case 3:
                        f += s * df;
                        break;
                }
                if (f == f_old_old) return (f);
            }
            throw new NonconvergenceException();

        }

    }

#endif

}
