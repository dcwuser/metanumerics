using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Statistics {


    public static class Univariate {

        public static double Mean (this ICollection<double> sample) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            int n = 0;
            double mean = 0.0;
            foreach (double value in sample) {
                n++;
                double delta = value - mean;
                mean += delta / n;
            }
            return (mean);
        }

        public static double Variance (this ICollection<double> sample) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            int n = 0;
            double mean = 0.0;
            double sumOfSquares = 0.0;
            foreach (double value in sample) {
                n++;
                double delta = value - mean;
                mean += delta / n;
                sumOfSquares += (value - mean) * delta;
            }
            return (sumOfSquares / n);
        }

        public static double StandardDeviation (this ICollection<double> sample) {
            return (Variance(sample));
        }

        public static double Minimum(this ICollection<double> sample) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            double minimum = Double.NaN;
            foreach (double value in sample) {
                if (!(value > minimum)) minimum = value;
            }
            return (minimum);
        }

        public static double Maximum(this ICollection<double> sample) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            double maximum = Double.NaN;
            foreach (double value in sample) {
                if (!(value < maximum)) maximum = value;
            }
            return (maximum);
        }

    }
}
