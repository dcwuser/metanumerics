using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Meta.Numerics.Statistics {

    
    /// <summary>
    /// Tracks summary statistics for a stream of data points.
    /// </summary>
    public sealed class StreamingSampleSummary {

        public StreamingSampleSummary () { }

        public StreamingSampleSummary (IEnumerable<double> values) : this() {
            Add(values);
        }

        private int n = 0;
        private double m = 0.0;
        private double s2 = 0.0;
        private double s3 = 0.0;
        private double s4 = 0.0;
        private double min = Double.MaxValue;
        private double max = Double.MinValue;

        public void Add (double value) {

            double s2_old = s2;
            double s3_old = s3;

            n++;
            double d = (value - m) / n;
            m += d;
            double e = d * d * n * (n - 1);
            s2 += e;
            s3 += d * (e * (n - 2) - 3.0 * s2_old);
            s4 += d * (d * (e * (n * n - 3 * n + 3) + 6.0 * s2_old) - 4.0 * s3_old);

            if (value < min) min = value;
            if (value > max) max = value;
        }

        public void Add (IEnumerable<double> values) {
            if (values == null) throw new ArgumentNullException(nameof(values));
            foreach (double value in values) Add(value);
        }

        public int Count {
            get {
                return (n);
            }
        }

        public double Mean {
            get {
                return (m);
            }
        }

        public double SecondRawMoment {
            get {
                return (m * m + s2 / n);
            }
        }

        public double ThirdRawMoment {
            get {
                return (s3 / n + m * (3.0 * s2 / n + m * m));
            }
        }

        public double Variance {
            get {
                return (s2 / n);
            }
        }

        public double ThirdCentralMoment {
            get {
                return (s3 / n);
            }
        }

        public double FourthCentralMoment {
            get {
                return (s4 / n);
            }
        }

        public double StandardDeviation {
            get {
                return (Math.Sqrt(Variance));
            }
        }

        public double Minimum {
            get {
                return (min);
            }
        }

        public double Maximum {
            get {
                return (max);
            }
        }

        public UncertainValue PopulationMean {
            get {
                return (new UncertainValue(m, Math.Sqrt(s2 / n / (n - 1))));
            }
        }

        public UncertainValue PopulationVariance {
            get {
                return (new UncertainValue(s2 / (n - 1), Math.Sqrt((s4 / n - MoreMath.Sqr(s2 / n)) / (n - 1))));
            }
        }

        public UncertainValue PopulationStandardDeviation { get; private set; }


        public static StreamingSampleSummary Combine (StreamingSampleSummary a, StreamingSampleSummary b) {

            double d = b.m - a.m;
            double d2 = d * d;
            double d3 = d2 * d;
            double d4 = d2 * d2;

            StreamingSampleSummary ab = new StreamingSampleSummary();
            ab.n = a.n + b.n;
            ab.m = (a.n * a.m + b.n * b.m) / ab.n;
            ab.s2 = a.s2 + b.s2 + d2 * a.n * b.n / ab.n;
            ab.s3 = a.s3 + b.s3 + (d / ab.n) * (d2 * a.n * b.n / ab.n * (a.n - b.n) + 3.0 * (a.n * b.s2 - b.n * a.s2));
            ab.s4 = a.s4 + b.s4 + (d / ab.n) * (d2 * a.n * b.n / ab.n * (d / ab.n) * (a.n * a.n - a.n * b.n + b.n * b.n) + 6.0 * (d / ab.n) * (a.n * a.n * b.s2 + b.n * b.n * a.s2) + 4.0 * (a.n * b.s3 - b.n * a.s3));

            ab.min = Math.Min(a.min, b.min);
            ab.max = Math.Max(a.max, b.max);

            return (ab);
        }


    }


}
