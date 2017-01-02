using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a triangular distribution.
    /// </summary>
    /// <remarks>
    /// <para>Like a uniform distribution, a triangular distribution is confined to a finite interval. Unlike a
    /// uniform distribution, a triangular distribution is not uniform across the interval.</para>
    /// <para>Triangular distributions are often used in project planning, where a maximum, minimum, and most
    /// likely value for some quantity is known or supposed.</para>
    /// </remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/Triangular_distribution" />
    public sealed class TriangularDistribution : Distribution {

        /*
        public TriangularDistribution (Interval range) {
            r = range;
            c = range.Midpoint;
        }

        private TriangularDistribution (Interval range, double peak) {
            if (!range.ClosedContains(peak)) throw new ArgumentOutOfRangeException("peak");
            r = range;
            c = peak;
        }
        */

        /// <summary>
        /// Initializes a new triangular distribution.
        /// </summary>
        /// <param name="a">One inflection point of the distribution.</param>
        /// <param name="b">A second inflection point of the distribution.</param>
        /// <param name="c">A third inflection point of the distribution.</param>
        public TriangularDistribution (double a, double b, double c) {

            // check for validity
            if ((a == b) && (b == c)) throw new InvalidOperationException();

            // order points
            if (a > c) Global.Swap(ref a, ref c);
            if (a > b) Global.Swap(ref a, ref b);
            if (b > c) Global.Swap(ref b, ref c);

            // record values
            this.a = a;
            this.b = b;
            this.c = c;

            this.ab = b - a;
            this.bc = c - b;
            this.ac = c - a;

            this.h = 2.0 / ac;

        }

        private double a, b, c;
        private double ab, bc, ac;
        double h;

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if ((x <= a) || (x >= c)) {
                return (0.0);
            } else {
                if (x < b) {
                    return (h * (x - a) / ab);
                } else {
                    return (h * (c - x) / bc);
                }
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= a) {
                return (0.0);
            } else if (x >= c) {
                return (1.0);
            } else {
                if (x < b) {
                    return (LeftTriangleArea(x));
                } else {
                    return (1.0 - RightTrigangleArea(x));
                }
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x <= a) {
                return (1.0);
            } else if (x >= c) {
                return (0.0);
            } else {
                if (x < b) {
                    return (1.0 - LeftTriangleArea(x));
                } else {
                    return (RightTrigangleArea(x));
                }
            }
        }

        private double LeftTriangleArea (double x) {
            Debug.Assert(x >= a); Debug.Assert(x <= b);
            double ax = x - a;
            return (ax * ax / ab / ac);
        }

        private double RightTrigangleArea (double x) {
            Debug.Assert(x >= b); Debug.Assert(x <= c);
            double xc = c - x;
            return (xc * xc / bc / ac);
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            double Pb = ab / ac;
            if (P < Pb) {
                return (a + Math.Sqrt(ab * ac * P));
            } else {
                return (c - Math.Sqrt(bc * ac * (1.0 - P)));
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return ((a + b + c) / 3.0);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return ((ab * ab + bc * bc + ac * ac) / 36.0);
            }
        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(a, c));
                //return (r);
            }
        }

        private double MomentAboutMode (int n) {
            if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return ((bc - ab) / 3.0);
            } else {
                // ac = ab + bc is always a factor of the numerator; find a way to divide it out analytically
                if (n % 2 == 0) {
                    return (2.0 * (Math.Pow(bc, n + 1) + Math.Pow(ab, n + 1)) / (n + 1) / (n + 2) / ac);
                } else {
                    return (2.0 * (Math.Pow(bc, n + 1) - Math.Pow(ab, n + 1)) / (n + 1) / (n + 2) / ac);
                }
            }
        }

        /// <inheritdoc />
        public override double Moment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException("r");
            } else if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (Mean);
            } else {
                double M = MomentAboutMode(r);
                double t = 1.0;
                for (int k = r - 1; k >= 0; k--) {
                    t *= b;
                    M += AdvancedIntegerMath.BinomialCoefficient(r, k) * MomentAboutMode(k) * t;
                }
                return (M);
            }
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException("r");
            } else if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (0.0);
            } else {
                double M = MomentAboutMode(r);
                double s = -MomentAboutMode(1);
                double t = 1.0;
                for (int k = r - 1; k >= 0; k--) {
                    t *= s;
                    M += AdvancedIntegerMath.BinomialCoefficient(r, k) * MomentAboutMode(k) * t;
                }
                //Console.WriteLine("n={0}, M={1}", n, M);
                return (M);
            }
        }
    }

}