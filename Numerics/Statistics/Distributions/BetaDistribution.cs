﻿using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a beta distribution.
    /// </summary>
    /// <remarks>
    /// <para>The beta distribution is defined on the interval [0,1]. Depending on its two shape parameters, it can take on a
    /// variety of forms on this interval.</para>
    /// <para>The left shape parameter &#x3b1; controls the shape of the distribution near the left endpoint x = 0.
    /// The right shapre paramater &#x3b2; controls the shape of the distribution near the right endpoint x = 1.
    /// If a shape parameter is less than one, the distribution is singular on that side. If a shape parameter is greater
    /// than one, the distribution goes to zero on that side. If a shape parameter is equal to one, the distribution
    /// goes to a constant on that side.</para>
    /// <para>If the two shape parameters are equal, the distribution is symmetric.</para>
    /// <para>When both shape parameters are one, the beta distribution reduces to a standard uniform distribution.</para>
    /// <para><img src="../images/UniformFromBeta.png" /></para>
    /// <para>Beta distributions describe the maximum and minimum values obtained from multiple, independent draws from a standard
    /// uniform distribution. For n draws, the maximum value is distributed as B(n,1).</para>
    /// <para><img src="../images/MaxUniformBetaRelation.png" /></para>
    /// <para>Similarly, the minimum value is distributed as B(1,n).</para>
    /// <para><img src="../images/MinUniformBetaRelation.png" /></para>
    /// <para>In fact, the <i>i</i>th order statistic (<i>i</i>th smallest value) in n draws from a uniform distribution
    /// is distributed as B(i, n - i + 1).</para>
    /// <para>Because of the wide variety of shapes it can take, the beta distribution is sometimes
    /// used as an ad hoc model to fit distributions observed on a finite interval.</para>
    /// </remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/Beta_distribution"/>
    /// <seealso href="http://mathworld.wolfram.com/BetaDistribution.html"/>
    public sealed class BetaDistribution : ContinuousDistribution {

        /// <summary>
        /// Initializes a new &#x3B2; distribution.
        /// </summary>
        /// <param name="alpha">The left shape parameter, which controls the form of the distribution near x=0.</param>
        /// <param name="beta">The right shape parameter, which controls the form of the distribution near x=1.</param>
        /// <remarks>
        /// <para>The <paramref name="alpha"/> shape parameter controls the form of the distribution near x=0. The
        /// <paramref name="beta"/> shape parameter controls the form of the distribution near z=1. If a shape parameter
        /// is less than one, the PDF diverges on the side of the distribution it controls. If a shape parameter
        /// is greater than one, the PDF goes to zero on the side of the distribution it controls. If the left and right
        /// shape parameters are equal, the distribution is symmetric about x=1/2.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="alpha"/> or <paramref name="beta"/> is non-positive.</exception>
        public BetaDistribution (double alpha, double beta) {
            if (alpha <= 0.0) throw new ArgumentOutOfRangeException(nameof(alpha));
            if (beta <= 0.0) throw new ArgumentOutOfRangeException(nameof(beta));
            this.a = alpha;
            this.b = beta;
            // cache value of B(alpha, beta) to avoid having to re-calculate it whenever needed
            this.bigB = AdvancedMath.Beta(alpha, beta);
            // get a beta generator
            if (alpha < 0.75 && beta < 0.75) {
                this.betaRng = new JoehnkBetaGenerator(alpha, beta);
            } else if (alpha > 1.0 && beta > 1.0) {
                this.betaRng = new ChengBetaGenerator(alpha, beta);
            } else {
                this.betaRng = DeviateGeneratorFactory.GetBetaGenerator(alpha, beta);
            }
            // get a beta inverter
            this.betaInverter = new BetaInverter(alpha, beta);
        }

        private readonly double a, b;
        private readonly double bigB;
        private readonly IDeviateGenerator<double> betaRng;
        private readonly BetaInverter betaInverter;

        /// <summary>
        /// Gets the left shape parameter.
        /// </summary>
        public double Alpha {
            get {
                return a;
            }
        }

        /// <summary>
        /// Gets the right shape parameter.
        /// </summary>
        public double Beta {
            get {
                return b;
            }
        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return Interval.Unit;
            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if ((x < 0.0) || (x > 1.0)) {
                return 0.0;
            } else {
                return Math.Pow(x, a - 1.0) * Math.Pow(1.0 - x, b - 1.0) / bigB;
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= 0.0) {
                return 0.0;
            } else if (x >= 1.0) {
                return 1.0;
            } else {
                // the value is I_x(a,b), which is LeftRegularizedBeta
                // since we have precomputed bigB, use it here instead of calling
                // that method, which would re-compute it
                return AdvancedMath.LeftRegularizedBeta(a, b, x);
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x <= 0.0) {
                return 1.0;
            } else if (x >= 1.0) {
                return 0.0;
            } else {
                // use 1.0 - I_x(a, b) = I_{1-x}(b, a), which is essentially the symmetry of
                // the beta distribution definition, to avoid calculating 1 - small number
                return AdvancedMath.LeftRegularizedBeta(b, a, 1.0 - x);
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return a / (a + b);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                double ab = a + b;
                return a * b / (ab + 1.0) / (ab * ab);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                double ab = a + b;
                return 2.0 * (b - a) / (ab + 2.0) * Math.Sqrt((ab + 1.0) / (a * b));
            }
        }

        /// <inheritdoc />
        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else {
                // this is just a recursive development of \Beta(\alpha + r, \beta) / \Beta(\alpha, \beta)
                double ab = a + b;
                double M = 1.0;
                for (int i = 0; i < r; i++) {
                    M = (a + i) / (ab + i) * M;
                }
                return M;
            }
        }

        // According to http://mathworld.wolfram.com/BetaDistribution.html, central moments are given by
        //   C_n = \left( -\frac{\alpha}{\alpha+\beta} \right)^n 2F1(\alpha, -n; \alpha + \beta; \frac{\alpha+\beta}{\alpha})
        // where 2F1 is a hypergeometric function. The Gauss recurrence for hypergeometric functions (Abromowitz and Stegun 15.2.22)
        //   [c - 2 b + (b - a) z] 2F1(a, b; c; z) + b (1 - z) 2F1(a, b + 1; c; z) - (c - b) 2F1(a, b - 1; c; z) = 0
        // implies
        //   C_{n+1} = \frac{n}{(\alpha + \beta + n) (\alpha + \beta)} \left[ (\beta - \alpha) C_{n} + \frac{\alpha \beta}{\alpha + \beta} C_{n - 1} \right]
        // This recurrence appears to work and not suffer from the cancelation errors that computation from the raw moments does.

        /// <inheritdoc />
        public override double CentralMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return 1.0;
            } else if (r == 1) {
                return 0.0;
            } else {

                // use recursion

                double ab = a + b;
                double ba = b - a;
                double u = a * b / ab;

                double CM = 1.0;
                double C0 = 0.0;
                for (int k = 1; k < r; k++) {
                    double CP = k / (ab + k) / ab * (ba * C0 + u * CM);
                    CM = C0;
                    C0 = CP;
                }
                return C0;

            }
        }

        internal override double[] CentralMoments (int rMax) {
            double[] central = new double[rMax + 1];
            IEnumerator<double> iterator = CentralMoments();
            for (int k = 0; k < central.Length; k++) {
                iterator.MoveNext();
                central[k] = iterator.Current;
            }
            return (central);
        }

        private IEnumerator<double> CentralMoments () {

            yield return (1.0);
            yield return (0.0);

            double ab = a + b;
            double ba = b - a;
            double u = a * b / ab;

            double CM = 1.0;
            double C0 = 0.0;
            for (int i = 1; true; i++) {
                double CP = i / (ab + i) / ab * (ba * C0 + u * CM);
                CM = C0;
                C0 = CP;
                yield return (C0);
            }

        }

        // inverse CDF

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));

            betaInverter.InverseRegularizedBeta(P, 1.0 - P, out double x, out double y);
            return x;
        }

        /// <inheritdoc />
        public override double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException(nameof(Q));

            betaInverter.InverseRegularizedBeta(1.0 - Q, Q, out double x, out double y);
            return x;
        }

        
        /// <inheritdoc />
        public override double GetRandomValue (Random rng) {
            if (rng is null) throw new ArgumentNullException(nameof(rng));
            return betaRng.GetNext(rng);
        }
        

        /// <summary>
        /// Computes the Beta distribution that best fits the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The best fit parameters.</returns>
        /// <remarks>
        /// <para>The returned fit parameters are the &#x3B1; (<see cref="Alpha"/>) and  &#x3B2; (<see cref="Beta"/>) parameters, in that order.
        /// These are the same parameters, in the same order, that are required by the <see cref="BetaDistribution(double,double)"/> constructor to
        /// specify a new Beta distribution.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than three values.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Not all the entries in <paramref name="sample" /> lie between zero and one.</exception>
        public static BetaFitResult FitToSample (IReadOnlyList<double> sample) {
            return Univariate.FitToBeta(sample);
        }

    }

    // The inversion of the Beta CDF turns out to be quite complex. We have moved the logic into an auxiliary class for tidiness. 

    internal class BetaInverter {

        public BetaInverter (double a, double b)
            : this(a, b, AdvancedMath.LogBeta(a, b)) {
        }

        public BetaInverter (double a, double b, double logB) {
            if (a <= 0.0) throw new ArgumentOutOfRangeException(nameof(a));
            if (b <= 0.0) throw new ArgumentOutOfRangeException(nameof(b));
            this.a = a;
            this.b = b;
            this.a1 = a - 1.0;
            this.b1 = b - 1.0;
            this.ma = a / (a + b);
            this.mb = b / (a + b);
            this.logB = logB;
        }

        private readonly double a, b;

        private readonly double a1, b1;

        private readonly double ma, mb;

        private readonly double logB;

        // We implement the inverse series with explicitly passed a and b so that we can reverse a and b
        // when expanding from the the right instead of the left

        private static double InverseRegularizedBetaSeriesArgument (double a, double logB, double P) {
            double w = Math.Exp((Math.Log(a * P) + logB) / a);
            return (w);
        }

        private static double InverseRegularizedBetaSeries (double a, double b, double w, out bool exact) {
            double t1 = (b - 1.0) / (a + 1.0) * w;
            double f1 = 1.0 + t1;
            double f2 = f1 + t1 * (a * a + 3.0 * b * a - a + 5.0 * b - 4.0) / (a + 1.0) / (a + 2.0) / 2.0 * w;
            exact = (f2 == f1);
            return (w * f2);
        }

        // For large a, b, I_x(a, b) approaches a normal distribution with m = a/(a+b) and s = \sqrt{\frac{ab}{a+b+1}} / (a + b).
        // We try to correct for the skewness and excess kuritosis using the Cornish-Fisher expansion.

        private bool NormalInverseRegularizedBeta (double P, double Q, out double x, out double y) {

            // Normal approximation
            double z = AdvancedMath.Probit(P, Q);

            // Cornish-Fisher correction
            double c3 = 2.0 * (b - a) / (a + b + 2.0) * Math.Sqrt((a + b + 1.0) / a / b);
            z = z + (z * z - 1.0) * c3 / 6.0;

            // Convert to x/y
            double s = Math.Sqrt(ma * mb / (a + b + 1.0));
            double d = s * z;

            x = ma + d;
            y = mb - d;

            // Accept only if not too close to boundaries
            return ((x > s / 4.0) && (y > s / 4.0));

        }

        private void CrudeInverseRegularizedBeta (double P, double Q, out double x, out double y) {
            double Ba = Math.Pow(ma, a) / a;
            double Bb = Math.Pow(mb, b) / b;
            double B = Ba + Bb;
            double I0 = Ba / B;
            if (P < I0) {
                x = Math.Pow(a * B * P, 1.0 / a);
                y = 1.0 - x;
            } else {
                y = Math.Pow(b * B * Q, 1.0 / b);
                x = 1.0 - y;
            }


        }

        public bool ApproximateRegularizedInverseBeta (double P, double Q, out double x, out double y) {

            Debug.Assert(P + Q == 1.0);

            // Check for trivial cases

            if (P == 0.0) {
                x = 0.0;
                y = 1.0;
                return (true);
            }

            if (Q == 0.0) {
                y = 0.0;
                x = 1.0;
                return (true);
            }

            // Try to use the inverse series from each side. This is most likely to be useful if the relevant shape parameter a
            // is less than one, so the power (1/a) is greater than one, making the expansion parameter small.

            // Since the last term gives an error estimate, we need not further refine the result if
            // the series result was already accurate to double precision.

            bool exact = false;
            if (a <= 2.0) {
                double w = InverseRegularizedBetaSeriesArgument(a, logB, P);
                if ((w < 0.25) && (b * w < 0.25)) {
                    x = InverseRegularizedBetaSeries(a, b, w, out exact);
                    y = 1.0 - x;
                    return (exact);
                }
            }

            if (b <= 2.0) {
                double w = InverseRegularizedBetaSeriesArgument(b, logB, Q);
                if ((w < 0.25) && (a * w < 0.25)) {
                    y = InverseRegularizedBetaSeries(b, a, w, out exact);
                    x = 1.0 - y;
                    return (exact);
                }
            }

            // If the parameters are large enough that the shape is vaguely normal, try
            // a normal approximation. But reject if the result is too close to the edge
            // of the allowed range.

            if ((a >= 2.0) && (b >= 2.0)) {
                bool wellInBounds = NormalInverseRegularizedBeta(P, Q, out x, out y);
                if (wellInBounds) {
                    return (false);
                }
            }

            // If no other approximations were available, use the "crude approximation" described in NR.
            // This approximation usually isn't too bad (within a few %) and also has the advantage
            // of being guaranteed never to go out of bounds.

            CrudeInverseRegularizedBeta(P, Q, out x, out y);
            return (false);

        }

        private void IncompleteBetaAndDerivative (double x, double y, out double P, out double Q, out double D) {
            D = Math.Exp(a1 * Math.Log(x) + b1 * MoreMath.LogOnePlus(-x) - logB);
            double xtp = (a + 1.0) / (a + b + 2.0);
            if (x < xtp) {
                P = D * x * y * AdvancedMath.RegularizedBeta_ContinuedFraction(a, b, x);
                Q = 1.0 - P;
            } else {
                Q = D * x * y * AdvancedMath.RegularizedBeta_ContinuedFraction(b, a, y);
                P = 1.0 - Q;
            }
        }

        // Ideally we would treat x and y symmetrically but that turns out to be quite hard.
        // Since in current applications we always return x, for now we are content to
        // ensure x is accurate for small x but not give the same guarantee for small y.

        private void RefineInverseRegularizedBeta (double P0, double Q0, ref double x, ref double y) {

            Debug.Assert((0 < P0) && (P0 <= 1.0));
            Debug.Assert((0 < Q0) && (Q0 <= 1.0));
            Debug.Assert(P0 + Q0 == 1.0);

            // maintain bounds to ensure Newton does not shoot out of range
            double xMin = 0.0;
            double xMax = 1.0;

            for (int k = 0; k < 32; k++) {

                // Remember old coordinates
                double x_old = x;
                //double y_old = y;

                // Evaluate I and I' at the coordinates
                IncompleteBetaAndDerivative(x, y, out double P, out double Q, out double D);
                Debug.Assert(P + Q == 1.0);

                // Compute the deficit from the smaller of P and Q, to reduce inaccuracies when P or Q is near 1
                double z = (P0 < Q0) ? P - P0 : Q0 - Q;

                // Adjust the bounds based on the deficit
                if (z < 0.0) {
                    // we are to the left of the sought point
                    xMin = x;
                } else if (z > 0.0) {
                    // we are to the right of the sought point
                    xMax = x;
                } else {
                    // we found it, so return
                    return;
                }

                // Compute the Halley step, or just the Newton step if the Halley correction dominates
                double DC = z / 2.0 * (a1 / x - b1 / y);
                double dx = (Math.Abs(DC) < Math.Abs(D) / 2.0) ? -z / (D - DC) : -z / D;
                x += dx;

                //Debug.WriteLine("k={0} x={1:R} y={2:R} [x0={3:R} x1={4:R}] P={5} Q={6} z={7} D={8} DC={9} dx={10}", k, x, y, xMin, xMax, P, Q, z, D, DC, dx);
                
                // If the new point is out of bounds, use bisection instead.
                if (x <= xMin) {
                    dx = (xMin - x_old) / 2.0;
                    x = x_old + dx;
                } else if (x >= xMax) {
                    dx = (xMax - x_old) / 2.0;
                    x = x_old + dx;
                }

                y = 1.0 - x;

                // If the change is small enough, we're finished.
                // Change x -> Min(x, y) to treat x and y symmetrically, but that causes problems with bounds enforcement.
                if (Math.Abs(dx) <= x * Global.Accuracy) return;

            }

            throw new NonconvergenceException();
        }

        public void InverseRegularizedBeta (double P, double Q, out double x, out double y) {

            // Get an initial approximation
            bool exact = ApproximateRegularizedInverseBeta(P, Q, out x, out y);

            // If the initial approximation was already full precision, just return it 
            if (exact) return;

            // Otherwise, refine it via Newton-Halley iteration
            RefineInverseRegularizedBeta(P, Q, ref x, ref y);

        }


    }

    // http://math.ubbcluj.ro/~tradu/Randvargen.pdf

    // R. C. H. Cheng, Generating Beta Variates with Nonintegral Shape Parameters,
    // Communications of the ACM 21 (1978) 317

    internal class ChengBetaGenerator : IDeviateGenerator<double> {

        public ChengBetaGenerator (double a, double b) {
            Debug.Assert(a > 1.0);
            Debug.Assert(b > 1.0);
            this.a = a;
            this.b = b;
            this.alpha = a + b;
            this.beta = Math.Sqrt((alpha - 2.0) / (2.0 * a * b - alpha));
            this.gamma = a + 1.0 / beta;
        }

        // Cheng uses the names alpha and beta for some values he calculates. They are
        // not the shape parametes of the distribution, which he calls a and b.
        private readonly double a, b;
        private readonly double alpha, beta, gamma;

        private static readonly double log4 = Math.Log(4.0);

        public double GetNext (Random rng) {
            while (true) {
                double u = rng.NextDouble();
                double v = beta * Math.Log(u / (1.0 - u));
                double w = a * Math.Exp(v);
                double z = u * u * rng.NextDouble();
                if (gamma * v - log4 + alpha * Math.Log(alpha / (b + w)) > Math.Log(z)) {
                    return w / (b + w);
                }
            }
        }
    }

    // Johnk M.D., Erzeugung von Betaverteilten und Gammaverteilten Zufallszahlen; Metrika 8 (1964) 5-15

    internal class JoehnkBetaGenerator : IDeviateGenerator<double> {

        public JoehnkBetaGenerator (double a, double b) {
            Debug.Assert(a < 1.0);
            Debug.Assert(b < 1.0);
            this.ai = 1.0 / a;
            this.bi = 1.0 / b;
        }

        private readonly double ai, bi;

        public double GetNext (Random rng) {
            while (true) {
                double u = rng.NextDouble();
                double v = rng.NextDouble();
                double x = Math.Pow(u, ai);
                double y = Math.Pow(v, bi);
                double z = x + y;
                if (z < 1.0) return x / z;
            }

        }
    }


}