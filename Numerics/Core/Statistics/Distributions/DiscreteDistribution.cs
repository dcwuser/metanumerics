using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics {

    // really, the inheritance hierarchy should be:
    //   Distribution
    //     ContinuousDistribution
    //       Concrete continuous distributions
    //     DiscreteDistribution
    //       Concrete distrete distributions
    // but what should have been ContinuousDistribution was defined as Distribution
    // we should correct this in a future new major release

    /// <summary>
    /// Represents all discrete distrubtions.
    /// </summary>
    /// <remarks>
    /// <para>A discrete distribution is a distribution over the integers.</para>
    /// </remarks>
    public abstract class DiscreteDistribution {

        // replace ProbabilityDensity with ProbabilityMass

        /// <summary>
        /// Returns the probability of the obtaining the given value.
        /// </summary>
        /// <param name="k">The value.</param>
        /// <returns>The probability of obtaining the value.</returns>
        public abstract double ProbabilityMass (int k);

        /// <summary>
        /// Gets the interval over which the distribution is nonvanishing.
        /// </summary>
        public abstract DiscreteInterval Support { get; }

        /// <summary>
        /// Computes the probability of obtaining a value less than or equal to the given value.
        /// </summary>
        /// <param name="k">The value.</param>
        /// <returns>The total probability of obtaining a value less than or equal to <paramref name="k"/>.</returns>
        /// <seealso cref="RightProbability"/>
        public virtual double LeftProbability (int k) {
            DiscreteInterval support = Support;
            if (k < support.LeftEndpoint) {
                return (0.0);
            } if (k >= support.RightEndpoint) {
                return (1.0);
            } else {
                double P = 0.0;
                for (int j = support.LeftEndpoint; j <= k; j++) {
                    P += ProbabilityMass(j);
                }
                return (P);
            }
        }

        /// <summary>
        /// Computes the probability of obtaining a value greater than the given value.
        /// </summary>
        /// <param name="k">The value.</param>
        /// <returns>The total probability of obtaining a value greater than <paramref name="k"/>.</returns>
        /// <seealso cref="LeftProbability"/>
        public virtual double RightProbability (int k) {
            DiscreteInterval support = Support;
            if (k < support.LeftEndpoint) {
                return (1.0);
            } else if (k >= support.RightEndpoint) {
                return (0.0);
            } else {
                double Q = 0.0;
                for (int j = k + 1; j <= support.RightEndpoint; j++) {
                    Q += ProbabilityMass(j);
                }
                return (Q);
            }
        }

        public virtual int InverseLeftProbability (double P) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Computes the expectation value of an artibrary function.
        /// </summary>
        /// <param name="f">The function.</param>
        /// <returns>The expectation value of the function.</returns>
        public virtual double ExpectationValue (Function<int, double> f) {
            if (f == null) throw new ArgumentNullException("f");
            DiscreteInterval support = Support;

            int i0 = (int) Math.Round(Mean);
            Console.WriteLine("i0={0}", i0);
            double s_left = 0.0;
            for (int i = i0 - 1; i >= support.LeftEndpoint; i--) {
                double s_left_old = s_left;
                double ds = f(i) * ProbabilityMass(i);
                s_left += f(i) * ProbabilityMass(i);
                Console.WriteLine("i = {5}, ds = f * p = {0} * {1}, s1 = s0 + ds = {2} {3} = {4}", f(i), ProbabilityMass(i), s_left_old, ds, s_left, i);
                if (s_left == s_left_old) break;
            }
            double s_right = 0.0;
            for (int i = i0 + 1; i <= support.RightEndpoint; i++) {
                double s_right_old = s_right;
                double ds = f(i) * ProbabilityMass(i);
                s_right += f(i) * ProbabilityMass(i);
                Console.WriteLine("i = {5}, ds = f * p = {0} * {1}, s1 = s0 + ds = {2} {3} = {4}", f(i), ProbabilityMass(i), s_right_old, ds, s_right, i);
                if (s_right == s_right_old) break;
            }
            return (s_left + s_right + f(i0) * ProbabilityMass(i0));
            /*
            double s = 0.0;
            for (int j = support.LeftEndpoint; j <= support.RightEndpoint; j++) {
                s += f(j) * ProbabilityMass(j);
            }
            return (s);
            */
        }

        /// <summary>
        /// Gets the mean of the distribution.
        /// </summary>
        public virtual double Mean {
            get {
                Function<int, double> f = delegate (int k) {
                    return (k);
                };
                return (ExpectationValue(f));
            }

        }

        /// <summary>
        /// Gets the variance of the distribution.
        /// </summary>
        public virtual double Variance {
            get {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Gets the skewness of the distribution.
        /// </summary>
        public virtual double Skewness {
            get {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Getst the standard deviation of the distribution.
        /// </summary>
        public virtual double StandardDeviation {
            get {
                return (Math.Sqrt(Variance));
            }
        }

        /// <summary>
        /// Gets a raw moment of the distribution.
        /// </summary>
        /// <param name="n">The order of the moment.</param>
        /// <returns>The raw moment M<sub>n</sub>.</returns>
        public virtual double Moment (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else {
                // if m is large, this is slow; find a better way
                return (ExpectationValue(delegate(int k) { return (MoreMath.Pow(k, n)); }));
            } 
        }

        /// <summary>
        /// Gets a central moment of the distribution.
        /// </summary>
        /// <param name="n">The order of the moment.</param>
        /// <returns>The central moment C<sub>n</sub>.</returns>
        public virtual double MomentAboutMean (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (0.0);
            } else {
                double mu = Mean;
                return (ExpectationValue(delegate(int k) { return (MoreMath.Pow(k - mu, n)); }));
            }
        }

    }

    /// <summary>
    /// Represents a discrete distribution as a continous distribution.
    /// </summary>
    internal class DiscreteAsContinuousDistribution : Distribution {

        public DiscreteAsContinuousDistribution (DiscreteDistribution distribution) {
            if (distribution == null) throw new ArgumentNullException("distribution");
            this.d = distribution;
        }

        private DiscreteDistribution d;

        /// <summary>
        /// Not implemented for discrete distributions.
        /// </summary>
        /// <param name="x">The </param>
        /// <returns></returns>
        /// <remarks>
        /// <para>The probability density of a discrete distribution consists of Dirac delta functions
        /// at the integers, each with an integrated area equal to the probability mass of the discrete
        /// distribution at that integer. Because such a function cannot be represented as a
        /// numeric function over floating point types, this method is not implemented.</para>
        /// </remarks>
        public override double ProbabilityDensity (double x) {
            throw new NotImplementedException();
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return d.Mean;
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return d.Variance;
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (d.Skewness);
            }
        }

        /// <inheritdoc />
        public override double Moment (int n) {
            return(d.Moment(n));
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            return (d.MomentAboutMean(n));
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            return (d.LeftProbability((int)Math.Truncate(x)));
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            return (d.RightProbability((int)Math.Truncate(x)));
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            return(d.InverseLeftProbability(P));
        }

    }

    /// <summary>
    /// Represents an interval on the integers.
    /// </summary>
    public struct DiscreteInterval {

        private DiscreteInterval (int a, int b) {
            if (b < a) Global.Swap(ref a, ref b);
            this.a = a;
            this.b = b;
        }

        int a, b;

        /// <summary>
        /// Gets the lower (left) boundary.
        /// </summary>
        public int LeftEndpoint {
            get {
                return (a);
            }
        }

        /// <summary>
        /// Gets the upper (right) boundary.
        /// </summary>
        public int RightEndpoint {
            get {
                return (b);
            }
        }

        /// <summary>
        /// Gets the width of the interval.
        /// </summary>
        public int Width {
            get {
                return (b - a);
            }
        }

        /// <summary>
        /// Tests whether the given value is contained in the open interval.
        /// </summary>
        /// <param name="x">The value to test.</param>
        /// <returns>True if the value lies in the open interval, otherwise false.</returns>
        public bool OpenContains (double x) {
            return ((IntToDouble(a) < x) && (x < IntToDouble(b)));
        }

       /// <summary>
        /// Tests whether the given value is contained in the closed interval.
        /// </summary>
        /// <param name="x">The value to test.</param>
        /// <returns>True if the value lies in the open interval, otherwise false.</returns>
        public bool ClosedContains (double x) {
            return ((IntToDouble(a) <= x) && (x <= IntToDouble(b)));
        }

        /// <summary>
        /// Creates a new discrete interval, given the boundary integers.
        /// </summary>
        /// <param name="a">One boundary.</param>
        /// <param name="b">The other boundary.</param>
        /// <returns>A discrete interval structure representing the interval.</returns>
        public static DiscreteInterval FromEndpoints (int a, int b) {
            return new DiscreteInterval(a, b);
        }

        // convert int to double, replacing minimum and maxiumum integer values by negative and positivel infinity

        private static double IntToDouble (int i) {
            if (i == Int32.MinValue) {
                return(Double.NegativeInfinity);
            } else if (i == Int32.MaxValue) {
                return(Double.PositiveInfinity);
            } else {
                return((double) i);
            }
        }

        public static implicit operator Interval (DiscreteInterval i) {
            return (Interval.FromEndpoints(IntToDouble(i.a), IntToDouble(i.b)));
        }

    }


}