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

        public abstract DiscreteInterval Support { get; }

        // 

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
        public virtual double ExpectationValue (Function<double, double> f) {
            if (f == null) throw new ArgumentNullException("f");
            DiscreteInterval support = Support;
            double s = 0.0;
            for (int j = support.LeftEndpoint; j <= support.RightEndpoint; j++) {
                s += f(j) * ProbabilityMass(j);
            }
            return (s);
        }

        /// <summary>
        /// Gets the mean of the distribution.
        /// </summary>
        public virtual double Mean {
            get {
                Function<double, double> f = delegate (double k) {
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
            throw new NotImplementedException();
        }

        /// <summary>
        /// Gets a central moment of the distribution.
        /// </summary>
        /// <param name="n">The order of the moment.</param>
        /// <returns>The central moment C<sub>n</sub>.</returns>
        public virtual double MomentAboutMean (int n) {
            throw new NotImplementedException();
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
    /// Represented a Poisson distribution.
    /// </summary>
    public class PoissonDistribution : DiscreteDistribution {

        private double mu;

        public override double ProbabilityMass (int k) {
            if (k < 0) {
                return (0.0);
            } else {
                return (Math.Exp(
                    k * Math.Log(mu) - AdvancedIntegerMath.LogFactorial(k) - mu
                ));
            }
        }

        public override DiscreteInterval Support {
            get {
                return (DiscreteInterval.FromEndpoints(0, Int32.MaxValue));
            }
        }

        public override double Mean {
            get {
                return (mu);
            }
        }

        public override double Variance {
            get {
                return (mu);
            }
        }

    }

    

    public struct Interval<T> where T : struct, IComparable<T> {

        private T a, b;

        public Interval (T a, T b) {
            if (a.CompareTo(b) <= 0) {
                this.a = a;
                this.b = b;
            } else {
                this.a = b;
                this.b = a;
            }
        }

        public T LeftEndpoint {
            get {
                return(a);
            }
        }

        public T RightEndpoint {
            get {
                return(b);
            }
        }

        public bool OpenContains (T x) {
            return( (a.CompareTo(x) < 0) && (b.CompareTo(x) > 0) );
        }

        public bool ClosedContains (T x) {
            return( (a.CompareTo(x) <= 0) && (b.CompareTo(x) >= 0) );
        }

    }

    public struct DiscreteInterval {

        private DiscreteInterval (int a, int b) {
            if (b < a) Global.Swap(ref a, ref b);
            this.a = a;
            this.b = b;
        }

        int a, b;

        public int LeftEndpoint {
            get {
                return (a);
            }
        }

        public int RightEndpoint {
            get {
                return (b);
            }
        }

        public int Width {
            get {
                return (b - a);
            }
        }

        public static DiscreteInterval FromEndpoints (int a, int b) {
            return new DiscreteInterval(a, b);
        }

    }


}