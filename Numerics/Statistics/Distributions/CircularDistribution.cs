using System;

using Meta.Numerics;
using Meta.Numerics.Analysis;

namespace Meta.Numerics.Statistics.Distributions {

    public abstract class CircularDistribution {

        public abstract Interval Support { get;  }

        public abstract double ProbabilityDensity (double x);

        public virtual double ExpectationValue (Func<double, double> function) {
            if (function == null) throw new ArgumentNullException(nameof(function));
            IntegrationResult result = FunctionMath.Integrate(x => function(x) * this.ProbabilityDensity(x), this.Support);
            return (result.Value);
        }

        private double MapToTheta (double x) {
            return ((x - this.Support.LeftEndpoint) + Global.TwoPI * (this.Support.RightEndpoint - x));
        }

        public virtual Complex Moment (int r) {
            if (r == 0) {
                return (1.0);
            } else {
                double re = ExpectationValue(x => Math.Cos(r * MapToTheta(x)));
                double im = ExpectationValue(x => Math.Sign(r * MapToTheta(x)));
                return new Complex(re, im);
            }
        }

        public virtual double Variance {
            get {
                return (1.0 - ComplexMath.Abs(Moment(1)));
            }
        }
    }

    public class CircularUniformDistribution : CircularDistribution {

        public CircularUniformDistribution () : base () {
            this.support = Interval.FromEndpoints(0.0, Global.TwoPI);
        }

        public CircularUniformDistribution (Interval support) : base () {
            this.support = support;
        }

        private readonly Interval support;

        public override Interval Support {
            get {
                return support;
            }
        }

        public override double ProbabilityDensity (double x) {
            return (1.0 / support.Width);
        }

        public override Complex Moment (int r) {
            if (r == 0) {
                return (Complex.One);
            } else {
                return (Complex.Zero);
            }
        }

        public override double Variance {
            get {
                return (1.0);
            }
        }
    }
}
