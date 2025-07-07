using Meta.Numerics.Extended;
using System;

namespace Meta.Numerics.Analysis {

    // This class is used to wrap a function, storing some associated state such as the evaluation count.
    // It isn't truly a functor in the C++ sense, since .NET doesn't allow () to be overloaded, but
    // it is a functor in the sense that it is a class used to represent a function.

    internal enum SpecialValueBehavior { Throw, Zero };

    internal class Functor {

        public Functor(Func<double, double> f)  {
            this.Function = f;
        }

        public Func<double, double> Function { get; protected set; }

        public int EvaluationCount { get; protected set; }

        public int EvaluationBudget { get; set; } = Int32.MaxValue;

        public SpecialValueBehavior NonFiniteBehavior { get; set; }

        public virtual double Evaluate(double x) {

            EvaluationCount++;
            if (EvaluationCount >= EvaluationBudget) {
                throw new NonconvergenceException();
            }

            double y = Function(x);

            if (ExtendedMath.IsNotFinite(y)) {
                switch (NonFiniteBehavior) {
                    case SpecialValueBehavior.Zero:
                        y = 0.0;
                        break;
                    default:
                        throw new InvalidOperationException();
                }
            }

            return y;

        }

    }
}
