using System;
using System.Diagnostics;
using System.Collections.Generic;


namespace Meta.Numerics.Analysis {

    // Strategies drive engines.

    // It's been a long road to get to this object model, and I'm still not convinced it's really
    // the best design, but it's better than what I started with.

    // One problem is that a lot of the details of each strategy depend on things like whether
    // there are multiple variables and whether the RHS is y' or y''. So a naive implementation
    // repeats a lot of strategy logic multiple times, one for each variation. The strategy/engine
    // split factors out all the variation-specific logic into the engine part, so the strategy
    // part appears only once.

    internal abstract class OdeStrategy {

        protected OdeStrategy (IOdeEngine engine) {
            Debug.Assert(engine != null);
            this.engine = engine;
        }

        private readonly IOdeEngine engine;

        public void IntegrateTo (double x1) {

            double x0 = engine.X;

            // Reverse direction, if necessary.
            if (Math.Sign(engine.DeltaX) != Math.Sign(x1 - x0)) engine.DeltaX = -engine.DeltaX;

            // We can't just check (X < X1) because sometimes we integrate the other way,
            // so instead check that "we are on the same side of X1 as X0".
            while (Math.Sign(engine.X - x1) == Math.Sign(x0 - x1)) {

                // If we would overshoot in the next step, reduce it.
                if (Math.Sign(engine.X + engine.DeltaX - x1) != Math.Sign(x0 - x1)) engine.DeltaX = x1 - engine.X;

                // Try to take a step.
                // This is the key method to be implemented by any stepper.
                // After it is called, if the step is sucessful, it should advance X (presumably by DeltaX but smaller is okay),
                // and update the value of Y to correspond to the new X. If the step is unsuccessful, it should leave X and Y
                // unchanged. If the stepper is adaptive, the method may also change DeltaX in preperation for the next step.

                Step();

            }

        }

        public abstract void Step ();

    }

    internal interface IOdeEngine {

        double DeltaX { get; set; }

        double X { get; }

        int SuccessCount { get; }

        int FailureCount { get; }

        void CompleteStep (bool accept);

    }

    internal class OdeEngine : IOdeEngine {

        public OdeEngine (double x) {
            this.x = x;
        }

        private double x;

        public double X {
            get {
                return (x);
            }
        }

        public double DeltaX {
            get {
                return (deltaX);
            }
            set {
                deltaX = value;
            }
        }

        public double deltaX;

        private int successCount = 0;

        private int failureCount = 0;

        public int SuccessCount {
            get {
                return (successCount);
            }
        }

        public int FailureCount {
            get {
                return (failureCount);
            }
        }

        protected virtual void AcceptStep () {
            x += DeltaX;
            successCount++;
        }

        protected virtual void RejectStep () {
            failureCount++;
        }

        public virtual void CompleteStep (bool accept) {
            if (accept) {
                AcceptStep();
            } else {
                RejectStep();
            }
        }

    }

}
