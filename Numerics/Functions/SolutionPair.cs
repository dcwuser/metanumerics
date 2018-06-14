using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Meta.Numerics.Functions {

    /// <summary>
    /// Contains a pair of solutions to a differential equation.
    /// </summary>
    /// <remarks>
    /// <para>Any linear second order differential equation has two independent solutions. For example,
    /// the Bessel differential equation (<see cref="AdvancedMath.Bessel"/>) has solutions J and Y,
    /// the Coulomb wave equation (<see cref="AdvancedMath.Coulomb"/>) has solutions F and G,
    /// and the Airy differential equation (<see cref="AdvancedMath.Airy"/>) has solutions Ai and Bi.</para>
    /// <para>A solution pair structure contains values for both solutions and for their derivatives. It is often useful to
    /// have all this information together when fitting boundary conditions.</para>
    /// <para>Which solution is considered the first and which is considered the second is
    /// a matter of convention. When one solution is regular (finite) at the origin and the other is not, we take the regular solution
    /// to be the first.</para>
    /// </remarks>
    public struct SolutionPair : IEquatable<SolutionPair> {

        private readonly double j, jPrime, y, yPrime;

        /// <summary>
        /// Gets the value of the first solution.
        /// </summary>
        public double FirstSolutionValue {
            get {
                return (j);
            }
        }

        /// <summary>
        /// Gets the derivative of the first solution.
        /// </summary>
        public double FirstSolutionDerivative {
            get {
                return (jPrime);
            }
        }

        /// <summary>
        /// Gets the value of the second solution.
        /// </summary>
        public double SecondSolutionValue {
            get {
                return (y);
            }
        }

        /// <summary>
        /// Gets the derivative of the second solution.
        /// </summary>
        public double SecondSolutionDerivative {
            get {
                return (yPrime);
            }
        }

        // Leaving out the Wronskian for now because it can be subject to extreme cancelation error.
        /*
        /// <summary>
        /// Gets the Wronsikan of the solution pair.
        /// </summary>
        /// <remarks>
        /// <para>The Wronskian of a solution pair is the product of the first solution value and the second solution derivative minus the
        /// product of the second solution value and the first solution derivative.</para>
        /// </remarks>
        public double Wronskian {
            get {
                return (j * yPrime - y * jPrime);
            }
        }
        */

        /// <summary>
        /// Initializes a new solution pair with the given values.
        /// </summary>
        /// <param name="firstSolutionValue">The value of the first solution.</param>
        /// <param name="firstSolutionDerivative">The derivative of the first solution.</param>
        /// <param name="secondSolutionValue">The value of the second solution.</param>
        /// <param name="secondSolutionDerivative">The derivative of the second solution.</param>
        public SolutionPair (double firstSolutionValue, double firstSolutionDerivative, double secondSolutionValue, double secondSolutionDerivative) {
            this.j = firstSolutionValue;
            this.jPrime = firstSolutionDerivative;
            this.y = secondSolutionValue;
            this.yPrime = secondSolutionDerivative;
        }

        // Equality

        public override bool Equals (object obj) {
            if (obj is SolutionPair other) {
                return (Equals(this, other));
            } else {
                return (false);
            }
        }

        public bool Equals (SolutionPair other) {
            return (Equals(this, other));
        }

        public static bool operator == (SolutionPair a, SolutionPair b) {
            return (Equals(a, b));
        }

        public static bool operator != (SolutionPair a, SolutionPair b) {
            return (!Equals(a, b));
        }

        private static bool Equals (SolutionPair a, SolutionPair b) {
            return ((a.j == b.j) && (a.jPrime == b.jPrime) && (a.y == b.y) && (a.yPrime == b.yPrime));
        }

        public override int GetHashCode () {
            unchecked {
                return (j.GetHashCode() + 5 * jPrime.GetHashCode() + 7 * y.GetHashCode() + 13 * yPrime.GetHashCode()); 
            }
        }

        // Fitting
        /*
        public Tuple<double, double> GetCoefficients(double value, double derivative) {
            return ((y * derivative - yPrime * value) / (y * jPrime - j * yPrime)); 
        }
        */
    }

}
