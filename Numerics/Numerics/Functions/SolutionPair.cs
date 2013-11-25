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
    /// the Coulomb wave equation has solutions F and G,
    /// and the Airy differential equation has solutions Ai and Bi.</para>
    /// <para>A solution pair structure contains values for both solutions and for their derivatives. It is often useful to
    /// have all this information together when fitting boundary conditions.</para>
    /// <para>Which solution is considered the first and which is considered the second is
    /// a matter of convention. When one solution is regular (finite) at the origin and the other is not, we take the regular solution
    /// to be the first.</para>
    /// </remarks>
    public struct SolutionPair {

        private double j, jPrime, y, yPrime;

        /// <summary>
        /// Gets the value of the first solution.
        /// </summary>
        public double FirstSolutionValue {
            get {
                return (j);
            }
            internal set {
                j = value;
            }
        }

        /// <summary>
        /// Gets the derivative of the first solution.
        /// </summary>
        public double FirstSolutionDerivative {
            get {
                return (jPrime);
            }
            internal set {
                jPrime = value;
            }
        }

        /// <summary>
        /// Gets the value of the second solution.
        /// </summary>
        public double SecondSolutionValue {
            get {
                return (y);
            }
            internal set {
                y = value;
            }
        }

        /// <summary>
        /// Gets the derivative of the second solution.
        /// </summary>
        public double SecondSolutionDerivative {
            get {
                return (yPrime);
            }
            internal set {
                yPrime = value;
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

        internal SolutionPair (double j, double jPrime, double y, double yPrime) {
            this.j = j;
            this.jPrime = jPrime;
            this.y = y;
            this.yPrime = yPrime;
        }

    }

}
