using System;
using System.Diagnostics;

namespace Meta.Numerics.Matrices
{
    internal class AnyRectangularMatrixDebuggerTypeProxy
    {

        public AnyRectangularMatrixDebuggerTypeProxy(AnyRectangularMatrix A) {
            this.A = A;
        }
    
        private AnyRectangularMatrix A;

        [DebuggerBrowsable(DebuggerBrowsableState.RootHidden)]
        public double[,] Entries {
            get {
                return (A.ToArray());
            }
        }

    }
}
