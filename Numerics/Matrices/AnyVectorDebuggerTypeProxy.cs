using System;
using System.Diagnostics;

namespace Meta.Numerics.Matrices
{
    internal class AnyVectorDebuggerTypeProxy
    {
        public AnyVectorDebuggerTypeProxy(AnyVector v) {
            this.v = v;
        }

        private AnyVector v;

        [DebuggerBrowsable(DebuggerBrowsableState.RootHidden)]
        public double[] Entries {
            get {
                return (v.ToArray());
            }
        }
    }
}
