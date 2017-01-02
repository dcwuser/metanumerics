using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Threading;

namespace Meta.Numerics.Analysis {

    internal class MultiFunctor {

        private readonly Func<IList<double>, double> function;

        private readonly CoordinateTransform[] map;

        private readonly bool negate = false;

        private int count = 0;

        public MultiFunctor (Func<IList<double>, double> function, bool negate) : this(function) {
            this.negate = negate;
        }

        public MultiFunctor (Func<IList<double>, double> function) : this(function, null) { }

        public MultiFunctor (Func<IList<double>, double> function, CoordinateTransform[] map) {
            this.function = function;
            this.map = map;
        }

        public bool IgnoreInfinity { get; set; }

        public bool IgnoreNaN { get; set; }

        public double Evaluate (double[] x) {
            Interlocked.Increment(ref count);
            if (map == null) {
                double z = function(new ReadOnlyCollection<double>(x));
                if (IgnoreInfinity && Double.IsInfinity(z)) z = 0.0;
                if (IgnoreNaN && Double.IsNaN(z)) z = 0.0;
                if (negate) z = -z;
                return (z);
            } else {
                double j = 1.0;
                for (int i = 0; i < x.Length; i++) {
                    CoordinateTransform transform = map[i];
                    if (transform != null) transform.TransformInPlace(ref x[i], ref j);
                }
                double z = function(new ReadOnlyCollection<double>(x));
                if (IgnoreInfinity && Double.IsInfinity(z)) z = 0.0;
                if (IgnoreNaN && Double.IsNaN(z)) z = 0.0;
                return (j * z);
            }

        }

        public int EvaluationCount {
            get {
                return (count);
            }
        }

        public bool IsNegated {
            get {
                return (negate);
            }
        }

    }


}
