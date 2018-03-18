using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {
    /// <summary>
    /// Represents a collection of fit parameters.
    /// </summary>
    public sealed class ParameterCollection : IReadOnlyCollection<Parameter> {

        // One-parameter constructor
        internal ParameterCollection (string name, double p, double vp) {
            Debug.Assert(name != null);
            this.best = new ColumnVector(p);
            this.covariance = new SymmetricMatrix(1);
            this.covariance[0, 0] = vp;
            this.best.IsReadOnly = true;
            this.covariance.IsReadOnly = true;
            this.names = new string[] { name };
            this.map = new Dictionary<string, int>() { { name, 0 } };
        }

        // Two-parameter constructor
        internal ParameterCollection (string name1, double p1, double v1, string name2, double p2, double v2, double c12) {
            Debug.Assert(name1 != null);
            Debug.Assert(name2 != null);
            this.best = new ColumnVector(p1, p2);
            this.covariance = new SymmetricMatrix(2);
            this.covariance[0, 0] = v1;
            this.covariance[1, 1] = v2;
            this.covariance[0, 1] = c12;
            this.best.IsReadOnly = true;
            this.covariance.IsReadOnly = true;
            this.names = new string[] { name1, name2 };
            this.map = new Dictionary<string, int>() {
                {name1, 0 }, {name2, 1 }
            };
        }

        internal ParameterCollection (IReadOnlyList<string> names, ColumnVector best, SymmetricMatrix covariance) {
            Debug.Assert(names != null);
            Debug.Assert(best != null);
            Debug.Assert(covariance != null);
            Debug.Assert(best.Dimension == names.Count);
            Debug.Assert(covariance.Dimension == names.Count);
            best.IsReadOnly = true;
            covariance.IsReadOnly = true;
            this.best = best;
            this.covariance = covariance;
            this.names = names;
            this.map = new Dictionary<string, int>(names.Count);
            for(int i = 0; i < names.Count; i++) {
                this.map.Add(names[i], i);
            }
        }

        internal readonly IReadOnlyList<string> names;
        internal readonly Dictionary<string, int> map;
        internal readonly ColumnVector best;
        internal readonly SymmetricMatrix covariance;

        /// <summary>
        /// Gets the number of parameters.
        /// </summary>
        public int Count {
            get {
                return (map.Count);
            }
        }

        /// <summary>
        /// Gets the set of best-fit parameters, as a vector.
        /// </summary>
        public ColumnVector Best {
            get {
                return (best);
            }
        }

        /// <summary>
        /// Gets the covariance matrix of the fit parameters.
        /// </summary>
        public SymmetricMatrix Covariance {
            get {
                return (covariance);
            }
        }

        /// <summary>
        /// Gets the variance of the named parameter.
        /// </summary>
        /// <param name="name">The parameter name.</param>
        /// <returns>The variance of the named parameter.</returns>
        public double VarianceOf (string name) {
            if (name == null) throw new ArgumentNullException(nameof(name));
            int index = IndexOf(name);
            return (covariance[index, index]);
        }

        /// <summary>
        /// Gets the covariance of the two named parameters.
        /// </summary>
        /// <param name="name1">The name of the first parameter.</param>
        /// <param name="name2">The name of the second parameter.</param>
        /// <returns>The covariance of the two named parameters.</returns>
        public double CovarianceOf (string name1, string name2) {
            if (name1 == null) throw new ArgumentNullException(nameof(name1));
            if (name2 == null) throw new ArgumentNullException(nameof(name2));
            int index1 = IndexOf(name1);
            int index2 = IndexOf(name2);
            return (covariance[index1, index2]);
        }

        /// <summary>
        /// Gets the parameter with the given index.
        /// </summary>
        /// <param name="index">The index of the paramter.</param>
        /// <returns>The requested parameter.</returns>
        public Parameter this[int index] {
            get {
                if ((index < 0) || (index >= Count)) throw new ArgumentOutOfRangeException(nameof(index));
                return (new Parameter(this, index));
            }
        }

        /// <summary>
        /// Gets the parameter with the given name.
        /// </summary>
        /// <param name="name">The name of the parameter.</param>
        /// <returns>The requested parameter.</returns>
        public Parameter this[string name] {
            get {
                if (name == null) throw new ArgumentNullException(nameof(name));
                int index = IndexOf(name);
                if (index < 0) throw new ArgumentOutOfRangeException(nameof(name));
                return (new Parameter(this, index));
            }
        }

        /// <summary>
        /// Gets the index of the parameter with the given name.
        /// </summary>
        /// <param name="name">The name of the parameter.</param>
        /// <returns>The index of the parameter, or -1 if no such parameter exists.</returns>
        public int IndexOf (string name) {
            if (name == null) throw new ArgumentNullException(nameof(name));
            int index;
            if (map.TryGetValue(name, out index)) {
                Debug.Assert((0 <= index) && (index < map.Count));
                return (index);
            } else {
                return (-1);
            }
        }

        IEnumerator<Parameter> IEnumerable<Parameter>.GetEnumerator () {
            for (int i = 0; i < Count; i++) {
                yield return new Parameter(this, i);
            }
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (((IEnumerable<Parameter>) this).GetEnumerator());
        }
    }

    /// <summary>
    /// Represents a parameter from a fit.
    /// </summary>
    public sealed class Parameter {

        internal Parameter (ParameterCollection parameters, int index) {
            Debug.Assert(parameters != null);
            Debug.Assert(0 <= index);
            Debug.Assert(index <= parameters.map.Count);
            this.parameters = parameters;
            this.index = index;
        }

        private readonly ParameterCollection parameters;
        private readonly int index;

        /// <summary>
        /// Gets the index of the parameter.
        /// </summary>
        public int Index { get { return (index); } }

        /// <summary>
        /// Gets the name of the parameter.
        /// </summary>
        public string Name { get { return (parameters.names[index]); } }

        /// <summary>
        /// Gets the estimated value and uncertainty of the parameter.
        /// </summary>
        public UncertainValue Estimate {
            get {
                return (new UncertainValue(parameters.best[index], Math.Sqrt(parameters.covariance[index, index])));
            }
        }

        /// <summary>
        /// Returns a string representation of the parameter.
        /// </summary>
        /// <returns>A string representation of the parameter.</returns>
        public override string ToString () {
            return ($"{Name} = {Estimate}");
        }

    }

}
