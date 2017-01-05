﻿using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Functions {

    /// <summary>
    /// Represents an integer partition.
    /// </summary>
    public class IntegerPartition : IEquatable<IntegerPartition> {

        internal IntegerPartition (int[] values) {
            Debug.Assert(values != null);
            this.values = values;
        }

        private int[] values;
        private List<Element> elements;

        // Values is assumed to be in standard, ascending form, e.g.
        // {1, 1, 1, 2, 4, 4} represents 1 + 1 + 1 + 2 + 4 + 4 = 13

        private void ComputeElements () {
            Debug.Assert(values != null);
            List<Element> elements = new List<Element>();
            if (values.Length > 0) {
                int u = values[0];
                int m = 1;
                for (int i = 1; i < values.Length; i++) {
                    int v = values[i];
                    if (v == u) {
                        m++;
                    } else {
                        elements.Add(new Element(u, m));
                        u = v;
                        m = 1;
                    }
                }
                elements.Add(new Element(u, m));
            }
            this.elements = elements;
        }

        /// <summary>
        /// Gets the values in the partition.
        /// </summary>
        /// <value>A read-only list of values that add up to the partitioned integer.</value>
        public IList<int> Values {
            get {
                return (new ReadOnlyCollection<int>(values));
            }
        }

        /// <summary>
        /// Gets the elements of the partition.
        /// </summary>
        /// <value>A read-only list of elements that make up the partition.</value>
        /// <remarks>
        /// <para>This is the multiplicity representation of the partition.</para>
        /// </remarks>
        public IList<Element> Elements {
            get {
                if (elements == null) ComputeElements();
                return (new ReadOnlyCollection<Element>(elements));
            }
        }

        /// <summary>
        /// Gets the rank of the partition.
        /// </summary>
        public int Rank {
            get {
                int l = values.Length;
                return (values[l - 1] - l);
            }
        }

        /// <summary>
        /// Computes the conjugate partition.
        /// </summary>
        /// <returns>The conjugate of the partition.</returns>
        public IntegerPartition Conjugate () {

            // The conjugate of a partition is produced by reading off the Ferrers diagrams along the opposite axis.
            // For example, the conjugate of 1+1+1+2+4+4 is 2+2+3+6

            //   442111
            // 6 XXXXXX
            // 3 XXX
            // 2 XX
            // 2 XX

            // The largest value of the conjugate is thus always the number of values of the original. We then
            // proceed to fill in the remainder by reducing the value so that we can still reach the remaining
            // blocks.

            Debug.Assert(values != null);
            Debug.Assert(values.Length > 0);

            int[] cValues = new int[values[values.Length - 1]];

            int cValue = values.Length;
            for (int i = 1; i <= cValues.Length; i++) {
                cValues[cValues.Length - i] = cValue;
                while (values[values.Length - cValue] < (i + 1)) {
                    cValue--;
                    if (cValue == 0) break;
                }
            }

            return (new IntegerPartition(cValues));
        }

        // Equality

        private static bool InternalEquals (IntegerPartition a, IntegerPartition b) {
            Debug.Assert(!Object.ReferenceEquals(a, null));
            Debug.Assert(!Object.ReferenceEquals(b, null));
            if (a.values.Length != b.values.Length) return (false);
            for (int i = 0; i < a.values.Length; i++) {
                if (a.values[i] != b.values[i]) return (false);
            }
            return (true);
        }

        /// <summary>
        /// Determines whether two partitions are equal.
        /// </summary>
        /// <param name="a">The first partition.</param>
        /// <param name="b">The second partition.</param>
        /// <returns><see langword="true"/> if <paramref name="a"/> and <paramref name="b"/> represent the same
        /// partition, otherwise <see langword="false"/>.</returns>
        public static bool operator == (IntegerPartition a, IntegerPartition b) {
            if (Object.ReferenceEquals(a, b)) {
                return (true);
            } else if (Object.ReferenceEquals(a, null) || Object.ReferenceEquals(b, null)) {
                return (false);
            } else {
                return (InternalEquals(a, b));
            }
        }

        /// <summary>
        /// Determines whether two partitions are not equal.
        /// </summary>
        /// <param name="a">The first partition.</param>
        /// <param name="b">The second partition.</param>
        /// <returns><see langword="true"/> if <paramref name="a"/> and <paramref name="b"/> represent different
        /// partitions, otherwise <see langword="false"/>.</returns>
        public static bool operator != (IntegerPartition a, IntegerPartition b) {
            return (!(a == b));
        }

        /// <summary>
        /// Determines whether another partition is equal to this one.
        /// </summary>
        /// <param name="other">The other partition.</param>
        /// <returns><see langword="true"/> is <paramref name="other"/> represents the same partition, otherwise
        /// <see langword="false"/>.</returns>
        public bool Equals (IntegerPartition other) {
            if (Object.ReferenceEquals(other, null)) {
                return (false);
            } else {
                return (InternalEquals(this, other));
            }
        }

        /// <inheritdoc />
        public override bool Equals (object obj) {
            return (Equals(obj as IntegerPartition));
        }

        /// <inheritdoc />
        public override int GetHashCode () {
            int hash = values.Length;
            int max = Math.Min(values.Length, 16);
            for(int i = 0; i < max; i++) {
                hash = values[i] + 17 * hash;
            }
            return (hash);
        }

        /// <inheritdoc />
        public override string ToString () {
            StringBuilder text = new StringBuilder();
            text.Append(values[0]);
            for (int i = 1; i < values.Length; i++) {
                text.AppendFormat("+{0}", values[i]);
            }
            return (text.ToString());
        }
    }

    /// <summary>
    /// Describes the multiplicity of an integer in a set.
    /// </summary>
    public struct Element {

        internal Element (int value, int multiplicity) {
            Debug.Assert(multiplicity > 0);
            this.value = value;
            this.multiplicity = multiplicity;
        }

        private readonly int value;

        private readonly int multiplicity;

        /// <summary>
        /// Gets the value of the integer.
        /// </summary>
        public int Value {
            get {
                return (value);
            }
        }

        /// <summary>
        /// Gets the multiplicity of the integer.
        /// </summary>
        public int Multiplicity {
            get {
                return (multiplicity);
            }
        }

    }
}
