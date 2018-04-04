// Meta.Numerics Library
// Copyright 2014 by David Wright.
// All rights reserved.

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Text;

namespace Meta.Numerics.Functions {

    // A permutation can be represented as a set of cycles, e.g. (1 2)(3), or as a "map" that shows which positions are mapped to which.
    // Some properties are easier to compute in the cycle representation, some easier in the map representation. The public Permutation
    // class can store either representation internally, and generate one from the other if necessary.

    internal class PermutationAsCycles {

        public PermutationAsCycles (int[][] cycles) {
            this.cycles = cycles;
        }

        internal int[][] cycles;

        public int Dimension () {
            // the dimension is the sum of the length of all cycles
            int dimension = 0;
            foreach (int[] cycle in cycles) {
                dimension += cycle.Length;
            }
            return (dimension);
        }

        public bool IsEven () {
            // each even-length cycle contributes a negative sign to the signature
            int count = 0;
            foreach (int[] cycle in cycles) {
                if (cycle.Length % 2 == 0) count++;
            }
            return (count % 2 == 0);
        }

        public bool IsIdentity () {
            // a permutation is the identity if every cycle is of length 1
            foreach (int[] cycle in cycles) {
                if (cycle.Length != 1) return (false);
            }
            return (true);
        }

        public bool IsDerangement () {
            // a permutation is a derangement if there are no length-1 cycles
            foreach (int[] cycle in cycles) {
                if (cycle.Length == 1) return (false);
            }
            return (true);
        }

        public bool IsInvolution () {
            // a permutation is an involution if all cycles are of length 1 or 2
            foreach (int[] cycle in cycles) {
                if (cycle.Length > 2) return (false);
            }
            return (true);
        }

        public long Order () {
            // Since the identity is obtained when every cycle completes at the same time,
            // this is just the LCM of all the cycle lengths.
            // There is no known simple closed expression for this, but it is known to grow
            // like e^{\sqrt{n \ln n}} for large n.
            // See http://mathworld.wolfram.com/LandausFunction.html.
            long order = 1;
            foreach (int[] cycle in cycles) {
                order = (int) AdvancedIntegerMath.LCM(order, cycle.Length);
            }
            // Will overflow for Int32 around n ~ 100, for Int64 around n ~ 300.
            return (order);
        }

        public void Apply<T> (IList<T> x) {
            Debug.Assert(x != null);
            Debug.Assert(x.Count == this.Dimension());
            // knowing cycles allows us to apply in-place
            foreach (int[] cycle in cycles) {
                T t = x[cycle[cycle.Length - 1]];
                for (int i = cycle.Length - 1; i > 0; i--) {
                    x[cycle[i]] = x[cycle[i - 1]];
                }
                x[cycle[0]] = t;
            }
        }

        public override string ToString() {
            StringBuilder text = new StringBuilder();
            foreach (int[] cycle in cycles) {
                text.AppendFormat("({0}", cycle[0]);
                for (int i = 1; i < cycle.Length; i++) {
                    text.AppendFormat(" {0}", cycle[i]);
                }
                text.Append(")");
            }
            return (text.ToString());
        }

        public static bool TryParse (string text, out PermutationAsCycles result) {

            result = null; 

            // Construct cycles list, adding up dimension as we go
            List<int[]> cyclesList = new List<int[]>();
            int start = 0;
            int dimension = 0;
            while (start < text.Length) {
                if (text[start] != '(') return (false);
                int end = text.IndexOf(')', start);
                if (end < 0) return (false);
                string[] cycleText = text.Substring(start + 1, end - start - 1).Split((char[]) null, StringSplitOptions.RemoveEmptyEntries);
                int[] cycle = new int[cycleText.Length];
                dimension += cycle.Length;
                for (int k = 0; k < cycle.Length; k++) {
                    int value;
                    if (!Int32.TryParse(cycleText[k], out value)) return (false);
                    cycle[k] = value;
                }
                cyclesList.Add(cycle);
                start = end + 1;
            }
            int[][] cycles = cyclesList.ToArray();

            // Validate cycles. Each integer should appear once.
            bool[] flags = new bool[dimension];
            for (int i = 0; i < cycles.Length; i++) {
                for (int j = 0; j < cycles[i].Length; j++) {
                    int value = cycles[i][j];
                    if ((value < 0) || (value >= dimension)) return (false);
                    if (flags[value]) return (false);
                    flags[value] = true;
                }
            }

            // At this point the cycles are valid but not necessarily canonical.

            // At this point we want to cannonicalize the cycles. We can do this by
            // shuffling and sorting rows, but a quick-and-dirty solution is to
            // just use the non-canonical cycles to generate a map and let the
            // map to cycles logic generate canonical cycles.

            result = new PermutationAsCycles(cycles);

            return (true);
        }

    }

    internal class PermutationAsMap {

        public PermutationAsMap (int[] map) {
            this.map = map;
        }

        internal int[] map;

        public int Dimension {
            get {
                return (map.Length);
            }
        }

        public bool IsIdentity () {
            // a permutation is the identity if every entry maps to itself
            for (int i = 0; i < map.Length; i++) {
                if (map[i] != i) return (false);
            }
            return (true);
        }

        public bool IsDerangement () {
            // a permutation is a derangement if no entry maps to itself
            for (int i = 0; i < map.Length; i++) {
                if (map[i] == i) return (false);
            }
            return (true);
        }

        public bool IsInvolution () {
            // a permutation is an involution if all entries map back to themselves when applied twice
            for (int i = 0; i < map.Length; i++) {
                if (map[map[i]] != i) return (false);
            }
            return (true);
        }

        public PermutationAsMap Inverse () {
            int[] inverseMap = new int[map.Length];
            for (int i = 0; i < map.Length; i++) {
                inverseMap[map[i]] = i;
            }
            return (new PermutationAsMap(inverseMap));
        }

        public override string ToString () {
            StringBuilder text = new StringBuilder("[");
            text.Append(map[0]);
            for (int i = 1; i < map.Length; i++) {
                text.AppendFormat(" {0}", map[i]);
            }
            text.Append("]");
            return (text.ToString());
        }

        public static bool TryParse (string text, out PermutationAsMap result) {
            result = null;
            if (text.Length < 2) return (false);
            if ((text[0] != '[') || (text[text.Length - 1] != ']')) return (false);
            string[] mapText = text.Substring(1, text.Length - 2).Split((char[]) null, StringSplitOptions.RemoveEmptyEntries);
            int[] map = new int[mapText.Length];
            bool[] flags = new bool[mapText.Length];
            for (int i = 0; i < map.Length; i++) {
                int value;
                if (!Int32.TryParse(mapText[i], out value)) return (false);
                if ((value < 0) || (value > (map.Length - 1))) return (false);
                if (flags[value]) return (false);
                map[i] = value;
                flags[value] = true;
            }
            result = new PermutationAsMap(map);
            return (true);
        }

        public static PermutationAsMap operator * (PermutationAsMap a, PermutationAsMap b) {
            Debug.Assert(a.map.Length == b.map.Length);
            int[] map = new int[a.map.Length];
            for (int i = 0; i < map.Length; i++) {
                map[i] = a.map[b.map[i]];
            }
            return (new PermutationAsMap(map));
        }

        public static PermutationAsMap GetRandomPermutation (int n, Random rng) {

            // Do a Fisher-Yeats shuffle on the integers 0 to n-1

            int[] map = new int[n];
            for (int i = 0; i < map.Length; i++) map[i] = i;

            for (int i = 0; i <map.Length; i++) {
                int j = rng.Next(i, map.Length);
                Global.Swap(ref map[i], ref map[j]);
            }

            return (new PermutationAsMap(map));
        }

        public static PermutationAsMap Identity (int n) {
            int[] map = new int[n];
            for (int i = 0; i < map.Length; i++) map[i] = i;
            return (new PermutationAsMap(map));
        }

        public override int GetHashCode () {
            // Okay to ignore one value, since it is determined by the other values.
            int hash = 19;
            for (int i = 1; i < map.Length; i++) {
                hash = 31 * hash + map[i];
            }
            return (hash);
        }

    }

    /// <summary>
    /// Represents a permutation.
    /// </summary>
    /// <seealso href="http://en.wikipedia.org/wiki/Permutation"/>
    public sealed class Permutation : IEquatable<Permutation>, IFormattable {

        internal Permutation (PermutationAsMap map) {
            this.map = map;
        }

        internal Permutation (PermutationAsCycles cycles) {
            this.cycles = cycles;
        }

        // Some operations work best on the map representation, others best on the cycles representations.
        // We accept either and compute the other as needed.

        private PermutationAsMap map;
        private PermutationAsCycles cycles;

        /// <summary>
        /// Gets the dimension of the permutation.
        /// </summary>
        /// <remarks>
        /// <para>The dimension of a permutation is the number of elements to which the permutation applies.</para>
        /// </remarks>
        public int Dimension {
            get {
                // can get dimension from either map or cycles, but map is faster
                if (map != null) {
                    return (map.Dimension);
                } else {
                    return (cycles.Dimension());
                }
            }
        }

        /// <summary>
        /// Gets a Boolean value that is true if the permutation is even and false if the permutation is odd.
        /// </summary>
        /// <remarks>
        /// <para>An even permutation moves an even number of elements; an odd permutation moves an odd number of elements.</para>
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Even_and_odd_permutations"/>
        public bool IsEven {
            get {
                if (this.cycles == null) ComputeCycles();
                return (cycles.IsEven());
            }
        }

        private string ToStringAsMap () {
            if (map == null) ComputeMap();
            return (map.ToString());
        }

        private string ToStringAsCycles () {
            if (cycles == null) ComputeCycles();
            return (cycles.ToString());
        }

        /// <inheritdoc />
        public override string ToString () {
            return (ToStringAsMap());
        }

        /// <summary>
        /// Converts the permutation to its string representation in the given format.
        /// </summary>
        /// <param name="format">A standard or custom permutation format string.</param>
        /// <returns>The requested string representation of the permutation.</returns>
        /// <remarks>
        /// <para>The standard permutation format strings are "M", which produces a map representation, and "C", which produces a cycle representation.
        /// For explanations of the map and cycle representations of a permutation, see <see cref="Parse"/>.</para>
        /// </remarks>
        public string ToString (string format) {
            return (ToString(format, null));
        }

        /// <summary>
        /// Converts the permutation to its string representation in the given format.
        /// </summary>
        /// <param name="format">A standard or custom permutation format string.</param>
        /// <param name="formatProvider">An object that provides culture-specific formatting information.</param>
        /// <returns>The requested string representation of the permutation.</returns>
        /// <remarks>
        /// <para>The standard permutation format strings are "M", which produces a map representation, and "C",
        /// which produces a cycle representation.
        /// For explanations of the map and cycle representations of a permutation, see <see cref="Parse"/>.</para>
        /// </remarks>
        public string ToString (string format, IFormatProvider formatProvider) {

            if (String.IsNullOrEmpty(format)) format = "M";
            if (formatProvider == null) formatProvider = CultureInfo.CurrentCulture;

            switch (format.ToUpperInvariant()) {
                case "M":
                    return (ToStringAsMap());
                case "C":
                    return (ToStringAsCycles());
                default:
                    throw new FormatException();
            }
        }

        /// <summary>
        /// Converts a text representation into a permutation.
        /// </summary>
        /// <param name="text">A text representation of the permutation.</param>
        /// <returns>The corresponding permutation.</returns>
        /// <remarks>
        /// <para>This method is able to parse both map representations and cycle representations of permutations.</para>
        /// <para>A map representation of an n-dimensional permutation is a space-separated list of all integers between 0 and n-1,
        /// enclosed in square brackets. Each number indicates the index of the location to which the object that appears at
        /// that location is mapped by the permutation. For example, [2 1 0] denotes the permutation that moves the object
        /// at index 0 to index 2, does not move the object at index 1, and moves the object at index 2 to index 0. Note
        /// that the numbers in the map representation are the same as the numbers on the second line of Cauchy's two-line
        /// notation.</para>
        /// <para>A cycle representation of an n-dimensional representation is a space-separated list of all integers between 0 and n-1,
        /// grouped into cycles by parenthesis. Each cycle indicates that the element at the location with the first index in the cycle is moved to
        /// the location with the second index in the cycle, the element at the location with the second index in the cycle is moved
        /// to the location with the third index in the cycle, and so on, until the element at the location with the last index
        /// is moved to the location with the first index. Thus (0 2)(1) indicates that the elements at locations 0 and 2 change
        /// places and the element at location 1 is left there. So (0 2)(1) and [2 1 0] represent the same permutation.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="text"/> is null.</exception>
        /// <exception cref="FormatException"><paramref name="text"/> is not a valid text representation of a permutation.</exception>
        public static Permutation Parse (string text) {
            if (text == null) throw new ArgumentNullException("text");
            Permutation output;
            if (TryParse(text, out output)) {
                return (output);
            } else {
                throw new FormatException();
            }
        }

        /// <summary>
        /// Attempts to convert a text representation into a permutation.
        /// </summary>
        /// <param name="text">A text representation of the permutation.</param>
        /// <param name="output">The corresponding permutation.</param>
        /// <returns>True if the conversion succeeded, otherwise false.</returns>
        /// <remarks>
        /// <para>For information on supported text representations, see <see cref="Parse"/>.</para>
        /// </remarks>
        public static bool TryParse (string text, out Permutation output) {

            output = null;

            if (text == null) return (false);

            bool success = false;
            if ((text.Length > 0) && (text[0] == '[')) {
                // try to parse as ordering
                PermutationAsMap map;
                success = PermutationAsMap.TryParse(text, out map);
                if (success) output = new Permutation(map);
            } else {
                // try to parse as cycles
                PermutationAsCycles cycles;
                success = PermutationAsCycles.TryParse(text, out cycles);
                if (success) output = new Permutation(cycles);
            }

            return (success);

        }

        private static int[][] FromMapToCycles (int[] map) {
            List<int[]> cyclesList = new List<int[]>();
            for (int i = 0; i < map.Length; i++) {
                int j = map[i];
                while (j > i) j = map[j];
                if (j < i) continue;
                List<int> cycleList = new List<int>();
                do {
                    cycleList.Add(j);
                    j = map[j];
                } while (j != i);
                cyclesList.Add(cycleList.ToArray());
            }
            int[][] cycles = cyclesList.ToArray();
            return (cycles);
        }

        private static int[] FromCyclesToMap (int[][] cycles) {
            int dimension = 0;
            foreach (int[] cycle in cycles) dimension += cycle.Length;
            int[] map = new int[dimension];
            // This is not the same as starting with [0 1 2 \cdots n] and applying the permutation.
            // That would yield what Mathematica calls the "list representation". Instead we
            // write into each map[k] the index that value came from.
            foreach (int[] cycle in cycles) {
                for (int i = cycle.Length - 1; i > 0; i--) {
                    map[cycle[i - 1]] = cycle[i];
                }
                map[cycle[cycle.Length - 1]] = cycle[0];
            }
            return (map);
        }

        private void ComputeCycles () {
            Debug.Assert(map != null);
            cycles = new PermutationAsCycles(FromMapToCycles(map.map));
        }

        private void ComputeMap () {
            Debug.Assert(cycles != null);
            map = new PermutationAsMap(FromCyclesToMap(cycles.cycles));
        }


        /// <summary>
        /// Applies the permutation to a list.
        /// </summary>
        /// <typeparam name="T">The type of the list.</typeparam>
        /// <param name="x">The list.</param>
        public void Apply<T> (IList<T> x) {
            if (x == null) throw new ArgumentNullException("x");
            if (x.Count != this.Dimension) throw new DimensionMismatchException();
            if (cycles == null) ComputeCycles();
            cycles.Apply(x);
            // Obviously it's also possible to apply a permutation given its map representation, but doing so requires auxiluary storage.
        }

        /// <summary>
        /// Gets the inverse of the permutation.
        /// </summary>
        /// <returns>The inverse of the permutation.</returns>
        public Permutation Inverse () {
            if (map == null) ComputeMap();
            return (new Permutation(map.Inverse()));
        }

        /// <summary>
        /// Multiplies two permutations.
        /// </summary>
        /// <param name="a">The first permutation.</param>
        /// <param name="b">The second permutation.</param>
        /// <returns>The product permutation ab.</returns>
        /// <remarks>
        /// <para>The product ab means first applying b, then applying a. This right-to-left convention arises from the convention that operators are applied to the right.</para>
        /// </remarks>
        public static Permutation operator * (Permutation a, Permutation b) {
            if (a == null) throw new ArgumentNullException("a");
            if (b == null) throw new ArgumentNullException("b");
            if (a.Dimension != b.Dimension) throw new DimensionMismatchException();

            if (a.map == null) a.ComputeMap();
            if (b.map == null) b.ComputeMap();

            return (new Permutation(a.map * b.map));
        }

        /// <summary>
        /// Gets the order of the permutation.
        /// </summary>
        /// <remarks>
        /// <para>The order of a permutation is the number of times it must be applied in order to return all elements to their original position.
        /// Stated differently, the order of a permutation is the smallest power to which it must be raised to obtain the identity permutation.</para>
        /// <para>Some permutations with dimension greater than about 300 have an order larger than <see cref="Int64.MaxValue"/>; for these permutations
        /// the returned value will overflow.</para>
        /// <para>Note that the word order is also used to refer to the number of distinct permutations of a given dimension. That "order" is a property
        /// of the permutation group. This "order" is a property of each permutation.</para>
        /// </remarks>
        public long Order {
            get {
                if (cycles == null) ComputeCycles();
                return (cycles.Order());
            }
        }

        /// <summary>
        /// Gets a Boolean value indicating whether the permutation is the identity.
        /// </summary>
        /// <value>True if the permutation is the identity, otherwise false.</value>
        /// <remarks>
        /// <para>The identity permutation is the permutation that leaves all elements in their original positions.</para>
        /// </remarks>
        public bool IsIdentity {
            get {
                if (cycles != null) {
                    return(cycles.IsIdentity());
                } else {
                    return (map.IsIdentity());
                }
            }
        }

        /// <summary>
        /// Gets a Boolean flag indicating whether the permutation is a derangement.
        /// </summary>
        /// <value>True if the permutation is a derangement, otherwise false.</value>
        /// <remarks>
        /// <para>A derangement is a permutation that leaves no element in its original position.</para>
        /// </remarks>
        public bool IsDerangement {
            get {
                // can be computed from either representation, but is ever-so-slightly faster (fewer comparisons) from cycles representation
                if (cycles != null) {
                    return (cycles.IsDerangement());
                } else {
                    return (map.IsDerangement());
                }
            }
        }

        /// <summary>
        /// Gets a Boolean flag indicating whether the permutation is an involution.
        /// </summary>
        /// <value>True if the permutation is an involution, otherwise false.</value>
        /// <remarks>
        /// <para>An involution is a permutation that is its own inverse.</para>
        /// </remarks>
        public bool IsInvolution {
            get {
                if (cycles != null) {
                    return (cycles.IsInvolution());
                } else {
                    return (map.IsInvolution());
                }
            }
        }

        /// <summary>
        /// Generates all permutations of the given dimension.
        /// </summary>
        /// <param name="dimension">The number of elements on which the permutations act.</param>
        /// <returns>All permutations of the given dimension.</returns>
        /// <remarks>
        /// <para>The number of permutations of dimension n is n!, which increases very rapidly as n increases. Even in cases
        /// where n! would overflow a <see cref="Int32"/> or <see cref="Int64"/>, this method will successfully produce all permutations.
        /// Of course, in such cases, it will take a long time to enumerate them all.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="dimension"/> is negative.</exception>
        public static IEnumerable<Permutation> GetPermutations (int dimension) {
            if (dimension < 0) throw new ArgumentOutOfRangeException("dimension");

            // The permutations are enumerated using the Steinhaus-Johnson-Trotter algorithm with Even's speedup.
            // See http://en.wikipedia.org/wiki/Steinhaus%E2%80%93Johnson%E2%80%93Trotter_algorithm and Knuth 7.2.1.2.
            // The Wikipedia article is less than completely illuminating.
            // See http://tropenhitze.wordpress.com/2010/01/25/steinhaus-johnson-trotter-permutation-algorithm-explained-and-implemented-in-java/
            // and http://www.coderslexicon.com/johnson-trotter-algorithm-in-vc-net/ for more detailed descriptions.

            // Note that permutations are not produced in lexographic order. They do alternate between even and odd because there is a single transposition between each.

            // Initialize values in ascending order, all moving left.
            int[] values = new int[dimension];
            bool[] directions = new bool[dimension];
            for (int i = 0; i < dimension; i++) {
                values[i] = i;
                directions[i] = false;
            }

            while (true) {

                // Return the current permutations,
                // copying first so subsequent changes to values don't alter the returned Permutation object.
                int[] map = new int[dimension];
                Array.Copy(values, map, dimension);
                yield return (new Permutation(new PermutationAsMap(map)));

                // Find the largest mobile value, stopping if none is found.
                int index = FindLargestMobile(values, directions);
                if (index < 0) break;
                int value = values[index];

                // Switch it with the adjacent value in its direction.
                if (directions[value]) {
                    Global.Swap(ref values[index], ref values[index + 1]);
                } else {
                    Global.Swap(ref values[index], ref values[index - 1]);
                }

                // Switch the direction of all higher values.
                for (int j = value + 1; j < directions.Length; j++) {
                    directions[j] = !directions[j];
                }

            }

        }

        // This subroutine finds the largest mobile integer, which is the one that will be swapped in the SJT algorithm.
        // An integer is mobile if the adjacent value in its direction is less than itself.

        // Note directions[i] is the direction of the value i, not the direction of the value at values[i].
        // This convention allows us to reduce the number of required operations.

        private static int FindLargestMobile (int[] values, bool[] directions) {

            int largestIndex = -1;
            int largestValue = Int32.MinValue;

            for (int index = 0; index < values.Length; index++) {
                int value = values[index];
                if (value > largestValue) {
                    if (directions[value]) {
                        // right-moving value, compare number to right
                        if (((index + 1) < values.Length) && (values[index + 1] < value)) {
                            largestIndex = index;
                            largestValue = value;
                        }
                    } else {
                        // left-moving value, compare number to left
                        if ((index > 0) && (values[index - 1] < value)) {
                            largestIndex = index;
                            largestValue = value;
                        }
                    }
                }

            }

            return (largestIndex);

        }

        /// <summary>
        /// Returns the identity permutation of the given dimension.
        /// </summary>
        /// <param name="dimension">The number of elements on which the permutation acts.</param>
        /// <returns>The identity permutation of the requested dimension.</returns>
        public static Permutation Identity (int dimension) {
            if (dimension < 0) throw new ArgumentOutOfRangeException("dimension");
            return (new Permutation(PermutationAsMap.Identity(dimension)));
        }

        // This internal equals test should only be called with non-null arguments.

        private static bool InstanceEquals (Permutation a, Permutation b) {
            Debug.Assert(!Object.ReferenceEquals(a, null));
            Debug.Assert(!Object.ReferenceEquals(b, null));
            if (a.Dimension != b.Dimension) return (false);
            if (a.map == null) a.ComputeMap();
            if (b.map == null) b.ComputeMap();
            for (int i = 0; i < a.map.map.Length; i++) {
                if (a.map.map[i] != b.map.map[i]) return (false);
            }
            return (true);
        }

        /// <summary>
        /// Determines whether two permutations are equal.
        /// </summary>
        /// <param name="a">The first permutation.</param>
        /// <param name="b">The second permutation.</param>
        /// <returns>True if the two permutations are equal, otherwise false.</returns>
        public static bool operator == (Permutation a, Permutation b) {
            if (Object.ReferenceEquals(a, b)) {
                return (true);
            } else if (Object.ReferenceEquals(a, null) || Object.ReferenceEquals(b, null)) {
                return (false);
            } else {
                return (InstanceEquals(a, b));
            }
        }

        /// <summary>
        /// Determines whether two permutations are not equal.
        /// </summary>
        /// <param name="a">The first permutation.</param>
        /// <param name="b">The second permutation.</param>
        /// <returns>True if the two permutations are not equal, otherwise true.</returns>
        public static bool operator != (Permutation a, Permutation b) {
            return (!(a == b));
        }

        /// <summary>
        /// Determines whether the given permutation is equal to the permutation instance.
        /// </summary>
        /// <param name="other">The permutation to compare.</param>
        /// <returns>True if <paramref name="other"/> equals the permutation instance, otherwise false.</returns>
        public bool Equals (Permutation other) {
            if (Object.ReferenceEquals(other, null)) {
                return (false);
            } else {
                return (InstanceEquals(this, other));
            }
        }

        /// <summary>
        /// Determines whether the given object is equal to the permutation.
        /// </summary>
        /// <param name="obj">The object to compare.</param>
        /// <returns>True if <paramref name="obj"/> is a permutation equal to the permutation instance, otherwise false.</returns>
        public override bool Equals (object obj) {
            Permutation other = obj as Permutation;
            return (this.Equals(other));
        }

        /// <inheritdoc />
        public override int GetHashCode () {
            if (map == null) ComputeMap();
            return (map.GetHashCode());
        }

        /// <summary>
        /// Get a random permutation.
        /// </summary>
        /// <param name="dimension">The number of elements on which the permutation acts.</param>
        /// <param name="rng">A random number generator.</param>
        /// <returns>A random permutation of the specified dimension. All permutations of the specified dimension are equally likely.</returns>
        public static Permutation GetRandomPermutation (int dimension, Random rng) {
            if (dimension < 1) throw new ArgumentOutOfRangeException("dimension");
            if (rng == null) throw new ArgumentNullException("rng");
            return (new Permutation(PermutationAsMap.GetRandomPermutation(dimension, rng)));
        }

    }

}
