using System;
using System.Collections.Generic;
using System.Diagnostics;


namespace Meta.Numerics.Functions {
    public static partial class AdvancedIntegerMath {

        /// <summary>
        /// Enumerates all partitions of the given integer
        /// </summary>
        /// <param name="n">The integer to partition, which must be positive.</param>
        /// <returns>An enumeration of all partitions of the given integer.</returns>
        /// <remarks>
        /// <para>Integer partitions are ways to write an integer as a sum of smaller integers. For example, the integer 4 has 5 partitions: 4,
        /// 3 + 1, 2 + 2, 2 + 1 + 1, and 1 + 1 + 1 + 1.</para>
        /// <para>Integer partitions appear in combinatoric problems and solutions to problems that may be mapped into combinatoric problems.
        /// For example, the terms which appear in <a href="http://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno%27s_formula">Faà di Bruno's formula</a>
        /// correspond to integer partitions.</para>
        /// <para>The number of partitions grows very rapidly with n. Since enumerating through partitions does not require us to count them,
        /// no overflows will occur even for large values of <paramref name="n"/>. However, completing the enumeration of
        /// such a large number of partitions will take a long time, even though our algorithm produces each partition very quickly. For
        /// example, there are about two hundred million partitions of the integer 100.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is not positive.</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Integer_partition"/>
        internal static IEnumerable<IntegerPartition> Partitions (int n) {

            if (n < 1) throw new ArgumentOutOfRangeException(nameof(n));

            foreach (int[] p in InternalPartitions(n)) {
                yield return (new IntegerPartition(p));
            }

        }

        internal static IEnumerable<int[]> InternalPartitions (int n) {

            Debug.Assert(n > 0);

            // This is algorithm 5.3, called RuleAsc, from Jerome Kelleher's thesis
            // "Encoding Paritions As Ascending Compositions", University College Cork, 2005

            // He analyses and compares a large number of algorithms for the generation
            // of restricted and unrestricted partitions. This one is good.

            // He also describes another one that is considerably more complex and just
            // slightly better. We don't implement it.

            // Initialize the state
            int[] a = new int[n + 1];
            int k = 1;
            a[0] = 0;
            a[1] = n;

            while (k != 0) {

                // Advance to the next partition
                int y = a[k] - 1;
                k--;
                int x = a[k] + 1;
                while (x <= y) {
                    a[k] = x;
                    y -= x;
                    k++;
                }
                a[k] = x + y;

                // Return a copy so our state array can't be disturbed
                int[] p = new int[k + 1];
                Array.Copy(a, p, k + 1);
                yield return (p);
            }

        }

    }
}
