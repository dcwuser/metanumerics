using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

using Meta.Numerics.Functions;

namespace Examples
{
    public class Permutations
    {

        public static void EnumeratePermutations() {

            foreach (Permutation permutation in Permutation.GetPermutations(4)) {
                Console.WriteLine(permutation);
            }

        }

    }
}
