using System;
using System.Collections.Generic;
using System.Text;

namespace Meta.Numerics.Functions
{
    public static partial class AdvancedMath
    {

        /// <summary>
        /// Computes the requested zero of the Airy Ai function.
        /// </summary>
        /// <param name="k">The index of the zero.</param>
        /// <returns>The <paramref name="k"/>th value of x for which Ai(x) = 0.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> is less than 1.</exception>
        public static double AiryAiZero (int k) {
            return (BesselMath.AiryAiZero(k));
        }

        /// <summary>
        /// Computes the requested zero of the Airy Bi (Bairy) function.
        /// </summary>
        /// <param name="k">The index of the zero.</param>
        /// <returns>The <paramref name="k"/>th value of x for which Bi(x) = 0.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> is less than 1.</exception>
        public static double AiryBiZero (int k) {
            return (BesselMath.AiryBiZero(k));
        }

    }
}
