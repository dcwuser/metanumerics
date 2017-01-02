using System;
using System.Collections.Generic;
using System.Text;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Contains methods for the analysis of functions.
    /// </summary>
    /// <remarks>
    /// <para>Function analysis includes integration, finding maxima and minima, and finding roots.</para>
    /// <para>This class contains methods for the analysis of functions that both accept a single real argument and return a single real value.
    /// For the analysis of multi-dimensional functions, see the <see cref="MultiFunctionMath"/> class.</para>
    /// </remarks>
    public static partial class FunctionMath {

        // logic for finding function points is in other classes

        private static readonly double RelativePrecision = Math.Pow(2.0, -48);

        private static readonly double AbsolutePrecision = Math.Pow(2.0, -192);

    }

}
