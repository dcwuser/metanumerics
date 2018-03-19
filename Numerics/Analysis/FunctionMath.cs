using System;
using System.Collections.Generic;
using System.Text;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Contains methods for the analysis of functions.
    /// </summary>
    /// <remarks>
    /// <para>This is the primary class for the numerical analysis of functions of a single variable.
    /// Function analysis includes integration, optimization (finding maxima and minima), locating roots,
    /// and solving differential equations.</para>
    /// <para>The termination criteria and evaluation budget for numerical analysis can be controlled
    /// using the <see cref="EvaluationSettings"/> class and its child classes.</para>
    /// <para>For the analysis of multi-dimensional functions, see the <see cref="MultiFunctionMath"/> class.</para>
    /// </remarks>
    public static partial class FunctionMath {

        // logic for finding function points is in other classes

        private static readonly double RelativePrecision = Math.Pow(2.0, -48);

        private static readonly double AbsolutePrecision = Math.Pow(2.0, -192);

    }

}
