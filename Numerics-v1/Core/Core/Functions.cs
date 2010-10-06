using System;
using System.Collections.Generic;
using System.Text;

namespace Meta.Numerics {

    /// <summary>
    /// Represents a function of one input parameter.
    /// </summary>
    /// <typeparam name="TIn">The type of the input parameter.</typeparam>
    /// <typeparam name="TOut">The type of the output value.</typeparam>
    /// <param name="x">The input parameter.</param>
    /// <returns>The output value.</returns>
    public delegate TOut Function<TIn, TOut> (TIn x);

    /// <summary>
    /// Represents a function of two input parameters.
    /// </summary>
    /// <typeparam name="TIn1">The type of the first input aparameter.</typeparam>
    /// <typeparam name="TIn2">The type of the second input parameter.</typeparam>
    /// <typeparam name="TOut">The type of the output value.</typeparam>
    /// <param name="x1">The first input parameter.</param>
    /// <param name="x2">The second input parameter.</param>
    /// <returns>The output value.</returns>
    public delegate TOut Function<TIn1, TIn2, TOut> (TIn1 x1, TIn2 x2);

}
