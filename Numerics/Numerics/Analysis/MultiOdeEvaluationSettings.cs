using System;
using System.Collections.Generic;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Contains settings used to solve a multi-dimensional ordinary differential equation.
    /// </summary>
    public class MultiOdeEvaluationSettings : EvaluationSettings {

        /// <summary>
        /// Initializes a new instance of evaluation settings for multiple ODEs.
        /// </summary>
        public MultiOdeEvaluationSettings () : base(null) { }


        /// /// <summary>
        /// Gets or sets the handler that is called to report on progress toward the solution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// For each successful step, this delegate is called to report the solution at the
        /// new abcissa. If you have the listener record these results, we can interpolate
        /// in order to approximate the solution at intermediate points.
        /// </para>
        /// </remarks>
        public Action<MultiOdeResult> Listener { get; set; }

    }

}
