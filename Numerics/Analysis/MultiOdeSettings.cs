using System;
using System.Collections.Generic;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Contains settings used to solve a set of coupled ordinary differential equations.
    /// </summary>
    public class MultiOdeSettings : EvaluationSettings {

        /// <summary>
        /// Initializes a new instance of evaluation settings for coupled ODEs.
        /// </summary>
        public MultiOdeSettings () : base() { }


        /// /// <summary>
        /// Gets or sets the handler that is called to report on progress toward the solution.
        /// </summary>
        /// <remarks>
        /// <para>As the ODE is integrated, the specified handler is called for each integration
        /// step. This allows the handler to monitor progress toward the solution and to
        /// obtain intermediate values via interpolation.</para>
        /// </remarks>
        public Action<MultiOdeResult> Listener { get; set; }

    }

}
