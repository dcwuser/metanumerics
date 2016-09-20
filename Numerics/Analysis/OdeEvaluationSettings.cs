using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Contains settings controlling the solution of an ordinary differential equation.
    /// </summary>
    public class OdeEvaluationSettings : EvaluationSettings {

        /// <summary>
        /// Gets or sets the handler that is called to report on the progress toward the solution.
        /// </summary>
        /// <remarks>
        /// <para>As the ODE is integrated, the specified handler is called for each integration
        /// step. This allows the handler to monitor progress toward the solution and to
        /// obtain intermediate solution values.</para>
        /// </remarks>
        public Action<OdeResult> UpdateHandler { get; set; }

    }
}
