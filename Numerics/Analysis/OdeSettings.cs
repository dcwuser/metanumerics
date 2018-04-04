using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Contains settings controlling the integration of an ordinary differential equation.
    /// </summary>
    public class OdeSettings : EvaluationSettings {

        /// <summary>
        /// Initializes a new instance of ODE evaluation settings.
        /// </summary>
        public OdeSettings () : base() {

        }

        /// <summary>
        /// Gets or sets the handler that is called to report on progress toward the solution.
        /// </summary>
        /// <remarks>
        /// <para>As the ODE is integrated, the specified handler is called for each integration
        /// step. This allows the handler to monitor progress toward the solution and to
        /// obtain intermediate values via interpolation.</para>
        /// </remarks>
        public Action<OdeResult> Listener { get; set; }

    }
}
