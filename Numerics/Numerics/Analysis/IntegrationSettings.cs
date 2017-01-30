using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Contains settings used to control the evaluation of integrals.
    /// </summary>
    public class IntegrationSettings : EvaluationSettings {

        /// <summary>
        /// A handler that will be called with updates as an integral is evaluated.
        /// </summary>
        public Action<IntegrationResult> Listener { get; set; }

    }
}
