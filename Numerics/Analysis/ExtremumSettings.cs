using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Contains settings used to control optimization.
    /// </summary>
    public class ExtremumSettings : EvaluationSettings {

        /// <summary>
        /// Gets or sets the handler that is called to report on progress toward the optimum.
        /// </summary>
        public Action<Extremum> Listener { get; set; }

    }

    /// <summary>
    /// Contains settings used to control multi-dimensional optimization.
    /// </summary>
    public class MultiExtremumSettings : EvaluationSettings {

        /// <summary>
        /// Gets or sets the handler than is called to report on progress toward the optimum.
        /// </summary>
        public Action<MultiExtremum> Listener { get; set; }

    }
}
