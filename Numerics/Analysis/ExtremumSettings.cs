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
}
