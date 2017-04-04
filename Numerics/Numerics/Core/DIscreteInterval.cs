using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics {
    public struct DiscreteInterval {

        public DiscreteInterval (int min, int max) {
            this.min = min;
            this.max = max;
        }

        private readonly int min;

        private readonly int max;

        public int LeftEndpoint {
            get {
                return (min);
            }
        }

        public int RightEndpoint {
            get {
                return (max);
            }
        }

    }
}
