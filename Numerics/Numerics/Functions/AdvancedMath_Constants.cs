using System;

namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {

        // Mathematical constants

        /// <summary>
        /// The golden ratio.
        /// </summary>
        /// <remarks><para>The golden ratio &#x3C6; = 1.1618...</para></remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Golden_ratio"/>
        /// <seealso href="http://mathworld.wolfram.com/GoldenRatio.html" />
        public static readonly double GoldenRatio = (1.0 + Math.Sqrt(5.0)) / 2.0;

        /// <summary>
        /// The Euler constant.
        /// </summary>
        /// <remarks><para>The Euler constant &#x3B3; = 0.5772...</para></remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Euler_gamma"/>
        /// <seealso href="http://mathworld.wolfram.com/Euler-MascheroniConstant.html" />
        public const double EulerGamma = 0.577215664901532860606512;

        /// <summary>
        /// Catalan's constant.
        /// </summary>
        /// <remarks><para>Catalan's constant 0.9159...</para></remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Catalan_constant"/>
        public const double Catalan = 0.915965594177219015054604;


        /// <summary>
        /// Apéry's constant.
        /// </summary>
        /// <remarks><para>Apéry's constant is the value of the Riemann Zeta function for the argument 3.</para></remarks>
        /// <seealso cref="AdvancedMath.RiemannZeta(double)"/>
        /// <seelaso href="https://en.wikipedia.org/wiki/Ap%C3%A9ry%27s_constant"/>
        public const double Apery = 1.20205690315959428540;
    }
}

