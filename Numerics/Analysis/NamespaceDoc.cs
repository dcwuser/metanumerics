using System.Runtime.CompilerServices;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Contains types used to solve, maximize, integrate, and otherwise perform analysis on user-supplied functions.
    /// </summary>
    /// <remarks>
    /// <para>This namespace contains types for the numerical analysis of user-supplied functions. Examples
    /// of such analysis include:</para>
    /// <list type="bullet">
    /// <item>evaluating a definite integral</item>
    /// <item>solving transcendental equations (i.e. finding roots)</item>
    /// <item>finding an optimum (minimum or maximum)</item>
    /// <item>solving a differential equation</item>
    /// </list>
    /// <para>The central types in the namespace are the static classes <see cref="FunctionMath"/>,
    /// which contains methods for the analysis of functions of a single variable, and
    /// <see cref="MultiFunctionMath"/>, which contains methods for the analysis of functions with
    /// multidimensional arguments or outputs.</para>
    /// <para>The termination criteria for most analysis methods can be controlled by providing
    /// an appropriate <see cref="EvaluationSettings"/> object. These settings objects allow the
    /// caller to specify the maximum number of function evaluations, the termination criteria
    /// for returning a result, and other evaluation settings.</para>
    /// </remarks>
    [CompilerGenerated]
    internal static class NamespaceDoc {
    }

}
