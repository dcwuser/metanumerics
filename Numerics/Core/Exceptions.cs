using System;
using System.Collections.Generic;

namespace Meta.Numerics {

    /// <summary>
    /// The exception that is thrown when an algorithm fails to converge.
    /// </summary>
    public class NonconvergenceException : Exception {

        /// <summary>
        /// Initializes a new non-convergence exception.
        /// </summary>
        public NonconvergenceException () : base("The algorithm did not converge within the allowed number of iterations.") { }

        /// <summary>
        /// Initializes a new non-convergence exception with the given exception message.
        /// </summary>
        /// <param name="message">The exception message.</param>
        public NonconvergenceException (String message) : base(message) { }

        /// <summary>
        /// Initializes a new non-convergence exception with the given exception message and inner exception.
        /// </summary>
        /// <param name="message">The exception message.</param>
        /// <param name="innerException">The inner exception.</param>
        public NonconvergenceException (String message, Exception innerException) : base(message, innerException) { }

    }

    /// <summary>
    /// The exception that is thrown when attempting an operation on objects with incompatible dimensions.
    /// </summary>
    public class DimensionMismatchException : InvalidOperationException {

        /// <summary>
        /// Initializes a new dimension mismatch exception.
        /// </summary>
        public DimensionMismatchException () : base("The object(s) did not have the expected dimension(s).") { }

        /// <summary>
        /// Initializes a new dimension mismatch exception with the given exception message.
        /// </summary>
        /// <param name="message">The exception message.</param>
        public DimensionMismatchException (String message) : base(message) { }

        /// <summary>
        /// Initializes a new dimension mismatch exception with the given exception message and inner exception.
        /// </summary>
        /// <param name="message">The exception message.</param>
        /// <param name="innerException">The inner exception.</param>
        public DimensionMismatchException (String message, Exception innerException) : base(message, innerException) { }

    }

}
