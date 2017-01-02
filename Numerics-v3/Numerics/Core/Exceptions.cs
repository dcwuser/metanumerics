using System;
using System.Collections.Generic;
#if !SILVERLIGHT
using System.Runtime.Serialization;
#endif

namespace Meta.Numerics {

    /// <summary>
    /// The exception that is thrown when an algorithm fails to converge.
    /// </summary>
#if !SILVERLIGHT
    [Serializable]
#endif
    public class NonconvergenceException : Exception {

        /// <summary>
        /// Initializes a new nonconvergence exception.
        /// </summary>
        public NonconvergenceException () : base() { }

        /// <summary>
        /// Inititalizes a new nonconvergence exception with the given exception message.
        /// </summary>
        /// <param name="message">The exception message.</param>
        public NonconvergenceException (String message) : base(message) { }

        /// <summary>
        /// Initializes a new nonconvergence exception with the given exception message and inner exception.
        /// </summary>
        /// <param name="message">The exeption message.</param>
        /// <param name="innerException">The inner exception.</param>
        public NonconvergenceException (String message, Exception innerException) : base(message, innerException) { }

#if !SILVERLIGHT
        /// <summary>
        /// Initalizes a new nonconvergence exception with the given serialization information and streaming context.
        /// </summary>
        /// <param name="info">The serialization information.</param>
        /// <param name="context">The streaming context.</param>
        protected NonconvergenceException (SerializationInfo info, StreamingContext context) : base(info, context) { }
#endif

    }

    /// <summary>
    /// The exception that is thrown when attempting an operation on objects with incompatible dimensions.
    /// </summary>
#if !SILVERLIGHT
    [Serializable]
#endif
    public class DimensionMismatchException : InvalidOperationException {

        /// <summary>
        /// Initializes a new dimension mismatch exception.
        /// </summary>
        public DimensionMismatchException () : base(Messages.DimensionMismatch) { }

        /// <summary>
        /// Inititalizes a new dimension mismatch exception with the given exception message.
        /// </summary>
        /// <param name="message">The exception message.</param>
        public DimensionMismatchException (String message) : base(message) { }

        /// <summary>
        /// Initializes a new dimension mismatch exception with the given exception message and inner exception.
        /// </summary>
        /// <param name="message">The exeption message.</param>
        /// <param name="innerException">The inner exception.</param>
        public DimensionMismatchException (String message, Exception innerException) : base(message, innerException) { }

#if !SILVERLIGHT
        /// <summary>
        /// Initalizes a new dimension mismatch exception with the given serialization information and streaming context.
        /// </summary>
        /// <param name="info">The serialization information.</param>
        /// <param name="context">The streaming context.</param>
        protected DimensionMismatchException (SerializationInfo info, StreamingContext context) : base(info, context) { }
#endif

    }

}
