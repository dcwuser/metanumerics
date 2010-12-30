using System;
using System.Runtime.Serialization;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// The exception that is thrown when an operation is attempted with less than the minimum required data.
    /// </summary>
    [Serializable]
    public class InsufficientDataException : InvalidOperationException {

        /// <summary>
        /// Initializes a new insufficient data exception.
        /// </summary>
        public InsufficientDataException () : base() { }

        /// <summary>
        /// Inititalizes a new insufficient data exception with the given exception message.
        /// </summary>
        /// <param name="message">The exception message.</param>
        public InsufficientDataException (String message) : base(message) { }

        /// <summary>
        /// Initializes a new insufficient data exception with the given exception message and inner exception.
        /// </summary>
        /// <param name="message">The exeption message.</param>
        /// <param name="innerException">The inner exception.</param>
        public InsufficientDataException (String message, Exception innerException) : base(message, innerException) { }

        /// <summary>
        /// Initalizes a new insufficient data exception with the given serialization information and streaming context.
        /// </summary>
        /// <param name="info">The serialization information.</param>
        /// <param name="context">The streaming context.</param>
        protected InsufficientDataException (SerializationInfo info, StreamingContext context) : base(info, context) { }


    }

}
