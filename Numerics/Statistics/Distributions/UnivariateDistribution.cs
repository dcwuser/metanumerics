using System;


namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a probability distribution over a single variable.
    /// </summary>
    public abstract class UnivariateDistribution {

        /// <summary>
        /// Gets the mean of the distribution.
        /// </summary>
        public virtual double Mean {
            get {
                return (Moment(1));
            }
        }

        /// <summary>
        /// Gets the variance of the distribution.
        /// </summary>
        public virtual double Variance {
            get {
                return (MomentAboutMean(2));
            }
        }

        /// <summary>
        /// Gets the standard deviation of the distribution.
        /// </summary>
        public virtual double StandardDeviation {
            get {
                return (Math.Sqrt(Variance));
            }
        }

        /// <summary>
        /// Gets the skewness of the distribution.
        /// </summary>
        /// <remarks>
        /// <para>The skewness of a distribution is a measurement of its asymmetry about its mean.
        /// It is the third moment about the mean, measured in units of the cubed standard deviation.</para>
        /// </remarks>
        public virtual double Skewness {
            get {
                return (MomentAboutMean(3) / Math.Pow(MomentAboutMean(2), 3.0 / 2.0));
            }
        }

        /// <summary>
        /// Computes a raw moment of the distribution.
        /// </summary>
        /// <param name="r">The order of the moment to compute.</param>
        /// <returns>The rth raw moment of the distribution.</returns>
        /// <seealso cref="Moment" />
        public abstract double Moment (int r);

        // Override this if sets of raw moments are needed and there is a more efficient way to compute a set of them than by
        // computing each individually. This can occur, for example, if they obey a recursion relation.

        internal virtual double[] RawMoments (int rMax) {
            double[] M = new double[rMax + 1];
            for (int r = 0; r < M.Length; r++) M[r] = Moment(r);
            return (M);
        }

        /// <summary>
        /// Computes a central moment of the distribution.
        /// </summary>
        /// <param name="r">The order of the moment to compute.</param>
        /// <returns>The rth central moment of the distribution.</returns>
        /// <seealso cref="Moment" />
        public virtual double MomentAboutMean (int r) {
            // This is a terrible way to compute central moments, subject to significant cancelations, so replace it if at all possible.
            double[] M = RawMoments(r);
            return(MomentMath.RawToCentral(M, r));
        }

        // Override this if sets of central moments are needed and there is a more efficient way to compute a set of them than by
        // computing each individually. This can occur, for example, if they obey a recursion relation.

        internal virtual double[] CentralMoments (int rMax) {
            double[] C = new double[rMax + 1];
            for (int r = 0; r < C.Length; r++) C[r] = MomentAboutMean(r);
            return (C);
        }

        /// <summary>
        /// Computes a cumulant of the distribution.
        /// </summary>
        /// <param name="r">The index of the cumulant to compute.</param>
        /// <returns>The rth cumulant of the distribution.</returns>
        /// <seealso href="http://en.wikipedia.org/wiki/Cumulant"/>
        public virtual double Cumulant (int r) {
            double[] C = CentralMoments(r);
            double[] K = MomentMath.CentralToCumulant(Mean, C);
            return (K[r]);
        }

        // Override this if sets of cumulants are needed and there is a more efficient way to compute a set of them than by
        // computing each individually. This can occur, for examplek if they obey a recursion relation.

        internal virtual double[] Cumulants (int rMax) {
            double[] K = new double[rMax + 1];
            for (int r = 0; r < K.Length; r++) K[r] = Cumulant(r);
            return (K);
        }

    }

}
