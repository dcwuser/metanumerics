using System;
using System.Collections.Generic;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Contains the result of a fit of time series data to an AR(1) model.
    /// </summary>
    public sealed class AR1FitResult : FitResult {

        internal AR1FitResult(UncertainValue mu, UncertainValue alpha, UncertainValue sigma, IReadOnlyList<double> residuals) : base() {
            this.mu = mu;
            this.alpha = alpha;
            this.sigma = sigma;
            this.residuals = residuals;
            goodnessOfFit = new Lazy<TestResult>(() => residuals.LjungBoxTest());
        }

        private readonly UncertainValue mu, alpha, sigma;
        private readonly IReadOnlyList<double> residuals;
        private readonly Lazy<TestResult> goodnessOfFit;

        /// <summary>
        /// Gets the &#x03bc;-parameter of the model.
        /// </summary>
        /// <value>An estimate, with uncertainty, of the &#x03bc;-parameter of the model, which is the
        /// mean of the underlying population from which the series was drawn.</value>
        public UncertainValue Mu {
            get {
                return (mu);
            }
        }

        /// <summary>
        /// Gets the &#x03b1;-parameter of the model.
        /// </summary>
        /// <value>An estimate, with uncertainty, of the &#x03b1;-parameter of the model, which is the
        /// lag-1 auto-correlation of the underlying population from which the series was drawn.</value>
        public UncertainValue Alpha {
            get {
                return (alpha);
            }
        }

        /// <summary>
        /// Gets the &#x03c3;-parameter of the model.
        /// </summary>
        /// <value>An estimate, with uncertainty, of the &#x03c3;-parameter of the model, which is the
        /// per-step Gaussian noise in the underlying population from which the series was drawn.</value>
        public UncertainValue Sigma {
            get {
                return (sigma);
            }
        }

        /// <summary>
        /// Gets the residuals of the time series fit.
        /// </summary>
        public IReadOnlyList<double> Residuals {
            get {
                return (residuals);
            }
        }

        /// <summary>
        /// Tests the goodness of the AR(1) fit.
        /// </summary>
        /// <returns>A Ljung-Box test for non-correlation of the residuals.</returns>
        public TestResult GoodnessOfFit {
            get {
                return (goodnessOfFit.Value);
            }
        }

        internal override ParameterCollection CreateParameters () {
            string[] names = new string[] { nameof(Mu), nameof(Alpha), nameof(Sigma) };
            ColumnVector b = new ColumnVector(mu.Value, alpha.Value, sigma.Value);
            SymmetricMatrix C = new SymmetricMatrix(3);
            C[0, 0] = MoreMath.Sqr(mu.Uncertainty);
            C[1, 1] = MoreMath.Sqr(alpha.Uncertainty);
            C[2, 2] = MoreMath.Sqr(sigma.Uncertainty);
            // The parameters co-variances vanish in the AR(1) fit.
            return (new ParameterCollection(names, b, C));
        }

    }
}
