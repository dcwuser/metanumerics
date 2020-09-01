using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {
    public static partial class Univariate {

        /// <summary>
        /// Finds the Beta distribution that best fits the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The best fit parameters.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than three values.</exception>
        /// <exception cref="InvalidOperationException">Not all the entries in <paramref name="sample" /> lie between zero and one.</exception>
        public static BetaFitResult FitToBeta (this IReadOnlyList<double> sample) {
            if (sample is null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 3) throw new InsufficientDataException();

            // maximum likelihood calculation
            //   \log L = \sum_i \left[ (\alpha-1) \log x_i + (\beta-1) \log (1-x_i) - \log B(\alpha,\beta) \right]
            // using \frac{\partial B(a,b)}{\partial a} = \psi(a) - \psi(a+b), we have
            //   \frac{\partial \log L}{\partial \alpha} = \sum_i \log x_i -     N \left[ \psi(\alpha) - \psi(\alpha+\beta) \right]
            //   \frac{\partial \log L}{\partial \beta}  = \sum_i \log (1-x_i) - N \left[ \psi(\beta)  - \psi(\alpha+\beta) \right]
            // set equal to zero to get equations for \alpha, \beta
            //   \psi(\alpha) - \psi(\alpha+\beta) = <\log x>
            //   \psi(\beta) - \psi(\alpha+\beta) = <\log (1-x)>

            // compute the mean log of x and (1-x)
            // these are the (logs of) the geometric means
            double ga = 0.0; double gb = 0.0;
            foreach (double value in sample) {
                if ((value <= 0.0) || (value >= 1.0)) throw new InvalidOperationException();
                ga += Math.Log(value); gb += Math.Log(1.0 - value);
            }
            ga /= sample.Count; gb /= sample.Count;

            // define the function to zero
            Func<IReadOnlyList<double>, IReadOnlyList<double>> f = delegate (IReadOnlyList<double> x) {
                double pab = AdvancedMath.Psi(x[0] + x[1]);
                return (new double[] {
                    AdvancedMath.Psi(x[0]) - pab - ga,
                    AdvancedMath.Psi(x[1]) - pab - gb
                });
            };

            // guess initial values using the method of moments
            //   M1 = \frac{\alpha}{\alpha+\beta} C2 = \frac{\alpha\beta}{(\alpha+\beta)^2 (\alpha+\beta+1)}
            // implies
            //   \alpha = M1 \left( \frac{M1 (1-M1)}{C2} - 1 \right)
            //   \beta = (1 - M1) \left( \frac{M1 (1-M1)}{C2} -1 \right)
            ComputeMomentsUpToSecond(sample, out int n, out double m, out double v);
            v /= n;

            double mm = 1.0 - m;
            double q = m * mm / v - 1.0;
            double[] x0 = new double[] { m * q, mm * q };

            // find the parameter values that zero the two equations
            ColumnVector ab = MultiFunctionMath.FindZero(f, x0);
            double a = ab[0]; double b = ab[1];

            // take more derivatives of \log L to get curvature matrix
            //   \frac{\partial^2 \log L}{\partial\alpha^2} = - N \left[ \psi'(\alpha) - \psi'(\alpha+\beta) \right]
            //   \frac{\partial^2 \log L}{\partial\beta^2}  = - N \left[ \psi'(\beta)  - \psi'(\alpha+\beta) \right]
            //   \frac{\partial^2 \log L}{\partial \alpha \partial \beta} = - N \psi'(\alpha+\beta)
            // covariance matrix is inverse of curvature matrix
            SymmetricMatrix C = new SymmetricMatrix(2);
            C[0, 0] = sample.Count * (AdvancedMath.Psi(1, a) - AdvancedMath.Psi(1, a + b));
            C[1, 1] = sample.Count * (AdvancedMath.Psi(1, b) - AdvancedMath.Psi(1, a + b));
            C[0, 1] = sample.Count * AdvancedMath.Psi(1, a + b);
            CholeskyDecomposition CD = C.CholeskyDecomposition();
            if (CD == null) throw new DivideByZeroException();
            C = CD.Inverse();

            // do a KS test on the result
            BetaDistribution distribution = new BetaDistribution(a, b);
            TestResult test = sample.KolmogorovSmirnovTest(distribution);

            return new BetaFitResult(ab, C, distribution, test);
        }

        /// <summary>
        /// Finds the exponential distribution that best fits the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The result of the fit.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than two values.</exception>
        /// <exception cref="InvalidOperationException"><paramref name="sample"/> contains non-positive values.</exception>
        public static ExponentialFitResult FitToExponential (this IReadOnlyList<double> sample) {

            if (sample is null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 2) throw new InsufficientDataException();

            // None of the data is allowed to be negative.
            foreach (double value in sample) {
                if (value < 0.0) throw new InvalidOperationException();
            }

            // It's easy to show that the MLE estimator of \mu is the sample mean and that its variance 
            // is \mu^2 / n, which is just the the variance of the mean, since the variance of the individual
            // values is \mu^2.

            // We can do better than an asymptotic result, though. Since we know that the sum
            // of exponential-distributed values is Gamma-distributed, we know the exact
            // distribution of the mean is Gamma(n, \mu / n). This has mean \mu and variance
            // \mu^2 / n, so the asymptotic results are actually exact.

            double lambda = sample.Mean();
            double dLambda = lambda / Math.Sqrt(sample.Count);

            ContinuousDistribution distribution = new ExponentialDistribution(lambda);
            TestResult test = sample.KolmogorovSmirnovTest(distribution);

            return new ExponentialFitResult(new UncertainValue(lambda, dLambda), test);

        }

        /// <summary>
        /// Finds the Gamma distribution that best fits the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The best fit parameters.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <exception cref="InvalidOperationException"><paramref name="sample"/> contains non-positive values.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than three values.</exception>
        public static GammaFitResult FitToGamma (this IReadOnlyList<double> sample) {

            if (sample is null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 3) throw new InsufficientDataException();

            // The log likelihood of a sample given k and s is
            //   \log L = (k-1) \sum_i \log x_i - \frac{1}{s} \sum_i x_i - N \log \Gamma(k) - N k \log s
            // Differentiating,
            //   \frac{\partial \log L}{\partial s} = \frac{1}{s^2} \sum_i x_i - \frac{N k}{s}
            //   \frac{\partial \log L}{\partial k} = \sum_i \log x_i - N \psi(k) - N \log s
            // Setting the first equal to zero gives
            //   k s = N^{-1} \sum_i x_i = <x>
            //   \psi(k) + \log s = N^{-1} \sum_i \log x_i = <log x>
            // Inserting the first into the second gives a single equation for k
            //   \log k - \psi(k) = \log <x> - <\log x>
            // Note the RHS need only be computed once.
            // \log k > \psi(k) for all k, so the RHS had better be positive. They get
            // closer for large k, so smaller RHS will produce a larger k.

            ComputeMomentsUpToSecond(sample, out int n, out double m, out double v);
            v /= n;

            double s = 0.0;
            foreach (double x in sample) {
                if (x <= 0.0) throw new InvalidOperationException();
                s += Math.Log(x);
            }
            s = Math.Log(m) - s / n;

            // We can get an initial guess for k from the method of moments
            //   \frac{\mu^2}{\sigma^2} = k

            double k0 = MoreMath.Sqr(m) / v;

            // Since 1/(2k) < \log(k) - \psi(k) < 1/k, we could get a bound; that
            // might be better to avoid the solver running into k < 0 territory

            double k1 = FunctionMath.FindZero(k => (Math.Log(k) - AdvancedMath.Psi(k) - s), k0);

            double s1 = m / k1;

            // Curvature of the log likelihood is straightforward
            //   \frac{\partial^2 \log L}{\partial s^2} = -\frac{2}{s^3} \sum_i x_i + \frac{Nk}{s^2} = - \frac{Nk}{s^2}
            //   \frac{\partial^2 \log L}{\partial k \partial s} = - \frac{N}{s}
            //   \frac{\partial^2 \log L}{\partial k^2} = - N \psi'(k)
            // This gives the curvature matrix and thus via inversion the covariance matrix.

            SymmetricMatrix C = new SymmetricMatrix(2);
            C[0, 0] = n * AdvancedMath.Psi(1, k1);
            C[0, 1] = n / s1;
            C[1, 1] = n * k1 / MoreMath.Sqr(s1);
            CholeskyDecomposition CD = C.CholeskyDecomposition();
            if (CD == null) throw new DivideByZeroException();
            C = CD.Inverse();

            // Do a KS test for goodness-of-fit
            GammaDistribution distribution = new GammaDistribution(k1, s1);
            TestResult test = sample.KolmogorovSmirnovTest(distribution);

            return new GammaFitResult(k1, s1, C, distribution, test);
        }

        /// <summary>
        /// Find the Gumbel distribution that best fit the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The fit result.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than three values.</exception>
        public static GumbelFitResult FitToGumbel (this IReadOnlyList<double> sample) {
            if (sample is null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 3) throw new InsufficientDataException();

            // To do a maximum likelihood fit, start from the log probability of each data point and aggregate to
            // obtain the log likelihood of the sample
            //   z_i = \frac{x_i - m}{s}
            //   -\ln p_i = \ln s + ( z_i + e^{-z_i})
            //   \ln L = \sum_i \ln p_i

            // Take derivatives wrt m and s.
            //   \frac{\partial \ln L}{\partial m} = \frac{1}{s} \sum_i ( 1 - e^{-z_i} )
            //   \frac{\partial \ln L}{\partial s} = \frac{1}{s} \sum_i ( -1 + z_i - z_i e^{-z_i} )

            // Set derivatives to zero to get a system of equations for the maximum.
            //    n = \sum_i e^{-z_i}
            //    n = \sum_i ( z_i - z_i e^{-z_i} )
            // that is, <e^z> = 1 and <z> - <z e^z> = 1.

            // To solve this system, pull e^{m/s} out of the sum in the first equation and solve for m
            //    n = e^{m / s} \sum_i e^{-x_i / s}
            //    m = -s \ln \left( \frac{1}{n} \sum_i e^{-x_i / s} \right) = -s \ln <e^{-x/s}>
            // Substituting this result into the second equation gets us to
            //    s = \bar{x} - \frac{ <x e^{-x/s}> }{ <e^{x/s}> }  
            // which involves only s. We can use a one-dimensional root-finder to determine s, then determine m
            // from the first equation.

            // To avoid exponentiating potentially large x_i, it's better to write the problem in terms
            // of d_i, where x_i = \bar{x} + d_i.
            //    m = \bar{x} - s \ln <e^{-d/s}>
            //    s = -\frac{ <d e^{-d/s}> }{ <e^{-d/s}> }

            // To get the covariance matrix, we need the curvature matrix at the minimum, so take more derivatives
            //    \frac{\partial^2 \ln L}{\partial m^2} = - \frac{1}{s} \sum_i e^{-z_i} = - \frac{n}{s^2}
            //    \frac{\partial^2 \ln L}{\partial m \partial s} = - \frac{n}{s^2} <z e^{-z}>
            //    \frac{\partial^2 \ln L}{\partial s^2} = - \frac{n}{s^2} ( <z^2 e^{-z}> + 1 )

            // Several crucial pieces of this analysis are taken from Mahdi and Cenac, "Estimating Parameters of Gumbel Distribution
            // "using the method of moments, probability weighted moments, and maximum likelihood", Revista de Mathematica:
            // Teoria y Aplicaciones 12 (2005) 151-156 (http://revistas.ucr.ac.cr/index.php/matematica/article/viewFile/259/239) 

            // We will be needed the sample mean and standard deviation
            Univariate.ComputeMomentsUpToSecond(sample, out int n, out double mean, out double stdDev);
            stdDev = Math.Sqrt(stdDev / n);

            // Use the method of moments to get an initial estimate of s.
            double s0 = Math.Sqrt(6.0) / Math.PI * stdDev;

            // Define the function to zero
            Func<double, double> fnc = (double s) => {
                double u, v;
                MaximumLikelihoodHelper(sample, n, mean, s, out u, out v);
                return (s + v / u);
            };

            // Zero it to compute the best-fit s
            double s1 = FunctionMath.FindZero(fnc, s0);

            // Compute the corresponding best-fit m
            double u1, v1;
            MaximumLikelihoodHelper(sample, n, mean, s1, out u1, out v1);
            double m1 = mean - s1 * Math.Log(u1);

            // Compute the curvature matrix
            double w1 = 0.0;
            double w2 = 0.0;
            foreach (double x in sample) {
                double z = (x - m1) / s1;
                double e = Math.Exp(-z);
                w1 += z * e;
                w2 += z * z * e;
            }
            w1 /= sample.Count;
            w2 /= sample.Count;
            SymmetricMatrix C = new SymmetricMatrix(2);
            C[0, 0] = (n - 2) / (s1 * s1);
            C[0, 1] = (n - 2) / (s1 * s1) * w1;
            C[1, 1] = (n - 2) / (s1 * s1) * (w2 + 1.0);
            SymmetricMatrix CI = C.CholeskyDecomposition().Inverse();
            // The use of (n-2) here in place of n is a very ad hoc attempt to increase accuracy.

            // Compute goodness-of-fit
            GumbelDistribution dist = new GumbelDistribution(m1, s1);
            TestResult test = sample.KolmogorovSmirnovTest(dist);

            return new GumbelFitResult(m1, s1, CI[0, 0], CI[1, 1], CI[0, 1], test);
        }

        // Compute u = <e^{-d/s}> and v = <d e^{-d/s}> for a given s and a given sample.

        private static void MaximumLikelihoodHelper (IEnumerable<double> sample, int n, double mean, double s, out double u, out double v) {
            u = 0.0;
            v = 0.0;
            foreach (double x in sample) {
                double d = x - mean;
                double e = Math.Exp(-d / s);
                u += e;
                v += d * e;
            }
            u /= n;
            v /= n;
        }

        /// <summary>
        /// Finds the log-normal distribution that best fits the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The best fit parameters.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than three values.</exception>
        /// <exception cref="InvalidOperationException"><paramref name="sample"/> contains non-positive values.</exception>
        public static LognormalFitResult FitToLognormal (this IReadOnlyList<double> sample) {
            if (sample is null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 3) throw new InsufficientDataException();

            // Writing out the log likelihood from p(x), taking its derivatives wrt mu and sigma, and setting them equal
            // to zero to find the minimizing values, you find that the results of the normal fit are reproduced exactly
            // with x -> log x, i.e. \mu = < \log x >, \sigma^2 = < (\log x - \mu)^2 >. So we just repeat the normal fit
            // logic with x -> log x.

            FitToNormalInternal(sample.Select(x => Math.Log(x)), out double m, out double dm, out double s, out double ds);

            LognormalDistribution distribution = new LognormalDistribution(m, s);
            TestResult test = sample.KolmogorovSmirnovTest(distribution);
            return new LognormalFitResult(new UncertainValue(m, dm), new UncertainValue(s, ds), distribution, test);

        }

        /// <summary>
        /// Finds the normal distribution that best fits the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The result of the fit.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than three values.</exception>
        public static NormalFitResult FitToNormal (this IReadOnlyList<double> sample) {

            if (sample is null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 3) throw new InsufficientDataException();

            FitToNormalInternal(sample, out double m, out double dm, out double s, out double ds);

            NormalDistribution distribution = new NormalDistribution(m, s);
            TestResult test = sample.KolmogorovSmirnovTest(distribution);

            return new NormalFitResult(new UncertainValue(m, dm), new UncertainValue(s, ds), distribution, test);

        }

        private static void FitToNormalInternal (IEnumerable<double> sample, out double m, out double dm, out double s, out double ds) {

            Debug.Assert(sample != null);

            // We factor out this method because it is used by both the normal and the log-normal fit methods.

            // Maximum likelihood estimate is straightforward.
            //   p_i = \frac{1}{\sqrt{2\pi}\sigma} \exp \left[ -\frac{1}{2} \left( \frac{x_i - \mu}{\sigma} \right)^2 \right]
            //   \ln p_i = -\ln (\sqrt{2\pi} \sigma) - \frac{1}{2} \left( \frac{x_i - \mu}{\sigma} \right)^2
            //   \ln L = \sum_i
            //  so
            //    \frac{\partial \ln L}{\partial \mu} = \sum_i \frac{x_i - \mu}{\sigma^2}
            //    \frac{\partial \ln L}{\partial \sigma} = -\frac{n}{\sigma} - \frac{1}{\sigma^2} \sum_i (x_i - \mu)^2
            //  Setting equal to zero and solving gives the unsurprising result
            //     \mu = n^{-1} \sum_i x_i
            //     \sigma^2 = n^{-1} \sum_i (x_i - \mu)^2
            //  that MLE says to estimate the model mean and variance by the sample mean and variance.

            // MLE estimators are guaranteed to be asymptotically unbiased, but they can be biased for finite n.
            // You can see that must be the case for \sigma because the denominator has n instead of n-1.

            // To un-bias our estimators, we will derive exact distributions for these quantities.

            // First the mean estimator. Start from x_i \sim N(\mu, \sigma). By the addition of normal deviates,
            //   \sum_i x_i \sim N(n \mu, \sqrt{n} \sigma). So
            //   m  = \frac{1}{n} \sum_i x_i \sim N(\mu, \sigma / \sqrt{n}).
            // which means the estimator m is normally distributed with mean \mu and standard deviation
            // \sigma / \sqrt{n}. Now we know that m is unbiased and we know its variance.

            Univariate.ComputeMomentsUpToSecond(sample, out int n, out m, out double ss);
            dm = Math.Sqrt(ss) / n;

            // Next the variance estimator. By the definition of the chi squared distribution and a bit of algebra that
            // reduces the degrees of freedom by one, u^2 = \sum_i ( \frac{x_i - m}{\sigma} )^2 \sim \chi^2(n - 1), which has
            // mean n - 1 and variance 2(n-1). Therefore the estimator
            //   v = \sigma^2 u^2 / (n-1) = \frac{1}{n-1} \sum_i ( x_i - m )^2
            // has mean \sigma^2 and variance 2 \sigma^4 / (n-1). 

            // If we consider \sigma^2 the parameter, we are done -- we have derived an estimator that is unbiased and
            // know its variance. But we don't consider the parameter \sigma^2, we consider it \sigma.
            // The mean of the square root is not the square root of the mean, so the square root of an unbiased
            // estimator of \sigma^2 will not be an unbiased estimator of \sigma. If we want an unbiased estimator of \sigma
            // itself, we need to go a bit further. Since u^2 ~ \chi^2(n-1), u ~ \chi(n-1). Its mean is a complicated ratio
            // of Gamma functions and it's variance is an even more complicated difference whose evaluation can be delicate,
            // but our machinery in the ChiDistribution class handles that. To get an unbiased estimator of \sigma, we just
            // need to apply the same principal of dividing by the mean of this distribution.
            //   s = \sigma u / <u> = \sqrt{\sum_i (x_i - m)^2} / <u>
            // to get an estimator with mean \sigma and known variance.

            ChiDistribution d = new ChiDistribution(n - 1);
            s = Math.Sqrt(ss) / d.Mean;
            ds = d.StandardDeviation / d.Mean * s;

        }

        /// <summary>
        /// Finds the Rayleigh distribution that best fits the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit, which must have at least 2 values.</param>
        /// <returns>The fit result.</returns>
        public static RayleighFitResult FitToRayleigh (this IReadOnlyList<double> sample) {

            if (sample is null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 2) throw new InsufficientDataException();

            // It's easy to follow maximum likelihood prescription, because there is only one
            // parameter and the functional form is fairly simple.

            // \ln L = \sum_i p_i = = \sum_i \left[ \ln x_i - 2 \ln \sigma - \frac{1}{2} \frac{x_i^2}{\sigma^2 \right]
            // \frac{\partial \ln L}{\partial \sigma} = \sum_i \left[ \frac{x_i^2}{\sigma^3} - \frac{2}{\sigma} \right]
            // \frac{\partial^2 \ln L}{\partial \sigma^2} = \sum_i \left[ \frac{2}{\sigma^2} - \frac{3 x_i^2}{\sigma^4} \right]

            // Set the first derivative to zero to obtain
            //   \hat{\sigma}^2 = \frac{1}{2n} \sum_i x_i^2
            //   \hat{\sigma} = \sqrt{\frac{1}{2} \left< x^2 \right>}
            // and plug this value into the second derivative to obtain
            //   \frac{\partial^2 \ln L}{\partial \sigma^2} = - \frac{4n}{\sigma^2}
            // at the minimum, so
            //   \delta \sigma = \frac{\sigma}{\sqrt{4n}}

            // Next consider bias. We know from the moments of the distribution that < x^2 > = 2 \sigma^2, so \hat{\sigma}^2 is
            // unbiased. But its square root \hat{\sigma} is not. To get exact distribution of \hat{\sigma},
            //    x_i \sim Rayleigh(\sigma)
            //    ( \frac{x_i}{\sigma} )^2 \sim \chi^2(2)
            //    \sum_i ( \frac{x_i}{\sigma} )^2 \sim \chi^2(2n)
            //    \left[ \sum_i ( \frac{x_i}{\sigma} )^2 \right]^{1/2} \sim \chi(2n)
            // And if z \sim \chi(k) then
            //    E(z) = \sqrt{2} \frac{\Gamma((k + 1)/2)}{\Gamma(k / 2)}
            //    V(z) = k - [ E(z) ]^2
            // Here k = 2n and our estimator \hat{\sigma} = z / sqrt{2n}, so
            //    E(\hat{\sigma}) = \frac{\Gamma(n + 1/2)}{\sqrt{n} \Gamma(n)} \sigma = \frac{(n)_{1/2}}{\sqrt{n}} \sigma
            //    V(\hat{\sigma}) = \left[ 1 - \frac{(n)_{1/2}^2}{n} \right] \sigma^2
            // We can use series expansion to verify that
            //    E(\hat{\sigma}) = \left[ 1 - \frac{1}{8n} + \cdots \right] \sigma
            //    V(\hat{\sigma}) = \left[ \frac{1}{4n} - \frac{1}{32n^2} + \cdots \right] \sigma^2
            // our estimator is asymptotically unbiased and its asymptotic variance agrees with our assessment.

            // We correct for the bias of \hat{\sigma} by multiplying by \frac{\sqrt{n}}{(n)_{1/2}}. We could
            // correct our variance estimate too, but in order to evaluate the correction at high n,
            // we need to do a series expansion, and I regard it as less important that the error estimate
            // be exact.

            int n = sample.Count;
            double s = Math.Sqrt(sample.RawMoment(2) / 2.0) * (Math.Sqrt(n) / AdvancedMath.Pochhammer(n, 1.0 / 2.0));
            double ds = s / Math.Sqrt(4.0 * n);

            RayleighDistribution distribution = new RayleighDistribution(s);
            TestResult test = sample.KolmogorovSmirnovTest(distribution);

            return new RayleighFitResult(new UncertainValue(s, ds), distribution, test);
        }

        /// <summary>
        /// Finds the Wald distribution that best fits a sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The fit.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <exception cref="InvalidOperationException"><paramref name="sample"/> contains non-positive values.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than three values.</exception>
        public static WaldFitResult FitToWald (this IReadOnlyList<double> sample) {

            if (sample is null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 3) throw new InsufficientDataException();

            // For maximum likelihood estimation, take logs of pdfs and sum:
            //    \ln p = \frac{1}{2} \ln \lambda - \frac{1}{2} \ln (2\pi) - \frac{3}{2} \ln x
            //            - \frac{\lambda x}{2\mu^2} + \frac{\lambda}{\mu} - \frac{\lambda}{2x}
            //    \ln L = \sum_i p_i

            // Take derivative wrt \mu
            //    \frac{\partial \ln L}{\partial \mu} = \sum_i \left[ \frac{\lambda x_i}{\mu^3} - \frac{\lambda}{\mu^2} \right]
            // and set equal to zero to obtain
            //    \mu = \frac{1}{n} \sum_i x_i = <x>
            // which agrees with method of moments.

            // Take derivative wrt \lambda
            //    \frac{\partial \ln L}{\partial \lambda} = \sum_i \left[ \frac{1}{2 \lambda} -\frac{x_i}{2\mu^2} + \frac{1}{\mu} - \frac{1}{2 x_i} \right]
            // Set equal to zero, plug in our expression for \mu, and solve for \lambda to get
            //    \frac{n}{\lambda} = \sum_i \left( \frac{1}{x_i} - \frac{1}{\mu} \right)
            //  i.e. \lambda^{-1} = <(x^{-1} - \mu^{-1})>

            ComputeMomentsUpToFirst(sample, out int n, out double mu);
            double mui = 1.0 / mu;
            double lambda = 0.0;
            foreach (double value in sample) {
                if (value <= 0.0) throw new InvalidOperationException();
                lambda += (1.0 / value - mui);
            }
            lambda = (n - 3) / lambda;

            // If x ~ IG(\mu, \lambda), then \sum_i x_i ~ IG(n \mu, n^2 \lambda), so \hat{\mu} ~ IG (\mu, n \lambda). This gives us
            // not just the exact mean and variance of \hat{\mu}, but its entire distribution. Since its mean is \mu, \hat{\mu} is an
            // unbiased estimator. And by the variance formula for IG, the variance of \hat{\mu} is \frac{\mu^3}{n \lambda}.

            // Tweedie, "Statistical Properties of Inverse Gaussian Distributions" (http://projecteuclid.org/download/pdf_1/euclid.aoms/1177706964)
            // showed that \frac{n \lambda}{\hat{\lambda}} ~ \chi^2_{n-1}. Since the mean of \chi^2_{k} is k, the MLE estimator of
            // \frac{1}{\lambda} can be made unbiased by replacing n by (n-1). However, we are estimating \lambda, not \frac{1}{\lambda}.
            // By the relation between chi squared and inverse chi squared distributions, \frac{\hat{\lambda}}{n \lambda} ~ I\chi^2_{n-1}.
            // The mean of I\chi^2_{k} is \frac{1}{n-2}, so to get an unbiased estimator of \lambda, we need to replace n by (n-3). This is
            // what we have done above. Furthermore, the variance of I\chi^2_{k} is \frac{2}{(k-2)^2 (k-4)}, so the variance of \hat{\lambda}
            // is \frac{2 \lambda^2}{(n-5)}.

            // We can also get covariances from the MLE approach. To get a curvature matrix, take additional derivatives
            //   \frac{\partial^2 \ln L}{\partial \mu^2} = \sum_i \left[ -\frac{3 \lambda x_i}{\mu^4} + \frac{2 \lambda}{\mu^3} \right]
            //   \frac{\partial^2 \ln L}{\partial \mu \partial \lambda} = \sum_i \left[ \frac{x_i}{\mu^3} - \frac{1}{\mu^2} \right]
            //   \frac{\partial^2 \ln L}{\partial \lambda^2} =\sum_i \left[ - \frac{1}{2 \lambda^2} \right]
            // and substitutue in best-fit values of \mu and \lambda
            //   \frac{\partial^2 \ln L}{\partial \mu^2} = - \frac{n \lambda}{\mu^3}
            //   \frac{\partial^2 \ln L}{\partial \mu \partial \lambda} = 0
            //   \frac{\partial^2 \ln L}{\partial \lambda^2} = - \frac{n}{2 \lambda^2}
            // Mixed derivative vanishes, so matrix is trivially invertible to obtain covariances. These results agree with the
            // results from exact distributions in the asymptotic regime.

            double v_mu_mu = mu * mu * mu / lambda / n;
            double v_lambda_lambda = 2.0 * lambda * lambda / (n - 5);
            //double v_mu_lambda = 0.0;

            WaldDistribution dist = new WaldDistribution(mu, lambda);
            TestResult test = sample.KolmogorovSmirnovTest(dist);
            return new WaldFitResult(mu, lambda, v_mu_mu, v_lambda_lambda, dist, test);
        }

        /// <summary>
        /// Finds the Weibull distribution that best fits the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The best fit parameters.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <exception cref="InvalidOperationException"><paramref name="sample"/> contains non-positive values.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than three values.</exception>
        public static WeibullFitResult FitToWeibull (this IReadOnlyList<double> sample) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 3) throw new InsufficientDataException();
            foreach (double value in sample) {
                if (value <= 0.0) throw new InvalidOperationException();
            }

            // The log likelihood function is
            //   \log L = N \log k + (k-1) \sum_i \log x_i - N K \log \lambda - \sum_i \left(\frac{x_i}{\lambda}\right)^k
            // Taking derivatives, we get
            //   \frac{\partial \log L}{\partial \lambda} = - \frac{N k}{\lambda} + \sum_i \frac{k}{\lambda} \left(\frac{x_i}{\lambda}\right)^k
            //   \frac{\partial \log L}{\partial k} =\frac{N}{k} + \sum_i \left[ 1 - \left(\frac{x_i}{\lambda}\right)^k \right] \log \left(\frac{x_i}{\lambda}\right)
            // Setting the first expression to zero and solving for \lambda gives
            //   \lambda = \left( N^{-1} \sum_i x_i^k \right)^{1/k} = ( < x^k > )^{1/k}
            // which allows us to reduce the problem from 2D to 1D.
            // By the way, using the expression for the moment < x^k > of the Weibull distribution, you can show there is
            // no bias to this result even for finite samples.
            // Setting the second expression to zero gives
            //   \frac{1}{k} = \frac{1}{N} \sum_i \left[ \left( \frac{x_i}{\lambda} \right)^k - 1 \right] \log \left(\frac{x_i}{\lambda}\right)
            // which, given the equation for \lambda as a function of k derived from the first expression, is an implicit equation for k.
            // It cannot be solved in closed form, but we have now reduced our problem to finding a root in one-dimension.

            // We need a starting guess for k.
            // The method of moments equations are not solvable for the parameters in closed form
            // but the scale parameter drops out of the ratio of the 1/3 and 2/3 quantile points
            // and the result is easily solved for the shape parameter
            //   k = \frac{\log 2}{\log\left(\frac{x_{2/3}}{x_{1/3}}\right)}
            double x1 = sample.InverseLeftProbability(1.0 / 3.0);
            double x2 = sample.InverseLeftProbability(2.0 / 3.0);
            double k0 = Global.LogTwo / Math.Log(x2 / x1);
            // Given the shape paramter, we could invert the expression for the mean to get
            // the scale parameter, but since we have an expression for \lambda from k, we
            // don't need it.
            //double s0 = sample.Mean / AdvancedMath.Gamma(1.0 + 1.0 / k0);

            // Simply handing our 1D function to a root-finder works fine until we start to encounter large k. For large k,
            // even just computing \lambda goes wrong because we are taking x_i^k which overflows. Horst Rinne, "The Weibull
            // Distribution: A Handbook" describes a way out. Basically, we first move to variables z_i = \log(x_i) and
            // then w_i = z_i - \bar{z}. Then lots of factors of e^{k \bar{z}} cancel out and, even though we still do
            // have some e^{k w_i}, the w_i are small and centered around 0 instead of large and centered around \lambda.

            //Sample transformedSample = sample.Copy();
            //transformedSample.Transform(x => Math.Log(x));
            double[] transformedSample = new double[sample.Count];
            for (int j = 0; j < sample.Count; j++) transformedSample[j] = Math.Log(sample[j]);
            double zbar = transformedSample.Mean();
            for (int j = 0; j < transformedSample.Length; j++) transformedSample[j] -= zbar;

            // After this change of variable the 1D function to zero becomes
            //   g(k) = \sum_i ( 1 - k w_i ) e^{k w_i}
            // It's easy to show that g(0) = n and g(\infinity) = -\infinity, so it must cross zero. It's also easy to take
            // a derivative
            //   g'(k) = - k \sum_i w_i^2 e^{k w_i}
            // so we can apply Newton's method.

            int i = 0;
            double k1 = k0;
            while (true) {
                i++;
                double g = 0.0;
                double gp = 0.0;
                foreach (double w in transformedSample) {
                    double e = Math.Exp(k1 * w);
                    g += (1.0 - k1 * w) * e;
                    gp -= k1 * w * w * e;
                }
                double dk = -g / gp;
                k1 += dk;
                if (Math.Abs(dk) <= Global.Accuracy * Math.Abs(k1)) break;
                if (i >= Global.SeriesMax) throw new NonconvergenceException();
            }

            // The corresponding lambda can also be expressed in terms of zbar and w's.

            double t = 0.0;
            foreach (double w in transformedSample) {
                t += Math.Exp(k1 * w);
            }
            t /= transformedSample.Length;
            double lambda1 = Math.Exp(zbar) * Math.Pow(t, 1.0 / k1);

            // We need the curvature matrix at the minimum of our log likelihood function
            // to determine the covariance matrix. Taking more derivatives...
            //    \frac{\partial^2 \log L} = \frac{N k}{\lambda^2} - \sum_i \frac{k(k+1) x_i^k}{\lambda^{k+2}}
            //    = - \frac{N k^2}{\lambda^2}
            // The second expression follows by inserting the first-derivative-equal-zero relation into the first.
            // For k=1, this agrees with the variance formula for the mean of the best-fit exponential.

            // Derivatives involving k are less simple.

            // We end up needing the means < (x/lambda)^k log(x/lambda) > and < (x/lambda)^k log^2(x/lambda) >

            double mpl = 0.0; double mpl2 = 0.0;
            foreach (double x in sample) {
                double r = x / lambda1;
                double p = Math.Pow(r, k1);
                double l = Math.Log(r);
                double pl = p * l;
                double pl2 = pl * l;
                mpl += pl;
                mpl2 += pl2;
            }
            mpl = mpl / sample.Count;
            mpl2 = mpl2 / sample.Count;

            // See if we can't do any better here. Transforming to zbar and w's looked ugly, but perhaps it
            // can be simplified? One interesting observation: if we take expectation values (which gives
            // the Fisher information matrix) the entries become simple:
            //   B_{\lambda \lambda} = \frac{N k^2}{\lambda^2}
            //   B_{\lambda k} = -\Gamma'(2) \frac{N}{\lambda}
            //   B_{k k } = [1 + \Gamma''(2)] \frac{N}{k^2}
            // Would it be bad to just use these directly?

            // Construct the curvature matrix and invert it.
            SymmetricMatrix C = new SymmetricMatrix(2);
            C[0, 0] = sample.Count * MoreMath.Sqr(k1 / lambda1);
            C[0, 1] = -sample.Count * k1 / lambda1 * mpl;
            C[1, 1] = sample.Count * (1.0 / MoreMath.Sqr(k1) + mpl2);
            CholeskyDecomposition CD = C.CholeskyDecomposition();
            if (CD == null) throw new DivideByZeroException();
            C = C.Inverse();

            // Do a KS test to compare sample to best-fit distribution
            WeibullDistribution distribution = new WeibullDistribution(lambda1, k1);
            TestResult test = sample.KolmogorovSmirnovTest(distribution);

            // return the result
            return (new WeibullFitResult(lambda1, k1, C, distribution, test));
        }

        /// <summary>
        /// Finds the parameters that make an arbitrary, parameterized distribution best fit the sample.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <param name="factory">A delegate that creates a distribution given its parameters.</param>
        /// <param name="start">An initial guess at the parameter values.</param>
        /// <returns>The best fit parameters.</returns>
        public static DistributionFitResult<ContinuousDistribution> MaximumLikelihoodFit(this IReadOnlyList<double> sample, Func<IReadOnlyDictionary<string,double>, ContinuousDistribution> factory, IReadOnlyDictionary<string, double> start) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (factory == null) throw new ArgumentNullException(nameof(factory));
            if (start == null) throw new ArgumentNullException(nameof(start));

            List<string> names = start.Keys.ToList();

            // Convert dictionary-based arguments from user to list-based arguments used internally.

            Func<IReadOnlyList<double>, ContinuousDistribution> listFactory = a => {
                Dictionary<string, double> d = new Dictionary<string, double>(names.Count);
                for (int i = 0; i < names.Count; i++) {
                    d.Add(names[i], a[i]);
                }
                return (factory(d));
            };

            List<double> listStart = new List<double>(start.Count);
            foreach (string name in names) {
                listStart.Add(start[name]);
            }

            return (MaximumLikelihoodFit(sample, listFactory, listStart, names));
        }

        internal static DistributionFitResult<ContinuousDistribution> MaximumLikelihoodFit (IReadOnlyList<double> sample, Func<IReadOnlyList<double>, ContinuousDistribution> factory, IReadOnlyList<double> start, IReadOnlyList<string> names) {

            Debug.Assert(sample != null);
            Debug.Assert(factory != null);
            Debug.Assert(start != null);
            Debug.Assert(names != null);
            Debug.Assert(start.Count == names.Count);

            // Define a log likelihood function
            Func<IReadOnlyList<double>, double> logL = (IReadOnlyList<double> a) => {
                ContinuousDistribution d = factory(a);
                double lnP = 0.0;
                foreach (double value in sample) {
                    double P = d.ProbabilityDensity(value);
                    if (P == 0.0) throw new InvalidOperationException();
                    lnP += Math.Log(P);
                }
                return (lnP);
            };

            // Maximize it
            MultiExtremum maximum = MultiFunctionMath.FindLocalMaximum(logL, start);
            ColumnVector b = maximum.Location;
            SymmetricMatrix C = maximum.HessianMatrix;
            CholeskyDecomposition CD = C.CholeskyDecomposition();
            if (CD == null) throw new DivideByZeroException();
            C = CD.Inverse();

            ContinuousDistribution distribution = factory(maximum.Location);
            TestResult test = sample.KolmogorovSmirnovTest(distribution);

            return (new ContinuousDistributionFitResult(names, b, C, distribution, test));
        }

    }

    internal class ContinuousDistributionFitResult : DistributionFitResult<ContinuousDistribution> {

        public ContinuousDistributionFitResult (IReadOnlyList<string> names, ColumnVector parameters, SymmetricMatrix covariance, ContinuousDistribution distribution, TestResult goodnessOfFit) : base (distribution, goodnessOfFit) {
            this.parameters = new ParameterCollection(names, parameters, covariance);
        }

        private readonly ParameterCollection parameters;

        internal override ParameterCollection CreateParameters () {
            return (parameters);
        }
    }

}
