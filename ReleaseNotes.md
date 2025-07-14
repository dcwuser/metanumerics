# Meta.Numerics 4.2

## More Distributions and Fits
We have added Skellam, Benford and SkewNormal distributions. We have added FitResult classes specific to distributions, with members that make it easier to access the parameters of that distribution.

## More Advanced Functions
We have addede Gegenbauer polynomials and Chebyshev polynomials of the second kind. We have added Complex Airy functions. We have made improvements to the performance and accuracy of the Gamma family of functions (Gamma, Beta, Pochhammer) and added a LogPochhammer function.

## Permutation Matrices
The Permutation class can now produce a PermutationMatrix, which can be used to permute rows and columns of vectors and matrices.

## Future Plans
The next release of the Meta.Numerics will either use .NET Standard 2.1.
We will eliminate classes like Sample, BivariateSample, and MultivariateSample in favor of the Univariate, Bivariate, and Multivariate extension method classes which can be used with any type implementing IReadOnlyList, including the built-in collection types and our own DataFrame types.
Since the AdvancedMath class is getting so big, we will move Bessel-type functions to their own static class.