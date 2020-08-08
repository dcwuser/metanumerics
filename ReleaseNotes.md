# Meta.Numerics 4.1

## More Extended Precision Types
We have made multiple improvements to the  the quadruple-precision (~32 decimal digit) DoubleDouble floating point type. It now handles infinities and NaNs appropriately. Its string rendering and parsing are improved. We have also added more DoubleDouble functions, including trig functions and LogGamma.

We have added the 128-bit integer types Int128 and UInt128. These behave like the built-in integer types Int64 and UInt64 (long and ulong), but support integer values up to ~10<sup>38</sup>. Arithmetic using these types is 1-4 times faster than using BigInteger, and unlike BigInteger, they behave like the other fixed-width register types with respect to overflow.

## More Advanced Functions
We have added a few more advanced functions. These include the complete elliptic integral of the third kind and computation of the elliptic nome (so now only the incomplete elliptic integral of the third kind remains unimplemented). We also added a scaled version of the incomplete Bessel function (allowing you to work with the function for arguments where the function value itself would over- or under-flow), functions that return the zeros of the Airy and Bessel functions, and the hyperbolic integral functions Cin and Shi.

We have also made many improvements to the internals of long-implemented functions to improve their speed, accuracy, and behavior at extreme arguments including infinities and NaN. 

## Other Improvements
We have added the RegressionResult type with exposes residuals and the sum of squared residuals on all FitResults that have them. We fixed a bug which could cause non-convergence in the multi-dimensional FindLocalMaximum and FindLocalMinimum methods. 

## Future Plans
The next release of the Meta.Numerics will use .NET Standard 2.1.
We will eliminate classes like Sample, BivariateSample, and MultivariateSample in favor of the Univariate, Bivariate, and Multivariate extension method classes which can be used with any type implementing IReadOnlyList, including the built-in collection types and our own DataFrame types.