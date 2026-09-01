# Arctan-Skewing: Flexible Skewness for Circular and Toroidal Data

A comprehensive R implementation of the arctan skewing transformation for circular and toroidal distributions. This mechanism enables flexible modeling of skewness in multivariate circular data while maintaining important distributional properties.

## Overview

The arctan skewing mechanism for toroidal data $\boldsymbol{\theta} \in [-\pi, \pi)^d$ is given by

$$f(\boldsymbol{\theta} - \boldsymbol{\mu}; \boldsymbol{\lambda}) := f_0 (\boldsymbol{\theta} - \boldsymbol{\mu}) \left( 1 + \frac{2}{\pi d} \sum_{i=1}^d \arctan \left( \frac{\lambda_i \sin (\theta_i - \mu_i)}{1 + \sin^2 (\theta_i - \mu_i)} \right) \right)$$

for symmetric densities $f_0 (\boldsymbol{\theta} - \boldsymbol{\mu})$, $\boldsymbol{\mu} \in [-\pi, \pi)^d$ and $\boldsymbol{\lambda} \in \mathbb{R}^d$.

This repository contains:

- **Density estimation** for skewed circular and toroidal distributions
- **Random sample generation** from arctan-skewed models
- **Maximum likelihood estimation** tools for parameter inference
- **Skewness and dependence measures** for bivariate distributions
- **Symmetry tests** for circular and toroidal data

## Project Structure

```
Rcode/
├── basic_functions.R                        # Core functions for density and sampling
├── mle.R                                    # Maximum likelihood estimation
├── skewness_dependence_measures.R          # Skewness and dependence calculations
└── symmetry_tests_arctan_skewing.R         # Symmetry tests for circular data
```

## Features

### Supported Distributions

**Univariate (d=1):**
- Von Mises
- Wrapped Cauchy
- Cardioid
- Wrapped Normal
- Kato-Jones

**Bivariate (d=2):**
- Sine
- Cosine
- Bivariate Wrapped Cauchy (BWC)

**Trivariate (d=3):**
- Trivariate Wrapped Cauchy Copula (TWCC)

### Key Functions

- `d_arctan()` - Density estimation with arctan skewing
- `r_arctan()` - Generate random samples from skewed distributions
- `mle_arctan()` - Main wrapper for maximum likelihood estimation (automatically detects dimension)
- `mle_d1()`, `mle_d2()`, `mle.copula_arctan_skewing()` - Dimension-specific MLE functions
- `plot_arctan_contour()` - Visualization of bivariate densities
- Skewness and dependence measure calculations

## Installation

```r
# Clone the repository
git clone https://github.com/Sophia-Loizidou/Arctan-skewing.git

# Load functions in R
source('Rcode/basic_functions.R')
source('Rcode/mle.R')
source('Rcode/skewness_dependence_measures.R')
source('Rcode/symmetry_tests_arctan_skewing.R')
```

### Dependencies

- `circular` - Circular statistics
- `Rsolnp` - Non-linear optimization
- `pracma` - Practical mathematical functions

## Usage Examples

### Univariate Sampling and Density

```r
# Generate samples from skewed von Mises distribution
samples <- r_arctan(n = 1000, d = 1, lambda = c(0.5), 
                    model = 'vonmises', 
                    params = list(mu = 0, kappa = 3))

# Evaluate density
theta <- seq(-pi, pi, length.out = 100)
density <- d_arctan(theta, lambda = c(0.5), model = 'vonmises',
                    params = list(mu = 0, kappa = 3))
```

### Bivariate Visualization

```r
# Plot contours of bivariate skewed sine distribution
plot_arctan_contour(lambda = c(1, 0.5), 
                   model = 'sine',
                   params = list(mu = c(0, 0), kappa = c(1, 1), rho = 0.3))
```

### Maximum Likelihood Estimation

```r
# Univariate: Fit von Mises model with arctan skewing
theta_1d <- r_arctan(n = 200, d = 1, lambda = c(0.5), model = 'vonmises',
                     params = list(mu = 0, kappa = 2))
result_1d <- mle_arctan(theta_1d, symmetric = FALSE, model = 'vonmises')

# Bivariate: Fit sine model with arctan skewing
theta_2d <- r_arctan(n = 200, d = 2, lambda = c(0.3, 0.6), model = 'sine',
                     params = list(mu = c(0, 0), kappa = c(2, 2), rho = 0.2))
result_2d <- mle_arctan(theta_2d, symmetric = FALSE, model = 'sine')

# Trivariate: Fit TWCC model with arctan skewing
theta_3d <- r_arctan(n = 200, d = 3, lambda = c(0.2, 0.3, 0.4), model = 'twcc',
                     params = list(rho12 = 0.5, rho13 = -0.2, rho23 = 0.3))
result_3d <- mle_arctan(theta_3d, symmetric = FALSE, model = 'twcc')
```

## Related Projects

This work builds upon and references:

- **[Symmetry Test on Hypertorus](https://github.com/Sophia-Loizidou/Symmetry_test_on_hypertorus)** - Symmetry testing methods for toroidal data
- **[Trivariate Wrapped Cauchy Copula](https://github.com/Sophia-Loizidou/Trivariate-wrapped-Cauchy-copula)** - Copula methods for multivariate circular data

## License

This project is licensed under the **Attribution 4.0 International (CC-BY-4.0)** license. See [license.md](license.md) for details.

## Author

**Sophia Loizidou**

Department of Mathematics  
University of Luxembourg

## Citation

If you use this code in your research, please cite:

```bibtex
@software{loizidou2026arctan,
  title={Arctan-Skewing: Flexible Skewness for Circular and Toroidal Data},
  author={Loizidou, Sophia},
  year={2026},
  url={https://github.com/Sophia-Loizidou/Arctan-skewing}
}
```

## Questions or Issues?

For questions about the code or to report issues, please open an issue on GitHub or contact the author directly.
