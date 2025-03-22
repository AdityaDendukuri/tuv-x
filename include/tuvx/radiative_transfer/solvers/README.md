# Delta-Eddington Radiative Transfer Solver

## Overview

This implementation solves the radiative transfer equation using the Delta-Eddington approximation as described in Toon et al. (1989). The implementation converts the Fortran algorithm to C++ and optimizes memory layout for vectorization. It follows the mathematical formulation from the original papers while leveraging tuv-x's existing components like Array3D and TridiagonalSolver.

## References

1. [Toon et al. (1989): "Rapid Calculation of Radiative Heating Rates and Photodissociation Rates in Inhomogeneous Multiple Scattering Atmospheres"](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JD094iD13p16287)
2. [Joseph et al. (1976): "The Delta-Eddington Approximation for Radiative Flux Transfer"](https://journals.ametsoc.org/view/journals/atsc/33/12/1520-0469_1976_033_2452_tdeafr_2_0_co_2.xml)

## Code Structure

The implementation uses a modular approach with separate inline files for each function:

- `delta_eddington.hpp`: Main header file with class definition
- `delta_eddington_constructor.inl`: Constructor implementation
- `delta_eddington_initialize.inl`: Initialization function
- `delta_eddington_calculate.inl`: Main calculation function
- `delta_eddington_apply_delta_scaling.inl`: Delta-Eddington scaling implementation
- `delta_eddington_calculate_column.inl`: Column calculation algorithm

## Mathematical Implementation

### Delta-Eddington Scaling (`delta_eddington_apply_delta_scaling.inl`)

The delta-scaling transformation adjusts optical properties to account for forward scattering:

```cpp
double f = g * g;
double scaled_tau = tau * (1.0 - omega * f);
double scaled_omega = omega * (1.0 - f) / (1.0 - omega * f);
double scaled_g = g / (1.0 + g);
```

The above implements the delta-Eddington scaling from Joseph et al. (1976).

### Gamma Coefficients (`delta_eddington_calculate_column.inl`)

The Eddington coefficients are calculated as:

```cpp
gamma1[i] = (7.0 - omega[i] * (4.0 + 3.0 * g[i])) / 4.0;  // γ₁ = (7 - ω₀(4+3g))/4
gamma2[i] = -(1.0 - omega[i] * (4.0 - 3.0 * g[i])) / 4.0; // γ₂ = -(1 - ω₀(4-3g))/4
gamma3[i] = (2.0 - 3.0 * g[i] * mu0) / 4.0;               // γ₃ = (2-3gμ₀)/4
double gamma4 = 1.0 - gamma3[i];                          // γ₄ = 1 - γ₃
```

These correspond to the γ coefficients in Step 1 of the algorithm.

### Lambda and Exponential Terms (`delta_eddington_calculate_column.inl`)

The λ parameter and exponential terms are calculated as:

```cpp
lambda[i] = std::sqrt(gamma1[i] * gamma1[i] - gamma2[i] * gamma2[i]); // λ = √(γ₁² - γ₂²)

e1[i] = 1.0 + gamma1[i] * optical_depths[i] / lambda[i];             // e₁ = 1 + Γe^(-λτ)
e2[i] = 1.0 - gamma1[i] * optical_depths[i] / lambda[i];             // e₂ = 1 - Γe^(-λτ)
e3[i] = std::exp(-lambda[i] * optical_depths[i]);                    // e₃ = e^(-λτ)
e4[i] = 1.0 / e3[i];                                                 // e₄ = 1/e₃
```

### Source Terms C+ and C- (`delta_eddington_calculate_column.inl`)

The source terms represent the contribution from the direct solar beam:

```cpp
double denom = lambda[i] * lambda[i] - 1.0 / (mu0 * mu0);
double exp_tau = std::exp(-optical_depths[i] / mu0);

c_plus[i] = omega[i] * exp_tau *
           ((gamma1[i] - 1.0/mu0) * gamma3[i] + gamma2[i] * gamma4) / denom;

c_minus[i] = omega[i] * exp_tau *
            ((gamma1[i] + 1.0/mu0) * gamma4 + gamma2[i] * gamma3[i]) / denom;
```

This implements the C+ and C- functions from Step 2.

### Tridiagonal System Setup (`delta_eddington_calculate_column.inl`)

The tridiagonal system coefficients (a, b, c, r) correspond to A, B, D, E in Steps 3-4:

```cpp
// Top boundary condition
a[0] = 0.0;                 // A₁ = 0
b[0] = 1.0;                 // B₁ = e₁₁
c[0] = -e3[0];              // D₁ = -e₂₁
r[0] = 0.0;                 // E₁ = F₀⁻(0) - C₁⁻(0)

// Interface conditions
for (size_t i = 0; i < n_layers - 1; ++i) {
    // Continuity of upward flux (n=odd)
    size_t j = 2 * i + 1;
    a[j] = e3[i] * (1.0 + lambda[i] * gamma2[i]);
    b[j] = e4[i] * (1.0 - lambda[i] * gamma2[i]);
    c[j] = -(1.0 + lambda[i+1] * gamma2[i+1]);
    r[j] = c_plus[i+1] - e3[i] * c_plus[i] - e4[i] * c_minus[i];
    
    // Continuity of downward flux (n=even)
    j = 2 * i + 2;
    a[j] = e3[i] * (1.0 - lambda[i] * gamma2[i]);
    b[j] = e4[i] * (1.0 + lambda[i] * gamma2[i]);
    c[j] = -(1.0 - lambda[i+1] * gamma2[i+1]);
    r[j] = c_minus[i+1] - e3[i] * c_minus[i] - e4[i] * c_plus[i];
}

// Bottom boundary condition
size_t j = 2 * n_layers - 1;
a[j] = 1.0 - lambda[n_layers-1] * gamma2[n_layers-1];
b[j] = albedo * (1.0 + lambda[n_layers-1] * gamma2[n_layers-1]);
c[j] = 0.0;
r[j] = albedo * direct_flux[n_layers] - c_minus[n_layers-1] - albedo * c_plus[n_layers-1];
```

### Tridiagonal System Solution (`delta_eddington_calculate_column.inl`)

We use tuv-x's existing tridiagonal solver:

```cpp
tri_solver_.solve(a, b, c, r, u);
```

This solves for the Y coefficients in Step 5.

### Flux Calculation (`delta_eddington_calculate_column.inl`)

Finally, calculate fluxes using the solution coefficients (Step 6):

```cpp
// Calculate diffuse fluxes using equations (31) and (32)
diffuse_up_flux[i] = y1 * std::exp(lambda[layer] * (tau - optical_depths[layer])) +
                     y2 * exp_lambda_tau +
                     c_plus[layer] * exp_beam_tau;

diffuse_down_flux[i] = y1 * gamma2[layer] / (gamma1[layer] - lambda[layer]) * 
                       std::exp(lambda[layer] * (tau - optical_depths[layer])) +
                       y2 * gamma2[layer] / (gamma1[layer] + lambda[layer]) *
                       exp_lambda_tau +
                       c_minus[layer] * exp_beam_tau;

// Calculate actinic flux (spherical irradiance)
spherical_flux[i] = direct_flux[i] / mu0 + 
                   2.0 * (diffuse_up_flux[i] + diffuse_down_flux[i]);
```

## Configuration Options

- `radiative_transfer.delta_eddington.use_delta_scaling`: Whether to use delta-scaling (default: true)
- `radiative_transfer.delta_eddington.solar_irradiance`: Solar irradiance at top of atmosphere (default: 1.0)

## Optimization Features

- Optimized memory layout: [column][wavelength][vertical layer] for better vectorization
- Uses Array3D for efficient multi-dimensional array access
- Leverages tuv-x's existing TridiagonalSolver for the linear system
- All inline implementation for compiler optimization opportunities
