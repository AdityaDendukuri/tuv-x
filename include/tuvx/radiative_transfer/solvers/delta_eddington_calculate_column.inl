// Copyright (C) 2020-2023 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace tuvx {
namespace radiative_transfer {

namespace {
  constexpr double TINY = 1.0e-10;
}

inline void DeltaEddington::calculateColumn(
    double mu0,
    double albedo,
    const std::vector<double>& optical_depths,
    const std::vector<double>& omega,
    const std::vector<double>& g,
    size_t n_layers,
    std::vector<double>& direct_flux,
    std::vector<double>& diffuse_down_flux,
    std::vector<double>& diffuse_up_flux,
    std::vector<double>& spherical_flux) const
{
    // Implementation of the Toon et al. (1989) algorithm
    
    // Initialize arrays for solver
    std::vector<double> gamma1(n_layers);
    std::vector<double> gamma2(n_layers);
    std::vector<double> gamma3(n_layers);
    std::vector<double> lambda(n_layers);
    std::vector<double> e1(n_layers);
    std::vector<double> e2(n_layers);
    std::vector<double> e3(n_layers);
    std::vector<double> e4(n_layers);
    
    // Arrays for source terms C+ and C-
    std::vector<double> c_plus(n_layers);
    std::vector<double> c_minus(n_layers);
    
    // Arrays for tridiagonal system
    std::vector<double> a(2 * n_layers);
    std::vector<double> b(2 * n_layers);
    std::vector<double> c(2 * n_layers);
    std::vector<double> r(2 * n_layers);
    std::vector<double> u(2 * n_layers);
    
    // Compute cumulative optical depth
    std::vector<double> cum_tau(n_layers + 1, 0.0);
    for (size_t i = 0; i < n_layers; ++i) {
        cum_tau[i + 1] = cum_tau[i] + optical_depths[i];
    }
    
    // Calculate direct solar beam (Beer-Lambert law)
    for (size_t i = 0; i <= n_layers; ++i) {
        direct_flux[i] = std::exp(-cum_tau[i] / mu0);
    }
    
    // Calculate the coefficients for each layer (Meador & Weaver, 1980)
    for (size_t i = 0; i < n_layers; ++i) {
        // Eddington coefficients from Table 1 in Toon et al. (1989)
        gamma1[i] = (7.0 - omega[i] * (4.0 + 3.0 * g[i])) / 4.0;
        gamma2[i] = -(1.0 - omega[i] * (4.0 - 3.0 * g[i])) / 4.0;
        gamma3[i] = (2.0 - 3.0 * g[i] * mu0) / 4.0;
        double gamma4 = 1.0 - gamma3[i];
        
        // Calculate lambda from equation (21)
        lambda[i] = std::sqrt(gamma1[i] * gamma1[i] - gamma2[i] * gamma2[i]);
        
        // Equation (44) - e1, e2, e3, e4 terms
        e1[i] = 1.0 + gamma1[i] * optical_depths[i] / lambda[i];
        e2[i] = 1.0 - gamma1[i] * optical_depths[i] / lambda[i];
        e3[i] = std::exp(-lambda[i] * optical_depths[i]);
        e4[i] = 1.0 / e3[i];
        
        // Calculate C+ and C- source terms from equations (23) and (24)
        double denom = lambda[i] * lambda[i] - 1.0 / (mu0 * mu0);
        
        // Avoid division by near-zero
        if (std::abs(denom) < TINY) {
            denom = TINY;
        }
        
        // Solar source terms (equation 23 and 24)
        double exp_tau = std::exp(-optical_depths[i] / mu0);
        
        c_plus[i] = omega[i] * exp_tau *
                   ((gamma1[i] - 1.0/mu0) * gamma3[i] + gamma2[i] * gamma4) / denom;
        
        c_minus[i] = omega[i] * exp_tau *
                    ((gamma1[i] + 1.0/mu0) * gamma4 + gamma2[i] * gamma3[i]) / denom;
    }
    
    // Set up the tridiagonal system following equations (41-43)
    
    // Top boundary condition - no diffuse flux from above
    a[0] = 0.0;
    b[0] = 1.0;
    c[0] = -e3[0];
    r[0] = 0.0;
    
    // Interface conditions between layers
    for (size_t i = 0; i < n_layers - 1; ++i) {
        // Continuity of upward flux
        size_t j = 2 * i + 1;
        a[j] = e3[i] * (1.0 + lambda[i] * gamma2[i]);
        b[j] = e4[i] * (1.0 - lambda[i] * gamma2[i]);
        c[j] = -(1.0 + lambda[i+1] * gamma2[i+1]);
        r[j] = c_plus[i+1] - e3[i] * c_plus[i] - e4[i] * c_minus[i];
        
        // Continuity of downward flux
        j = 2 * i + 2;
        a[j] = e3[i] * (1.0 - lambda[i] * gamma2[i]);
        b[j] = e4[i] * (1.0 + lambda[i] * gamma2[i]);
        c[j] = -(1.0 - lambda[i+1] * gamma2[i+1]);
        r[j] = c_minus[i+1] - e3[i] * c_minus[i] - e4[i] * c_plus[i];
    }
    
    // Bottom boundary condition with surface albedo
    size_t j = 2 * n_layers - 1;
    a[j] = 1.0 - lambda[n_layers-1] * gamma2[n_layers-1];
    b[j] = albedo * (1.0 + lambda[n_layers-1] * gamma2[n_layers-1]);
    c[j] = 0.0;
    r[j] = albedo * direct_flux[n_layers] - c_minus[n_layers-1] - albedo * c_plus[n_layers-1];
    
    // Solve the tridiagonal system using the existing solver
    tri_solver_.solve(a, b, c, r, u);
    
    // Calculate the diffuse fluxes using the solution from the tridiagonal system
    // This follows from equations (31) and (32)
    for (size_t i = 0; i <= n_layers; ++i) {
        // Determine which layer this level is in
        size_t layer = 0;
        if (i < n_layers) {
            layer = i;
        } else {
            layer = n_layers - 1;
        }
        
        // Calculate the diffuse upward and downward fluxes
        double y1 = u[2 * layer];       // Y1 in equation (29)
        double y2 = u[2 * layer + 1];   // Y2 in equation (30)
        
        // Depth within the layer (0 at top, optical_depths[layer] at bottom)
        double tau = i < n_layers ? 0.0 : optical_depths[layer];
        
        // Calculate the exponential terms for this depth
        double exp_lambda_tau = std::exp(-lambda[layer] * tau);
        double exp_beam_tau = std::exp(-tau / mu0);
        
        // Calculate diffuse fluxes using equations (31) and (32)
        diffuse_up_flux[i] = y1 * std::exp(lambda[layer] * (tau - optical_depths[layer])) +
                             y2 * exp_lambda_tau +
                             c_plus[layer] * exp_beam_tau;
        
        diffuse_down_flux[i] = y1 * gamma2[layer] / (gamma1[layer] - lambda[layer]) * 
                               std::exp(lambda[layer] * (tau - optical_depths[layer])) +
                               y2 * gamma2[layer] / (gamma1[layer] + lambda[layer]) *
                               exp_lambda_tau +
                               c_minus[layer] * exp_beam_tau;
        
        // Calculate actinic flux (4π times the mean intensity)
        // This is the sum of direct (converted from normal to spherical) and diffuse components
        spherical_flux[i] = direct_flux[i] / mu0 + 
                           2.0 * (diffuse_up_flux[i] + diffuse_down_flux[i]);
    }
}

} // namespace radiative_transfer
} // namespace tuvx
