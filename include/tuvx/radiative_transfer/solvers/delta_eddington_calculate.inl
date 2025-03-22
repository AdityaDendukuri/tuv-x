// Copyright (C) 2020-2023 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "tuvx/array3d.hpp"

namespace tuvx {
namespace radiative_transfer {

namespace {
  constexpr double PI = 3.14159265358979323846;
}

inline void DeltaEddington::calculate(const util::GridWarehouse& grid_warehouse,
                             const util::ProfileWarehouse& profile_warehouse,
                             util::PhotolysisRateWarehouse& photolysis_warehouse)
{
    // Get required grids
    auto& height_grid = grid_warehouse.get_grid("height");
    auto& wavelength_grid = grid_warehouse.get_grid("wavelength");
    auto& zenith_grid = grid_warehouse.get_grid("zenith");
    
    // Get required profiles
    auto& optical_depth = profile_warehouse.get_profile("optical_depth");
    auto& single_scattering_albedo = profile_warehouse.get_profile("single_scattering_albedo");
    auto& asymmetry_factor = profile_warehouse.get_profile("asymmetry_factor");
    auto& surface_albedo = profile_warehouse.get_profile("surface_albedo");
    
    // Initialize output arrays in the photolysis warehouse using Array3D
    // Use column (zenith), wavelength, height order for best vectorization
    auto& direct_flux = photolysis_warehouse.get_rate("direct_flux");
    auto& diffuse_down_flux = photolysis_warehouse.get_rate("diffuse_down_flux");
    auto& diffuse_up_flux = photolysis_warehouse.get_rate("diffuse_up_flux");
    auto& spherical_flux = photolysis_warehouse.get_rate("actinic_flux");
    
    // Ensure arrays are sized correctly
    if (direct_flux.size() != n_columns_ * n_wavelengths_ * n_heights_) {
        direct_flux.resize(n_columns_ * n_wavelengths_ * n_heights_);
        diffuse_down_flux.resize(n_columns_ * n_wavelengths_ * n_heights_);
        diffuse_up_flux.resize(n_columns_ * n_wavelengths_ * n_heights_);
        spherical_flux.resize(n_columns_ * n_wavelengths_ * n_heights_);
    }
    
    // Temporary vectors for a single column calculation
    std::vector<double> direct_flux_col(n_heights_);
    std::vector<double> diffuse_down_flux_col(n_heights_);
    std::vector<double> diffuse_up_flux_col(n_heights_);
    std::vector<double> spherical_flux_col(n_heights_);
    
    // Arrays for layer properties
    std::vector<double> layer_tau(n_layers_);
    std::vector<double> layer_omega(n_layers_);
    std::vector<double> layer_g(n_layers_);
    
    // Create Array3D views for each output array
    Array3D direct_flux_3d(direct_flux.data(), n_columns_, n_wavelengths_, n_heights_);
    Array3D diffuse_down_flux_3d(diffuse_down_flux.data(), n_columns_, n_wavelengths_, n_heights_);
    Array3D diffuse_up_flux_3d(diffuse_up_flux.data(), n_columns_, n_wavelengths_, n_heights_);
    Array3D spherical_flux_3d(spherical_flux.data(), n_columns_, n_wavelengths_, n_heights_);
    
    // Loop over all columns (zenith angles) first, then wavelengths, then vertical layers
    // This follows the memory layout for optimal vectorization
    for (size_t ic = 0; ic < n_columns_; ++ic) {
        // Calculate cosine of solar zenith angle (degrees to radians)
        double mu0 = std::cos(zenith_grid(ic) * PI / 180.0);
        
        for (size_t iw = 0; iw < n_wavelengths_; ++iw) {
            // Get surface albedo for this wavelength
            double albedo = surface_albedo(iw);
            
            // Extract layer properties for this wavelength
            for (size_t i = 0; i < n_layers_; ++i) {
                // Get optical properties from profiles
                size_t od_idx = optical_depth.map({ iw, i });
                size_t ssa_idx = single_scattering_albedo.map({ iw, i });
                size_t g_idx = asymmetry_factor.map({ iw, i });
                
                // Get the optical properties for the current layer
                double tau = optical_depth(od_idx);
                double omega = single_scattering_albedo(ssa_idx);
                double g = asymmetry_factor(g_idx);
                
                // Apply delta-Eddington scaling if requested
                if (use_delta_scaling_) {
                    std::tie(tau, omega, g) = applyDeltaScaling(tau, omega, g);
                }
                
                layer_tau[i] = tau;
                layer_omega[i] = omega;
                layer_g[i] = g;
            }
            
            // Calculate the radiation field for this column
            calculateColumn(
                mu0, 
                albedo, 
                layer_tau, 
                layer_omega, 
                layer_g, 
                n_layers_, 
                direct_flux_col, 
                diffuse_down_flux_col, 
                diffuse_up_flux_col, 
                spherical_flux_col
            );
            
            // Store the results using Array3D indexing
            // Using column, wavelength, height ordering for best vectorization
            for (size_t i = 0; i < n_heights_; ++i) {
                // Set the values using Array3D indexing
                direct_flux_3d(ic, iw, i) = direct_flux_col[i] * solar_irradiance_;
                diffuse_down_flux_3d(ic, iw, i) = diffuse_down_flux_col[i] * solar_irradiance_;
                diffuse_up_flux_3d(ic, iw, i) = diffuse_up_flux_col[i] * solar_irradiance_;
                spherical_flux_3d(ic, iw, i) = spherical_flux_col[i] * solar_irradiance_;
            }
        }
    }
}

} // namespace radiative_transfer
} // namespace tuvx
