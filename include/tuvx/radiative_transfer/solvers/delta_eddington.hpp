// Copyright (C) 2020-2023 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Delta-Eddington solver for radiative transfer
// Based on Toon et al. (1989): "Rapid Calculation of Radiative Heating Rates and
// Photodissociation Rates in Inhomogeneous Multiple Scattering Atmospheres"

#pragma once

#include <vector>
#include <memory>
#include <cmath>
#include <tuple>

#include "tuvx/radiative_transfer/radiation_field.hpp"
#include "tuvx/linear_algebra/tridiagonal_solver.hpp"
#include "tuvx/array3d.hpp"

namespace tuvx {
namespace radiative_transfer {

///
/// @brief Implementation of the Delta-Eddington radiative transfer solver
///
/// This implements the Delta-Eddington approximation for solving the radiative
/// transfer equation in multiple scattering atmospheres as described in
/// Toon et al. (1989) JGR.
///
class DeltaEddington : public RadiationField
{
  public:
    /// @brief Constructor
    DeltaEddington();
    
    /// @brief Destructor
    ~DeltaEddington() override = default;

    /// @brief Initialize the Delta-Eddington solver
    ///
    /// @param config Configuration parameters
    /// @param grid_warehouse Grid warehouse containing the height, wavelength, and zenith grids
    /// @param profile_warehouse Profile warehouse containing the optical properties
    void initialize(const util::ConfigMap& config,
                  const util::GridWarehouse& grid_warehouse,
                  const util::ProfileWarehouse& profile_warehouse) override;

    /// @brief Calculate the radiation field at each point in the grid
    ///
    /// @param grid_warehouse Grid warehouse containing the height, wavelength, and zenith grids
    /// @param profile_warehouse Profile warehouse containing the optical properties
    /// @param photolysis_warehouse Photolysis warehouse for storing results
    void calculate(const util::GridWarehouse& grid_warehouse,
                 const util::ProfileWarehouse& profile_warehouse,
                 util::PhotolysisRateWarehouse& photolysis_warehouse) override;

  private:
    /// @brief Apply Delta-Eddington scaling to optical properties
    ///
    /// @param tau Original optical depth
    /// @param omega Original single scattering albedo
    /// @param g Original asymmetry factor
    /// @return Tuple of scaled optical depth, single scattering albedo, and asymmetry factor
    std::tuple<double, double, double> applyDeltaScaling(double tau, double omega, double g) const;

    /// @brief Calculate the radiation field for a single column
    ///
    /// @param mu0 Cosine of solar zenith angle
    /// @param albedo Surface albedo
    /// @param optical_depths Vector of optical depths for each layer
    /// @param omega Vector of single scattering albedos for each layer
    /// @param g Vector of asymmetry factors for each layer
    /// @param nLayers Number of layers
    /// @param direct_flux Output direct flux at each level
    /// @param diffuse_down_flux Output diffuse downward flux at each level
    /// @param diffuse_up_flux Output diffuse upward flux at each level
    /// @param spherical_flux Output spherical (actinic) flux at each level
    void calculateColumn(
        double mu0,
        double albedo,
        const std::vector<double>& optical_depths,
        const std::vector<double>& omega,
        const std::vector<double>& g,
        size_t nLayers,
        std::vector<double>& direct_flux,
        std::vector<double>& diffuse_down_flux,
        std::vector<double>& diffuse_up_flux,
        std::vector<double>& spherical_flux) const;

    // Configuration parameters
    bool use_delta_scaling_;            ///< Whether to use delta scaling
    double solar_irradiance_;           ///< Solar irradiance at top of atmosphere [W/m^2]
    
    // Array dimensioning parameters
    size_t n_columns_;                  ///< Number of columns (zenith angles)
    size_t n_wavelengths_;              ///< Number of wavelength points
    size_t n_heights_;                  ///< Number of height levels (vertical layers + 1)
    size_t n_layers_;                   ///< Number of layers (n_heights_ - 1)
    
    // Using the existing tridiagonal solver from tuv-x
    linear_algebra::TridiagonalSolver tri_solver_;
};

} // namespace radiative_transfer
} // namespace tuvx

// implementations
#include "delta_eddington_constructor.inl"
#include "delta_eddington_initialize.inl"
#include "delta_eddington_calculate.inl"
#include "delta_eddington_apply_delta_scaling.inl"
#include "delta_eddington_calculate_column.inl"
