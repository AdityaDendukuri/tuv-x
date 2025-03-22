// Copyright (C) 2020-2023 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "tuvx/util/logging.hpp"

namespace tuvx {
namespace radiative_transfer {

inline void DeltaEddington::initialize(const util::ConfigMap& config,
                              const util::GridWarehouse& grid_warehouse,
                              const util::ProfileWarehouse& profile_warehouse)
{
    // Read configuration parameters
    use_delta_scaling_ = config.get_or_default("radiative_transfer.delta_eddington.use_delta_scaling", true);
    solar_irradiance_ = config.get_or_default("radiative_transfer.delta_eddington.solar_irradiance", 1.0);
    
    // Get grid sizes
    auto& height_grid = grid_warehouse.get_grid("height");
    auto& wavelength_grid = grid_warehouse.get_grid("wavelength");
    auto& zenith_grid = grid_warehouse.get_grid("zenith");
    
    n_heights_ = height_grid.size();
    n_wavelengths_ = wavelength_grid.size();
    n_columns_ = zenith_grid.size();
    n_layers_ = n_heights_ - 1;
    
    tuv_logger.info() << "Initialized Delta-Eddington radiative transfer solver with:";
    tuv_logger.info() << "  - " << n_columns_ << " columns (zenith angles)";
    tuv_logger.info() << "  - " << n_wavelengths_ << " wavelength points";
    tuv_logger.info() << "  - " << n_heights_ << " height levels";
    tuv_logger.info() << "  - Delta scaling: " << (use_delta_scaling_ ? "enabled" : "disabled");
}

} // namespace radiative_transfer
} // namespace tuvx
