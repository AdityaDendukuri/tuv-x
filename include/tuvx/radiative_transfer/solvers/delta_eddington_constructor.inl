// Copyright (C) 2020-2023 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace tuvx {
namespace radiative_transfer {

inline DeltaEddington::DeltaEddington()
    : use_delta_scaling_(true)
    , solar_irradiance_(1.0)
    , n_columns_(0)
    , n_wavelengths_(0)
    , n_heights_(0)
    , n_layers_(0)
    , tri_solver_()
{}

} // namespace radiative_transfer
} // namespace tuvx
