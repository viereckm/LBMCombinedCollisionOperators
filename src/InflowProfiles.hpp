#pragma once

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "geometry/all.h"
#include "lbm/all.h"
#include "timeloop/all.h"
#include "stencil/D2Q9.h"

#include <cmath>

//#define _USE_MATH_DEFINES
#include <math.h>

namespace walberla {
namespace lbm {

/**
 * @brief Sets the velocity at the inflow boundary to a constant value uMax
 * 
 * @param uMax prescribed velocity
 * @return Vector3< real_t > 
 */
Vector3< real_t > constantInflow ( Vector3< real_t > uMax, real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t )
{
    return uMax;
}

/**
 * @brief function setting the inflow velocity according to the Schäfer Turek benchmark in 2D with constant inflow
 * 
 * @param uMax Vector of maximal inflow velocity in each direction (only x coordinate relevant here)
 * @param x x coordinate (not used here since inflow depends on y only)
 * @param y y cordinate
 * @param z z coordinate (not used here)
 * @param t timestep (not used here)
 * @param Nx number of cells in x direction (not used here)
 * @param Ny number of cells in y direction
 * @param Nz number of cells in z direction (not used here)
 * @param Nt number of timesteps (not used here)
 * @return Vector3 < real_t > inflow velocity in ( x, y, z ) at time t
 */
Vector3< real_t > constantParabolicProfile2D( Vector3 < real_t > uMax, real_t /* x */, real_t y, real_t  /* z */, real_t /* t */, real_t /* Nx */, real_t Ny, real_t /* Nz */, real_t /* Nt */ )
{
    real_t scaling = 4 * uMax[ 0 ] / ( Ny * Ny );
    real_t u = scaling * y * ( Ny - y );

    return Vector3< real_t >( u, real_t( 0.0 ), real_t( 0.0 ) );
}

/**
 * @brief function setting the inflow velocity according to the Schäfer Turek benchmark in 2D with varying inflow
 * 
 * @param uMax Vector of maximal inflow velocity in each direction (only x coordinate relevant here)
 * @param x x coordinate (not used here since inflow depends on y only)
 * @param y y cordinate
 * @param z z coordinate (not used here)
 * @param t timestep
 * @param Nx number of cells in x direction (not used here)
 * @param Ny number of cells in y direction
 * @param Nz number of cells in z direction (not used here)
 * @param Nt number of timesteps
 * @return Vector3 < real_t > inflow velocity in ( x, y, z ) at time t
 */
Vector3< real_t > varyingParabolicProfile2D( Vector3 < real_t > uMax, real_t /* x */, real_t y, real_t /* z */, real_t t, real_t /* Nx */, real_t Ny, real_t /* Nz */, real_t Nt )
{
    real_t scaling = 4 * uMax[ 0 ] / ( Ny * Ny );
    real_t oscillation = std::sin( M_PI * t / Nt );
    real_t u = oscillation * scaling * y * (Ny - y );

    return Vector3< real_t >( u, real_t ( 0.0 ), real_t( 0.0 ) );
}


} // namespace lbm
} // namespace walberla
