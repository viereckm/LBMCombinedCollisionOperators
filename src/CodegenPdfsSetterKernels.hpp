#pragma once
#include "core/DataTypes.h"
#include "core/logging/Logging.h"
#include "core/Macros.h"

#include "field/GhostLayerField.h"
#include "field/SwapableCompare.h"
#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include <set>

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

namespace walberla {
namespace lbm {
namespace codegenKernels {

/**
 * @brief funcion setting BGK-LBM (SRT) PDF values to equilibrium corresponding to a velocity field
 * 
 * @param ctr_0 x coordinate of currently processed cell
 * @param ctr_1 y coordinate of currently processed cell
 * @param _data_pdfs PDF field to be initialized
 * @param _data_velocity velocity field PDF values are set according to
 * @param _stride_pdfs_0 x-stride of PDF field
 * @param _stride_pdfs_1 y-stride of PDF field
 * @param _stride_pdfs_2 population stride of PDF field
 * @param _stride_velocity_0 x-stride of velocity field
 * @param _stride_velocity_1 y-stride of velocity field
 * @param _stride_velocity_2 component stride of velocity field
 * @param rho_0 density
 */
void D2Q9SRTPdfsSetter(
    int64_t const ctr_0,
    int64_t const ctr_1, 
    double * RESTRICT _data_pdfs,
    double * RESTRICT const _data_velocity,
    int64_t const _stride_pdfs_0,
    int64_t const _stride_pdfs_1,
    int64_t const _stride_pdfs_2,
    int64_t const _stride_velocity_0, 
    int64_t const _stride_velocity_1, 
    int64_t const _stride_velocity_2, 
    double rho_0
);

/**
 * @brief funcion setting TRT-LBM PDF values to equilibrium corresponding to a velocity field (same as for SRT but sti implemented fr consistency reasons)
 * 
 * @param ctr_0 x coordinate of currently processed cell
 * @param ctr_1 y coordinate of currently processed cell
 * @param _data_pdfs PDF field to be initialized
 * @param _data_velocity velocity field PDF values are set according to
 * @param _stride_pdfs_0 x-stride of PDF field
 * @param _stride_pdfs_1 y-stride of PDF field
 * @param _stride_pdfs_2 population stride of PDF field
 * @param _stride_velocity_0 x-stride of velocity field
 * @param _stride_velocity_1 y-stride of velocity field
 * @param _stride_velocity_2 component stride of velocity field
 * @param rho_0 density
 */
void D2Q9TRTPdfsSetter(
    int64_t const ctr_0,
    int64_t const ctr_1, 
    double * RESTRICT _data_pdfs,
    double * RESTRICT const _data_velocity,
    int64_t const _stride_pdfs_0,
    int64_t const _stride_pdfs_1,
    int64_t const _stride_pdfs_2,
    int64_t const _stride_velocity_0, 
    int64_t const _stride_velocity_1, 
    int64_t const _stride_velocity_2, 
    double rho_0
);

/**
 * @brief funcion setting Cumulant-LBM PDF values to equilibrium corresponding to a velocity field
 * 
 * @param ctr_0 x coordinate of currently processed cell
 * @param ctr_1 y coordinate of currently processed cell
 * @param _data_pdfs PDF field to be initialized
 * @param _data_velocity velocity field PDF values are set according to
 * @param _stride_pdfs_0 x-stride of PDF field
 * @param _stride_pdfs_1 y-stride of PDF field
 * @param _stride_pdfs_2 population stride of PDF field
 * @param _stride_velocity_0 x-stride of velocity field
 * @param _stride_velocity_1 y-stride of velocity field
 * @param _stride_velocity_2 component stride of velocity field
 * @param rho_0 density
 */
void D2Q9CumulantPdfsSetter(
    int64_t const ctr_0,
    int64_t const ctr_1, 
    double * RESTRICT _data_pdfs,
    double * RESTRICT const _data_velocity,
    int64_t const _stride_pdfs_0,
    int64_t const _stride_pdfs_1,
    int64_t const _stride_pdfs_2,
    int64_t const _stride_velocity_0, 
    int64_t const _stride_velocity_1, 
    int64_t const _stride_velocity_2, 
    double rho_0
);

} // !namespace codegenKernels
} // !namespace lbm
} // !namespace walberla