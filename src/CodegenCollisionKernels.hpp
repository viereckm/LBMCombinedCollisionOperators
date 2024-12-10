#pragma once

#include "core/DataTypes.h"
#include "core/Macros.h"

#include <cmath>

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
 * @brief function performing BGK-LBM (SRT) collision in one cell of a D2Q9 lattice
 * 
 * @param ctr_0 x-coordinate
 * @param ctr_1 y-coodinate
 * @param _data_density density field
 * @param _data_pdfs pre-collision PDF field
 * @param _data_pdfs_tmp post-collision PDF field
 * @param _data_velocity macroscopic velocity field
 * @param _stride_density_0 stride in x-dimension in density field
 * @param _stride_density_1 stride in y-dimension in density field
 * @param _stride_pdfs_0 stride in x-dimension in pre-collision field
 * @param _stride_pdfs_1 stride in y-dimension in pre-collision field
 * @param _stride_pdfs_2 stride in f-dimension (populations) in pre-collision field
 * @param _stride_pdfs_tmp_0 stride in x-dimension in post-collision field
 * @param _stride_pdfs_tmp_1 stride in y-dimension in post-collision field
 * @param _stride_pdfs_tmp_2 stride in f-dimension (populations) in post-collision field
 * @param _stride_velocity_0 stride in x-dimension in velocity field
 * @param _stride_velocity_1 stride in y-dimension in velocity field
 * @param _stride_velocity_2 stride in component dimension in velocity field
 * @param omega relaxation parameter
 */
void D2Q9SRTKernel(
    int64_t ctr_0,
    int64_t ctr_1,
    double * RESTRICT const _data_density, 
    double * RESTRICT const _data_pdfs, 
    double * RESTRICT const _data_pdfs_tmp, 
    double * RESTRICT const _data_velocity,
    int64_t const _stride_density_0,
    int64_t const _stride_density_1,
    int64_t const _stride_pdfs_0, 
    int64_t const _stride_pdfs_1, 
    int64_t const _stride_pdfs_2, 
    int64_t const _stride_pdfs_tmp_0, 
    int64_t const _stride_pdfs_tmp_1, 
    int64_t const _stride_pdfs_tmp_2,
    int64_t const _stride_velocity_0, 
    int64_t const _stride_velocity_1, 
    int64_t const _stride_velocity_2, 
    double omega 
);

/**
 * @brief function performing TRT collision with magic number 3/16 in one cell of a D2Q9 lattice
 * 
 * @param ctr_0 x-coordinate
 * @param ctr_1 y-coodinate
 * @param _data_density density field
 * @param _data_pdfs pre-collision PDF field
 * @param _data_pdfs_tmp post-collision PDF field
 * @param _data_velocity macroscopic velocity field
 * @param _stride_density_0 stride in x-dimension in density field
 * @param _stride_density_1 stride in y-dimension in density field
 * @param _stride_pdfs_0 stride in x-dimension in pre-collision field
 * @param _stride_pdfs_1 stride in y-dimension in pre-collision field
 * @param _stride_pdfs_2 stride in f-dimension (populations) in pre-collision field
 * @param _stride_pdfs_tmp_0 stride in x-dimension in post-collision field
 * @param _stride_pdfs_tmp_1 stride in y-dimension in post-collision field
 * @param _stride_pdfs_tmp_2 stride in f-dimension (populations) in post-collision field
 * @param _stride_velocity_0 stride in x-dimension in velocity field
 * @param _stride_velocity_1 stride in y-dimension in velocity field
 * @param _stride_velocity_2 stride in component dimension in velocity field
 * @param omega relaxation parameter
 */
void D2Q9TRTKernel(
    int64_t ctr_0,
    int64_t ctr_1,
    double * RESTRICT const _data_density, 
    double * RESTRICT const _data_pdfs, 
    double * RESTRICT const _data_pdfs_tmp, 
    double * RESTRICT const _data_velocity,
    int64_t const _stride_density_0,
    int64_t const _stride_density_1,
    int64_t const _stride_pdfs_0, 
    int64_t const _stride_pdfs_1, 
    int64_t const _stride_pdfs_2, 
    int64_t const _stride_pdfs_tmp_0, 
    int64_t const _stride_pdfs_tmp_1, 
    int64_t const _stride_pdfs_tmp_2,
    int64_t const _stride_velocity_0, 
    int64_t const _stride_velocity_1, 
    int64_t const _stride_velocity_2, 
    double omega 
);

/**
 * @brief function performing cumulant collision in one cell of a D2Q9 lattice
 * 
 * @param ctr_0 x-coordinate
 * @param ctr_1 y-coodinate
 * @param _data_density density field
 * @param _data_pdfs pre-collision PDF field
 * @param _data_pdfs_tmp post-collision PDF field
 * @param _data_velocity macroscopic velocity field
 * @param _stride_density_0 stride in x-dimension in density field
 * @param _stride_density_1 stride in y-dimension in density field
 * @param _stride_pdfs_0 stride in x-dimension in pre-collision field
 * @param _stride_pdfs_1 stride in y-dimension in pre-collision field
 * @param _stride_pdfs_2 stride in f-dimension (populations) in pre-collision field
 * @param _stride_pdfs_tmp_0 stride in x-dimension in post-collision field
 * @param _stride_pdfs_tmp_1 stride in y-dimension in post-collision field
 * @param _stride_pdfs_tmp_2 stride in f-dimension (populations) in post-collision field
 * @param _stride_velocity_0 stride in x-dimension in velocity field
 * @param _stride_velocity_1 stride in y-dimension in velocity field
 * @param _stride_velocity_2 stride in component dimension in velocity field
 * @param omega relaxation parameter
 */
void D2Q9CumulantKernel(
    int64_t ctr_0,
    int64_t ctr_1,
    double * RESTRICT const _data_density, 
    double * RESTRICT const _data_pdfs, 
    double * RESTRICT const _data_pdfs_tmp, 
    double * RESTRICT const _data_velocity,
    int64_t const _stride_density_0,
    int64_t const _stride_density_1,
    int64_t const _stride_pdfs_0, 
    int64_t const _stride_pdfs_1, 
    int64_t const _stride_pdfs_2, 
    int64_t const _stride_pdfs_tmp_0, 
    int64_t const _stride_pdfs_tmp_1, 
    int64_t const _stride_pdfs_tmp_2,
    int64_t const _stride_velocity_0, 
    int64_t const _stride_velocity_1, 
    int64_t const _stride_velocity_2, 
    double omega 
);

} // !namespace codegenKernels
} // !namespace lbm
} // !namespace walberla