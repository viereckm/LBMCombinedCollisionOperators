//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \\file D2Q9CodegenNoSlip.cpp
//! \\author pystencils
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "D2Q9CodegenNoSlip.hpp"



#define FUNC_PREFIX

using namespace std;

namespace walberla {
namespace lbm {

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wconversion"
#endif

#ifdef __CUDACC__
#pragma push
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diag_suppress 177
#else
#pragma diag_suppress 177
#endif
#endif
//NOLINTBEGIN(readability-non-const-parameter*)
namespace internal_d2q9codegennoslip_boundary_D2Q9CodegenNoSlip {
static FUNC_PREFIX void d2q9codegennoslip_boundary_D2Q9CodegenNoSlip(uint8_t * RESTRICT const _data_indexVector, double * RESTRICT  _data_pdfs, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int32_t indexVectorSize)
{
   
   const int32_t f_in_inv_dir_idx [] = { 0,2,1,4,3,8,7,6,5 }; 
   const int32_t f_in_inv_offsets_x [] = { 0,0,0,-1,1,-1,1,-1,1 }; 
   const int32_t f_in_inv_offsets_y [] = { 0,1,-1,0,0,1,1,-1,-1 }; 
   
   for (int64_t ctr_0 = 0; ctr_0 < indexVectorSize; ctr_0 += 1)
   {
      const int32_t x = *((int32_t *  )(& _data_indexVector[12*ctr_0]));
      const int32_t y = *((int32_t *  )(& _data_indexVector[12*ctr_0 + 4]));
      const int32_t dir = *((int32_t *  )(& _data_indexVector[12*ctr_0 + 8]));
      _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_0*f_in_inv_offsets_x[dir] + _stride_pdfs_1*y + _stride_pdfs_1*f_in_inv_offsets_y[dir] + _stride_pdfs_2*f_in_inv_dir_idx[dir]] = _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*dir];
   }
}
}

//NOLINTEND(readability-non-const-parameter*)
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#ifdef __CUDACC__
#pragma pop
#endif


void D2Q9CodegenNoSlip::run_impl(IBlock * block, IndexVectors::Type type)
{
   auto * indexVectors = block->getData<IndexVectors>(indexVectorID);
   int32_t indexVectorSize = int32_c( indexVectors->indexVector(type).size() );
   if( indexVectorSize == 0)
      return;

   
   auto pointer = indexVectors->pointerCpu(type);
   

   uint8_t * _data_indexVector = reinterpret_cast<uint8_t*>(pointer);

   auto pdfs = block->getData< field::GhostLayerField<double, 9> >(pdfsID);

   
   
   WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_d2q9codegennoslip_boundary_D2Q9CodegenNoSlip::d2q9codegennoslip_boundary_D2Q9CodegenNoSlip(_data_indexVector, _data_pdfs, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, indexVectorSize);
}

void D2Q9CodegenNoSlip::run(IBlock * block)
{
   run_impl(block, IndexVectors::ALL);
}

void D2Q9CodegenNoSlip::inner(IBlock * block)
{
   run_impl(block, IndexVectors::INNER);
}

void D2Q9CodegenNoSlip::outer(IBlock * block)
{
   run_impl(block, IndexVectors::OUTER);
}

} // namespace lbm
} // namespace walberla

