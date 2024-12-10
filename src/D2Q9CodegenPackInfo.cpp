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
//! \\file D2Q9CodegenPackInfo.cpp
//! \\author pystencils
//======================================================================================================================

#include "stencil/Directions.h"
#include "core/cell/CellInterval.h"
#include "core/DataTypes.h"
#include "D2Q9CodegenPackInfo.hpp"

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wfloat-equal"
#   pragma GCC diagnostic ignored "-Wshadow"
#   pragma GCC diagnostic ignored "-Wconversion"
#   pragma GCC diagnostic ignored "-Wunused-variable"
#endif



namespace walberla {
namespace pystencils {

using walberla::cell::CellInterval;
using walberla::stencil::Direction;



namespace internal_pack_SW {
static FUNC_PREFIX void pack_SW(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
   {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
      {
         _data_buffer[_size_pdfs_0*ctr_1 + ctr_0] = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 7*_stride_pdfs_2];
      }
   }
}
}

namespace internal_pack_W {
static FUNC_PREFIX void pack_W(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
   {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
      {
         _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0] = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 3*_stride_pdfs_2];
         _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0 + 1] = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 5*_stride_pdfs_2];
         _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0 + 2] = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 7*_stride_pdfs_2];
      }
   }
}
}

namespace internal_pack_NW {
static FUNC_PREFIX void pack_NW(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
   {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
      {
         _data_buffer[_size_pdfs_0*ctr_1 + ctr_0] = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 5*_stride_pdfs_2];
      }
   }
}
}

namespace internal_pack_S {
static FUNC_PREFIX void pack_S(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
   {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
      {
         _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0] = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 2*_stride_pdfs_2];
         _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0 + 1] = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 7*_stride_pdfs_2];
         _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0 + 2] = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 8*_stride_pdfs_2];
      }
   }
}
}

namespace internal_pack_N {
static FUNC_PREFIX void pack_N(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
   {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
      {
         _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0] = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2];
         _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0 + 1] = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 5*_stride_pdfs_2];
         _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0 + 2] = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 6*_stride_pdfs_2];
      }
   }
}
}

namespace internal_pack_SE {
static FUNC_PREFIX void pack_SE(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
   {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
      {
         _data_buffer[_size_pdfs_0*ctr_1 + ctr_0] = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 8*_stride_pdfs_2];
      }
   }
}
}

namespace internal_pack_E {
static FUNC_PREFIX void pack_E(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
   {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
      {
         _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0] = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 4*_stride_pdfs_2];
         _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0 + 1] = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 6*_stride_pdfs_2];
         _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0 + 2] = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 8*_stride_pdfs_2];
      }
   }
}
}

namespace internal_pack_NE {
static FUNC_PREFIX void pack_NE(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
   {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
      {
         _data_buffer[_size_pdfs_0*ctr_1 + ctr_0] = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 6*_stride_pdfs_2];
      }
   }
}
}



namespace internal_unpack_SW {
static FUNC_PREFIX void unpack_SW(double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
   {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
      {
         _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 7*_stride_pdfs_2] = _data_buffer[_size_pdfs_0*ctr_1 + ctr_0];
      }
   }
}
}

namespace internal_unpack_W {
static FUNC_PREFIX void unpack_W(double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
   {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
      {
         _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 3*_stride_pdfs_2] = _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0];
         _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 5*_stride_pdfs_2] = _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0 + 1];
         _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 7*_stride_pdfs_2] = _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0 + 2];
      }
   }
}
}

namespace internal_unpack_NW {
static FUNC_PREFIX void unpack_NW(double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
   {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
      {
         _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 5*_stride_pdfs_2] = _data_buffer[_size_pdfs_0*ctr_1 + ctr_0];
      }
   }
}
}

namespace internal_unpack_S {
static FUNC_PREFIX void unpack_S(double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
   {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
      {
         _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 2*_stride_pdfs_2] = _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0];
         _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 7*_stride_pdfs_2] = _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0 + 1];
         _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 8*_stride_pdfs_2] = _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0 + 2];
      }
   }
}
}

namespace internal_unpack_N {
static FUNC_PREFIX void unpack_N(double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
   {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
      {
         _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2] = _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0];
         _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 5*_stride_pdfs_2] = _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0 + 1];
         _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 6*_stride_pdfs_2] = _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0 + 2];
      }
   }
}
}

namespace internal_unpack_SE {
static FUNC_PREFIX void unpack_SE(double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
   {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
      {
         _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 8*_stride_pdfs_2] = _data_buffer[_size_pdfs_0*ctr_1 + ctr_0];
      }
   }
}
}

namespace internal_unpack_E {
static FUNC_PREFIX void unpack_E(double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
   {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
      {
         _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 4*_stride_pdfs_2] = _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0];
         _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 6*_stride_pdfs_2] = _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0 + 1];
         _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 8*_stride_pdfs_2] = _data_buffer[3*_size_pdfs_0*ctr_1 + 3*ctr_0 + 2];
      }
   }
}
}

namespace internal_unpack_NE {
static FUNC_PREFIX void unpack_NE(double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
   {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
      {
         _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 6*_stride_pdfs_2] = _data_buffer[_size_pdfs_0*ctr_1 + ctr_0];
      }
   }
}
}



void D2Q9CodegenPackInfo::pack(Direction dir, unsigned char * byte_buffer, IBlock * block) const
{
    double * buffer = reinterpret_cast<double*>(byte_buffer);

    auto pdfs = block->getData< field::GhostLayerField<double, 9> >(pdfsID);

    CellInterval ci;
    pdfs->getSliceBeforeGhostLayer(dir, ci, 1, false);

    switch( dir )
    {
        case stencil::SW:
        {
            double * RESTRICT  _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
            double * RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
            const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
            const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
            internal_pack_SW::pack_SW(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        
        case stencil::W:
        {
            double * RESTRICT  _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
            double * RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
            const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
            const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
            internal_pack_W::pack_W(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        
        case stencil::NW:
        {
            double * RESTRICT  _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
            double * RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
            const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
            const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
            internal_pack_NW::pack_NW(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        
        case stencil::S:
        {
            double * RESTRICT  _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
            double * RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
            const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
            const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
            internal_pack_S::pack_S(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        
        case stencil::N:
        {
            double * RESTRICT  _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
            double * RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
            const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
            const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
            internal_pack_N::pack_N(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        
        case stencil::SE:
        {
            double * RESTRICT  _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
            double * RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
            const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
            const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
            internal_pack_SE::pack_SE(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        
        case stencil::E:
        {
            double * RESTRICT  _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
            double * RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
            const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
            const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
            internal_pack_E::pack_E(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        
        case stencil::NE:
        {
            double * RESTRICT  _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
            double * RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
            const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
            const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
            internal_pack_NE::pack_NE(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        

        default:
            WALBERLA_ASSERT(false);
    }
}


void D2Q9CodegenPackInfo::unpack(Direction dir, unsigned char * byte_buffer, IBlock * block) const
{
    double * buffer = reinterpret_cast<double*>(byte_buffer);

    auto pdfs = block->getData< field::GhostLayerField<double, 9> >(pdfsID);

    CellInterval ci;
    pdfs->getGhostRegion(dir, ci, 1, false);
    auto communciationDirection = stencil::inverseDir[dir];

    switch( communciationDirection )
    {
        case stencil::SW:
        {
            double * RESTRICT const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
            double * RESTRICT  _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
            const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
            const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
            internal_unpack_SW::unpack_SW(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        
        case stencil::W:
        {
            double * RESTRICT const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
            double * RESTRICT  _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
            const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
            const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
            internal_unpack_W::unpack_W(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        
        case stencil::NW:
        {
            double * RESTRICT const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
            double * RESTRICT  _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
            const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
            const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
            internal_unpack_NW::unpack_NW(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        
        case stencil::S:
        {
            double * RESTRICT const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
            double * RESTRICT  _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
            const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
            const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
            internal_unpack_S::unpack_S(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        
        case stencil::N:
        {
            double * RESTRICT const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
            double * RESTRICT  _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
            const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
            const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
            internal_unpack_N::unpack_N(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        
        case stencil::SE:
        {
            double * RESTRICT const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
            double * RESTRICT  _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
            const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
            const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
            internal_unpack_SE::unpack_SE(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        
        case stencil::E:
        {
            double * RESTRICT const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
            double * RESTRICT  _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
            const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
            const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
            internal_unpack_E::unpack_E(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        
        case stencil::NE:
        {
            double * RESTRICT const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
            double * RESTRICT  _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
            const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
            const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
            internal_unpack_NE::unpack_NE(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        

        default:
            WALBERLA_ASSERT(false);
    }
}


uint_t D2Q9CodegenPackInfo::size(stencil::Direction dir, const IBlock * block) const
{
    auto pdfs = block->getData< field::GhostLayerField<double, 9> >(pdfsID);

    CellInterval ci;
    pdfs->getGhostRegion(dir, ci, 1, false);

    uint_t elementsPerCell = 0;

    switch( dir )
    {
        case stencil::SW:
            elementsPerCell = 1;
            break;
        
        case stencil::W:
            elementsPerCell = 3;
            break;
        
        case stencil::NW:
            elementsPerCell = 1;
            break;
        
        case stencil::S:
            elementsPerCell = 3;
            break;
        
        case stencil::N:
            elementsPerCell = 3;
            break;
        
        case stencil::SE:
            elementsPerCell = 1;
            break;
        
        case stencil::E:
            elementsPerCell = 3;
            break;
        
        case stencil::NE:
            elementsPerCell = 1;
            break;
        
        default:
            elementsPerCell = 0;
    }
    return ci.numCells() * elementsPerCell * sizeof( double );
}



} // namespace pystencils
} // namespace walberla