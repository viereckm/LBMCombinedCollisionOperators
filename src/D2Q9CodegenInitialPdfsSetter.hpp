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
//! \\file D2Q9CodegenInitialPdfsSetter.h
//! \\author pystencils
//======================================================================================================================

#pragma once
#include "core/DataTypes.h"
#include "core/logging/Logging.h"

#include "field/GhostLayerField.h"
#include "field/SwapableCompare.h"
#include "field/FlagField.h"
#include "field/EvaluationFilter.h"

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include <set>
#include <functional>


#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wunused-parameter"
#   pragma GCC diagnostic ignored "-Wreorder"
#endif

namespace walberla {
namespace lbm {


class D2Q9CodegenInitialPdfsSetter
{
public:
    typedef walberla::uint16_t flag_t;
    typedef FlagField< flag_t > FlagField_T;
    typedef typename walberla::field::FlagFieldEvaluationFilter< FlagField_T > Filter_T;
    typedef typename std::function< void ( int64_t const /* x */, int64_t const /* y */,double * RESTRICT /* pdf field */, double * RESTRICT const /* velocity field */, int64_t const /* x stride pdfs */, int64_t const /* y stride pdfs */, int64_t const /* f stride pdfs */,  int64_t const /* x stride vel */, int64_t const /* y stride vel */, int64_t const /* component stride vel */, double /* rho_0 */ ) > Setter_T;

    D2Q9CodegenInitialPdfsSetter( BlockDataID pdfsID, BlockDataID velocityID, double rho_0, Filter_T fluidFilter, Filter_T filter1, Filter_T filter2, Setter_T setter1, Setter_T setter2 )
     : pdfsID_( pdfsID ), velocityID_( velocityID ), rho_0_( rho_0 ), fluidFilter_( fluidFilter ), filter1_( filter1 ), filter2_( filter2 ), setter1_( setter1 ), setter2_( setter2 )
    {};

    void run(IBlock * block);

    void runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block);

    void operator() (IBlock * block)
    {
        run(block);
    }
   
    static std::function<void (IBlock *)> getSweep(const shared_ptr<D2Q9CodegenInitialPdfsSetter> & kernel)
    {
        return [ kernel ] ( IBlock * b ) { kernel->run(b); };
    }

    static std::function<void (IBlock*)> getSweepOnCellInterval(const shared_ptr<D2Q9CodegenInitialPdfsSetter> & kernel, const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers=1)
    {
        return [ kernel, blocks, globalCellInterval, ghostLayers ] ( IBlock * b ) { kernel->runOnCellInterval( blocks, globalCellInterval, ghostLayers, b ); };
    }

    std::function<void (IBlock *)> getSweep()
    {
        return [ this ] ( IBlock * b ) { this->run( b ); };
    }

    std::function<void (IBlock *)> getSweepOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers=1)
    {
        return [ this, blocks, globalCellInterval, ghostLayers ] ( IBlock * b ) { this->runOnCellInterval( blocks, globalCellInterval, ghostLayers, b ); };
    }

    void configure( const shared_ptr<StructuredBlockStorage> & blocks, IBlock * block ){}
   

    private:
   
    BlockDataID pdfsID_;
    BlockDataID velocityID_;
    double rho_0_;

    Filter_T fluidFilter_;
    Filter_T filter1_;
    Filter_T filter2_;

    Setter_T setter1_;
    Setter_T setter2_;

};


} // namespace lbm
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif