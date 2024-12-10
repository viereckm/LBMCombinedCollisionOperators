#pragma once

#include "core/all.h"
#include "field/all.h"
#include "domain_decomposition/all.h"

#include <iostream>

#include <cmath>

namespace walberla {
namespace lbm {

class PressureDrop
{
    public:
    
    typedef walberla::uint16_t flag_t;
    typedef FlagField< flag_t > FlagField_T;
    typedef field::GhostLayerField< double, 1 > ScalarField_T;

    PressureDrop( const shared_ptr< StructuredBlockForest > & blocks, const BlockDataID densityFieldId, const BlockDataID flagFieldId, const FlagUID obstacleFlagUID, const real_t C_p, const uint_t frequency, cell_idx_t x0, cell_idx_t y0, cell_idx_t x1, cell_idx_t y1 )
        : blocks_( blocks ), densityFieldId_( densityFieldId ), C_p_( C_p ), t_( uint_t( 0 ) ), frequency_( frequency )
    {
        initialize( x0, y0, x1, y1, flagFieldId, obstacleFlagUID );
    }

    void operator()( IBlock * block )
    {
        if( ( t_ % frequency_ ) == 0 ) run( block );
        ++t_;
    }

    private:

    void initialize( cell_idx_t x0, cell_idx_t y0, cell_idx_t x1, cell_idx_t y1, BlockDataID flagFieldId, FlagUID obstacleFlagUID);
    void run( IBlock * block );

    const shared_ptr< StructuredBlockForest > & blocks_;
    BlockDataID densityFieldId_;

    real_t C_p_;
    uint_t t_, frequency_;
    cell_idx_t x0_, y0_, x1_, y1_;

}; // class PressureDrop

} // namespace lbm
} // namespace walberla