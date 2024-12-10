#pragma once

#include "core/all.h"
#include "field/all.h"
#include "domain_decomposition/all.h"
#include "stencil/all.h"

#include <vector>
#include <utility>
#include <map>
#include <iostream>


namespace walberla{ 
namespace lbm {

class LiftAndDrag
{
    public:
        typedef walberla::uint16_t flag_t;
        typedef FlagField< flag_t > FlagField_T;
        typedef field::GhostLayerField< double, 9 > PdfField_T;

//        typedef FlagField_T::iterator FieldIterator_T;
        typedef domain_decomposition::StructuredBlockStorage::iterator BlockIter_T;

        typedef stencil::D2Q9 Stencil_T;

        LiftAndDrag( const shared_ptr < StructuredBlockForest > & blocks,
                     const BlockDataID pdfFieldId, const BlockDataID flagFieldId,
                     const FlagUID fluid, const FlagUID obstacle,
                     real_t uMean, 
                     uint_t checkFrequency ) : 
            blocks_( blocks ),
            pdfFieldId_( pdfFieldId ), flagFieldId_( flagFieldId ),
            fluid_( fluid ), obstacle_( obstacle ),
            force_( Vector3< real_t >( 0.0, 0.0, 0.0 ) ),
            AD_( real_t( 0.0 ) ), AL_( real_t( 0.0 ) ), CD_( real_t( 0.0 ) ), CL_( real_t( 0.0 ) ), uMean_( uMean ), 
            executionCounter_( uint_t( 0 ) ), checkFrequency_( checkFrequency ) { initialize(); }

        ~LiftAndDrag() = default;

        void initialize();
        void run( IBlock * block );

        void operator()( IBlock * block )
        {
            if( ( executionCounter_ % checkFrequency_ ) == 0 ) run( block );
            ++executionCounter_;
        } 

    private:
        const shared_ptr< StructuredBlockForest > & blocks_;
        const BlockDataID pdfFieldId_, flagFieldId_;
        const FlagUID fluid_, obstacle_;
        std::map< IBlock *, std::vector< std::pair< Cell, stencil::Direction > > > boundaryLinks_;
        Vector3< real_t > force_;
        real_t AD_, AL_, CD_, CL_, uMean_, scalingCD_, scalingCL_;
        uint_t  executionCounter_, checkFrequency_;
        
};

} // namespace lbm
} // namesapce walberla