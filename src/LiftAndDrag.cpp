#include "LiftAndDrag.hpp"

namespace walberla {
namespace lbm {

void LiftAndDrag::initialize() 
{
    real_t AD( real_t( 0.0 ) ), AL( real_t( 0.0 ) );
    boundaryLinks_.clear();

    for( BlockIter_T blockIter = blocks_->begin(); blockIter != blocks_->end(); ++blockIter )
    {
        const FlagField_T * const flagField = blockIter->getData< FlagField_T >( flagFieldId_ );
        
        flag_t fluid    = flagField->getFlag( fluid_ );
        flag_t obstacle = flagField->getFlag( obstacle_ );

        auto xyz = flagField->xyzSize();

        for( cell_idx_t z = xyz.zMin(); z <= xyz.zMax(); ++z )
        {
            for( cell_idx_t y = xyz.yMin(); y <= xyz.yMax(); ++y )
            {
                for( cell_idx_t x = xyz.xMin(); x <= xyz.xMax(); ++x )
                {
                    if( flagField->isFlagSet( x, y, z, fluid ) )
                    {
                        for( auto it = Stencil_T::beginNoCenter(); it != Stencil_T::end(); ++it )
                        {
                            cell_idx_t nx = x + cell_idx_c( it.cx() );
                            cell_idx_t ny = y + cell_idx_c( it.cy() );
                            cell_idx_t nz = 0;

                            if( flagField->isFlagSet( nx, ny, nz, obstacle ) )
                            {
                                boundaryLinks_[ blockIter.get() ].push_back( std::make_pair( Cell( x, y, z ), *it ) );

                                if( it.cx() == 1 && it.cy() == 0 )
                                {
                                    AD += 1;
                                }
                                else if( it.cx() == 0 && it.cy() == 1 )
                                {
                                    AL += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    AD_ = AD;
    AL_ = AL;
    scalingCD_ = real_t( 2.0 ) / ( uMean_ * uMean_ * AD_ );
    scalingCL_ = real_t( 2.0 ) / ( uMean_ * uMean_ * AL_ );

    std::cout << "AD_ = " << AD_ <<"\tAL_ = " << AL_ << std::endl;
}

void LiftAndDrag::run( IBlock * block )
{
    // reset force vector:
    force_[ 0 ] = real_t( 0.0 );
    force_[ 1 ] = real_t( 0.0 );
    force_[ 2 ] = real_t( 0.0 );

    if( boundaryLinks_.find( block ) != boundaryLinks_.end() )
    {
        const PdfField_T * const pdfField = block->getData< PdfField_T >( pdfFieldId_ );

        const std::vector< std::pair< Cell, stencil::Direction > > & boundaryLinks = boundaryLinks_[ block ];

        for( auto pair = boundaryLinks.begin(); pair != boundaryLinks.end(); ++pair )
        {
            const Cell cell( pair->first );
            const stencil::Direction direction( pair->second );

            const real_t boundaryValue = pdfField->get( cell.x() + stencil::cx[ direction ], 
                                                        cell.y() + stencil::cy[ direction ], 
                                                        cell.z() + stencil::cz[ direction ], 
                                                        Stencil_T::idx[ stencil::inverseDir[ direction ] ] );
            
            const real_t fluidValue = pdfField->get( cell.x(), 
                                                     cell.y(), 
                                                     cell.z(), 
                                                     Stencil_T::idx[ direction ] );

            const real_t f = boundaryValue + fluidValue;

            force_[ 0 ] += real_c( stencil::cx[ direction ] ) * f;
            force_[ 1 ] += real_c( stencil::cy[ direction ] ) * f;
            force_[ 2 ] += real_c( stencil::cz[ direction ] ) * f;

        }
    }

    CD_ = scalingCD_ * force_[ 0 ];
    CL_ = scalingCL_ * force_[ 1 ];

    std::cout << "step: " << executionCounter_ << ",\tCD = " << CD_ << ",\tCL = " << CL_ << std::endl;
}

} // namespace lbm
} // namespace walberla