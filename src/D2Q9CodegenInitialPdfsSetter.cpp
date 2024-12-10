#include "core/Macros.h"
#include "D2Q9CodegenInitialPdfsSetter.hpp"

#include <cmath>

namespace walberla {
namespace lbm {

void D2Q9CodegenInitialPdfsSetter::run( IBlock * block )
{
    auto pdfs = block->getData< field::GhostLayerField< double, 9 > >( pdfsID_ );
    auto velocity = block->getData< field::GhostLayerField< double, 2 > >( velocityID_ );

    auto & rho_0 = this->rho_0_;
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT const _data_velocity = velocity->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(pdfs->xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(pdfs->xSize()) + 0);
    WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(pdfs->ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(pdfs->ySize()) + 0);
    WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
    const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
    const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
    const int64_t _stride_velocity_2 = int64_t(1 * int64_t(velocity->fStride()));

    fluidFilter_( *block );    
    filter1_( *block );
    filter2_( *block );

    for( int64_t ctr_1 = 1; ctr_1 < _size_pdfs_1; ++ctr_1 )
    {
        for( int64_t ctr_0 = 1; ctr_0 < _size_pdfs_0; ++ ctr_0 )
        {
            if( fluidFilter_ ( ctr_0, ctr_1, 0 ) )
            {
                if( filter1_( ctr_0, ctr_1, 0 ) )
                {
                    setter1_( ctr_0, ctr_1, _data_pdfs, _data_velocity, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, rho_0 );
                }
                else if ( filter2_( ctr_0, ctr_1, 0 ) )
                {
                    setter2_( ctr_0, ctr_1, _data_pdfs, _data_velocity, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, rho_0 );
                }
            }
        }
    }
}

void D2Q9CodegenInitialPdfsSetter::runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block)
{
   
    CellInterval ci = globalCellInterval;
    CellInterval blockBB = blocks->getBlockCellBB( *block);
    blockBB.expand( ghostLayers );
    ci.intersect( blockBB );
    blocks->transformGlobalToBlockLocalCellInterval( ci, *block );
    if( ci.empty() )
        return;

    auto pdfs = block->getData< field::GhostLayerField<double, 9> >(pdfsID_);
    auto velocity = block->getData< field::GhostLayerField<double, 2> >(velocityID_);

    auto & rho_0 = this->rho_0_;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT const _data_velocity = velocity->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
    const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
    const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
    const int64_t _stride_velocity_2 = int64_t(1 * int64_t(velocity->fStride()));
    
    filter1_( *block );
    filter2_( *block );

    for( int64_t ctr_1 = 1; ctr_1 < _size_pdfs_1; +ctr_1 )
    {
        for( int64_t ctr_0 = 1; ctr_0 < _size_pdfs_0; ++ ctr_0 )
        {
            if( filter1_( ctr_0, ctr_1, 0 ) )
            {
                setter1_( ctr_0, ctr_1, _data_pdfs, _data_velocity, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, rho_0 );
            }
            else if ( filter2_( ctr_0, ctr_1, 0 ) )
            {
                setter2_( ctr_0, ctr_1, _data_pdfs, _data_velocity, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, rho_0 );
            }
        }
    }

}


} // !namespace lbm
} // !namespace walberla