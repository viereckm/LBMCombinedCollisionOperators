# include <cmath>

#include "core/Macros.h"
#include "core/DataTypes.h"
#include "MixedCodegenCollisionSweep.hpp"

#include "CodegenCollisionKernels.hpp"


namespace walberla {
namespace lbm {

void D2Q9MixedCodegenCollisionSweep::run( IBlock * block )
{
    auto pdfs = block->getData< field::GhostLayerField< double, 9 > >( pdfsID_ );
    auto velocity = block->getData< field::GhostLayerField<double, 2> >( velocityID_ );
    auto density = block->getData< field::GhostLayerField< double, 1 > >( densityID_ );
    field::GhostLayerField< double, 9 > * pdfs_tmp;
    {
        auto it = cache_pdfs_.find( pdfs );
        if( it != cache_pdfs_.end() )
        { 
            pdfs_tmp = *it;
        }
        else
        {
            pdfs_tmp = pdfs->cloneUninitialized();
            cache_pdfs_.insert(pdfs_tmp );
        }
    }

    auto & omega = this->omega_;
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(density->nrOfGhostLayers()))
    double * RESTRICT  _data_density = density->dataAt(-1, -1, 0, 0);
    WALBERLA_ASSERT_EQUAL(density->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs->nrOfGhostLayers()))
    double * RESTRICT const _data_pdfs = pdfs->dataAt(-1, -1, 0, 0);
    WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_tmp->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_tmp = pdfs_tmp->dataAt(-1, -1, 0, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_tmp->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT  _data_velocity = velocity->dataAt(-1, -1, 0, 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(pdfs->xSize()) + 2))
    const int64_t _size_pdfs_0 = int64_t(int64_c(pdfs->xSize()) + 2);
    WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(pdfs->ySize()) + 2))
    const int64_t _size_pdfs_1 = int64_t(int64_c(pdfs->ySize()) + 2);
    WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
    const int64_t _stride_density_0 = int64_t(density->xStride());
    const int64_t _stride_density_1 = int64_t(density->yStride());
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
    const int64_t _stride_pdfs_tmp_0 = int64_t(pdfs_tmp->xStride());
    const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
    const int64_t _stride_pdfs_tmp_2 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
    const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
    const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
    const int64_t _stride_velocity_2 = int64_t(1 * int64_t(velocity->fStride()));

    fluidFilter_( *block );
    filter1_( *block );
    filter2_( *block );
    
    for( int64_t ctr_1 = 1; ctr_1 < _size_pdfs_1 - 1; ctr_1++ ) // iterate over y-dimension
    {
        for( int64_t ctr_0 = 1; ctr_0 < _size_pdfs_0 -1; ctr_0++ ) // iterate over x-direction
        {
            if( fluidFilter_( ctr_0 - 1, ctr_1 - 1, 0 ) )
            {
                if( filter1_( ctr_0 - 1, ctr_1 - 1, 0 ) )
                {
                    operator1_( ctr_0, ctr_1, _data_density, _data_pdfs, _data_pdfs_tmp, _data_velocity, _stride_density_0, _stride_density_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, omega );
                }
                else if ( filter2_( ctr_0 - 1, ctr_1 - 1, 0 ) )
                {
                    operator2_( ctr_0, ctr_1, _data_density, _data_pdfs, _data_pdfs_tmp, _data_velocity, _stride_density_0, _stride_density_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2,  _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, omega );
                }
                else
                {
                    std::cout << "No collision operator specified at ( " << ctr_0 << " , " << ctr_1 << " )" << std::endl;
                }
            }
        }
    }
    pdfs->swapDataPointers( pdfs_tmp );
}

void D2Q9MixedCodegenCollisionSweep::runOnCellInterval( const shared_ptr< StructuredBlockStorage > & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block )
{
    // Set up cell interval for simulation and check if it is part of the currently processed block:
    CellInterval ci = globalCellInterval;
    CellInterval blockBB = blocks->getBlockCellBB( *block );
    blockBB.expand( ghostLayers );
    ci.intersect( blockBB );
    blocks->transformGlobalToBlockLocalCellInterval( ci, *block );
    
    if( ci.empty() ) return;

    // From here: same as 'run' but with cell interval limits as loop limits...

    auto pdfs = block->getData< field::GhostLayerField< double, 9 > >( pdfsID_ );
    auto velocity = block->getData< field::GhostLayerField<double, 2> >( velocityID_ );
    auto density = block->getData< field::GhostLayerField< double, 1 > >( densityID_ );

    field::GhostLayerField<double, 9> * pdfs_tmp;
    {
        // Getting temporary field pdfs_tmp
        auto it = cache_pdfs_.find( pdfs );
        if( it != cache_pdfs_.end() )
        {
            pdfs_tmp = *it;
        }
        else
        {
            pdfs_tmp = pdfs->cloneUninitialized();
            cache_pdfs_.insert(pdfs_tmp);
        }
    }

    auto & omega = this->omega_;
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(density->nrOfGhostLayers()))
    double * RESTRICT  _data_density = density->dataAt(-1, -1, 0, 0);
    WALBERLA_ASSERT_EQUAL(density->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
    double * RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_tmp = pdfs_tmp->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_tmp->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT  _data_velocity = velocity->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 2))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 2);
    WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 2))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 2);
    WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
    const int64_t _stride_density_0 = int64_t(density->xStride());
    const int64_t _stride_density_1 = int64_t(density->yStride());
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(1 * int64_t(pdfs->fStride()));
    const int64_t _stride_pdfs_tmp_0 = int64_t(pdfs_tmp->xStride());
    const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
    const int64_t _stride_pdfs_tmp_2 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
    const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
    const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
    const int64_t _stride_velocity_2 = int64_t(1 * int64_t(velocity->fStride()));

    filter1_( *block );
    filter2_( *block );
    
    for( int64_t ctr_1 = 1; ctr_1 < _size_pdfs_1 - 1; ++ctr_1 ) // loop over y-dimension
    {
        for( int64_t ctr_0 = 1; ctr_0 < _size_pdfs_0 - 1; ++ctr_0 ) // loop over x-dimension
        {
            if( filter1_( ctr_0, ctr_1, 0 ) )
            { 
                operator1_( ctr_0, ctr_1, _data_density, _data_pdfs, _data_pdfs_tmp, _data_velocity, _stride_density_0, _stride_density_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, omega );
            }
            else if ( filter2_( ctr_0, ctr_1, 0 ) )
            {
                operator2_( ctr_0, ctr_1, _data_density, _data_pdfs, _data_pdfs_tmp, _data_velocity, _stride_density_0, _stride_density_1, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2,  _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, omega );
            }
        }
    }
}

} // !namespace lbm
} // !namespace walberla