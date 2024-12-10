/**
 * @file MixedCollisionSweep.h
 * @author Markus Viereck
 * @brief Definition of sweep class using different collision operators
 * @version 0.1
 * @date 2024-07-01
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#pragma once

#include "lbm/all.h"
#include "field/all.h"

#include <functional>
#include <iostream>

namespace walberla {
namespace lbm {

template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T >
class MixedCollisionSweep : public SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >
{
    public:
        typedef typename SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >::PdfField_T PdfField_T;
        typedef typename LatticeModel_T::Stencil Stencil_T;
        typedef typename std::function< void ( const LatticeModel_T & /* lm */, PdfField_T * /* src */, cell_idx_t /* x*/, cell_idx_t /* y*/, cell_idx_t /* z */, real_t /* rho */, Vector3< real_t > /* velocity */, DensityVelocityOut_T /* densityVelocityOut */ ) > CollisionOperator_T; 

        MixedCollisionSweep( const BlockDataID & pdfFieldId,
                             const Filter_T & _filter /* = walberla::field::DefaultEvaluationFilter() */, 
                             const DensityVelocityIn_T & _densityVelocityIn /* = DefaultDensityEquilibriumVelocityCalculation() */, 
                             const DensityVelocityOut_T & _densityVelocityOut /* = DefaultDensityVelocityCallback() */,
                             const Filter_T _filter1,
                             const Filter_T _filter2,
                             const CollisionOperator_T _operator1,
                             const CollisionOperator_T _operator2 ) :
            SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >( pdfFieldId, _filter, _densityVelocityIn, _densityVelocityOut ), filter1_( _filter1 ), filter2_( _filter2 ), operator1_( _operator1 ), operator2_( _operator2 ) {}

        MixedCollisionSweep( const BlockDataID & src, const BlockDataID & dst, 
                             const Filter_T & _filter /* = walberla::field::DefaultEvaluationFilter() */, 
                             const DensityVelocityIn_T & _densityVelocityIn /* = DefaultDensityEquilibriumVelocityCalculation() */, 
                             const DensityVelocityOut_T & _densityVelocityOut /* = DefaultDensityVelocityCallback() */,
                             const Filter_T _filter1,
                             const Filter_T _filter2,
                             const CollisionOperator_T _operator1,
                             const CollisionOperator_T _operator2 ) :
            SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >( src, dst, _filter, _densityVelocityIn, _densityVelocityOut ), filter1_( _filter1 ), filter2_( _filter2 ), operator1_( _operator1 ), operator2_( _operator2 ) {}

        /**
         * @brief method that is executed in every step of simulation timeloop
         * 
         * @param block currently processed block
         * @param numberOfGhostLayersToInclude number of ghost layers required for stencil in use
         */
        void operator()( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) )
        {
            streamCollide( block, numberOfGhostLayersToInclude );
        }

        void streamCollide( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) );
        void stream( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) );
        void collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) );       

    private:
        Filter_T filter1_;
        Filter_T filter2_;

        CollisionOperator_T operator1_;
        CollisionOperator_T operator2_;
}; // MixedCollisionSweep

template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T >
void MixedCollisionSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >::streamCollide( IBlock * const block, const uint_t numberOfGhostLayersToInclude /* = uint_t(0) */ )
{
    PdfField_T * src( NULL );
    PdfField_T * dst( NULL );
    
    // get block-local simulation data
    this->getFields( block, src, dst );

    WALBERLA_ASSERT_GREATER( src->nrOfGhostLayers(), numberOfGhostLayersToInclude );
    WALBERLA_ASSERT_GREATER_EQUAL( dst->nrOfGhostLayers(), numberOfGhostLayersToInclude );

    // required so that density and velocity calculations can be called on dst
    const LatticeModel_T & lm = src->latticeModel();
    dst->resetLatticeModel( lm );

    this->filter( *block );
    this->densityVelocityIn( *block );
    this->densityVelocityOut( *block );

    filter1_( *block );
    filter2_( *block );

    WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,
        {
            static bool print_z = true;
            if( print_z )
            {
                std::cout << "z_coordinate = " << z << std::endl;
                print_z = false;
            }
        } 
        if( this->filter( x, y, z ) )
        {
            for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
            {
                // streaming step
                dst->get( x, y, z, d.toIdx() ) = src->get( x-d.cx(), y - d.cy(), z-d.cz(), d.toIdx() );
            }

            Vector3< real_t > velocity;
            real_t rho = this->densityVelocityIn( velocity, dst, x, y, z );
            //this->densityVelocityOut( x, y, z, lm, velocity, rho );

            if( filter1_( x, y, z ) )
            {
                operator1_( lm, dst, x, y, z, rho, velocity, this->densityVelocityOut_ );
            }
            else if( filter2_( x, y, z ) )
            {
                operator2_( lm, dst, x, y, z, rho, velocity, this->densityVelocityOut_ );
            }
            else
            { 
                std::cout << "No collision operator specified at ( " << x << ", " << y << ", "  << z << " )" << std::endl;
            }
        }
    ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ

    src->swapDataPointers( dst );
}

template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T >
void MixedCollisionSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >::stream( IBlock * const block, const uint_t numberOfGhostLayersToInclude /* = uint_t(0) */ )
{
    PdfField_T * src( NULL );
    PdfField_T * dst( NULL );
    this->getFields( block, src, dst );
    StreamPull< LatticeModel_T >::execute( src, dst, block, this->filter_, numberOfGhostLayersToInclude );
}

template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T >
void MixedCollisionSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >::collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude /* = uint_t(0) */ )
{
    PdfField_T * src = this->getSrcField( block );
    WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers(), numberOfGhostLayersToInclude );

    const auto & lm = src->latticeModel();

    this->filter( *block );
    this->densityVelocityIn( *block );
    this->densityVelocityOut( *block );

    filter1_( *block );
    filter2_( *block );

    WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,
        if( this->filter( x, y, z ) )
        {
            Vector3< real_t > velocity;
            real_t rho = this->densityVelocityIn( velocity, src, x, y, z );
            this->densityVelocityOut( x, y, z, lm, velocity, rho );

            if( filter1_( x, y, z ) )
            {
                operator1_( lm, src, x, y, z );
            }
            else if( filter2_( x, y, z ) )
            {
                operator2_( lm, src, x, y, z );
            }
            else
            { 
                std::cout << "No collision operator specified at ( " << x << ", " << y << ", "  << z << " )" << std::endl;
            }
        }
    ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}


template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T, typename CollisionOperator_T >
shared_ptr< MixedCollisionSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T > >
makeMixedCollisionSweep( const BlockDataID & pdfFieldId, 
                         const BlockDataID & flagFieldId, 
                         const FlagUID & fluidFlag, 
                         const CollisionOperator_T operator1,
                         const CollisionOperator_T operator2,
                         const FlagUID & flag1,
                         const FlagUID & flag2 )
{
    using MCS_T = MixedCollisionSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >;
    walberla::Set< FlagUID > cellsToEvaluate( fluidFlag );
    walberla::Set< FlagUID > cellsToEvaluate1( flag1 );
    walberla::Set< FlagUID > cellsToEvaluate2( flag2 );

    return shared_ptr< MCS_T > ( new MCS_T( pdfFieldId, 
                                            Filter_T( flagFieldId, cellsToEvaluate ),
                                            DefaultDensityEquilibriumVelocityCalculation(),
                                            DefaultDensityVelocityCallback(),
                                            Filter_T( flagFieldId, cellsToEvaluate1 ),
                                            Filter_T( flagFieldId, cellsToEvaluate2 ),
                                            operator1,
                                            operator2 ) );
}
} // namespace lbm
} // namespace walberla