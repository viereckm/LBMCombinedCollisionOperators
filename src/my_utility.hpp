#pragma once

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "geometry/all.h"
#include "lbm/all.h"
#include "timeloop/all.h"
#include "stencil/D2Q9.h"

#include <cmath>

#include <math.h>

#include <string>
#include <fstream>
#include <cstdio>
#include <filesystem>

namespace walberla {
namespace lbm {

typedef domain_decomposition::StructuredBlockStorage::iterator BlockIter_T;

template< typename VectorField_T >
void initZeroVelocityField( const shared_ptr< StructuredBlockForest > & blocks, const BlockDataID & velocityFieldId )
{
    for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
    {
        auto u = blockIt->getData< VectorField_T >( velocityFieldId );

        for( auto cellIt = u->beginWithGhostLayerXYZ(); cellIt != u->end(); ++cellIt )
        {
            u->get( cellIt.cell(), 0 ) = real_t( 0 );
            u->get( cellIt.cell(), 1 ) = real_t( 0 );
        }
    }
}


template< typename VectorField_T >
void initConstantVelocityField( const shared_ptr< StructuredBlockForest > & blocks, const BlockDataID & velocityFieldId, Vector3< real_t > uInit )
{
    using FieldIterator_T = typename VectorField_T::iterator;

    auto uInData = uInit.data();
    
    for( BlockIter_T iter = blocks->begin(); iter != blocks->end(); ++iter )
    {
        VectorField_T * u = iter->getData< VectorField_T > ( velocityFieldId );
        uint_t dim = u->fSize();

        for( FieldIterator_T fieldIter = u->beginWithGhostLayerXYZ(); fieldIter != u->end(); ++fieldIter )
        {
            for( uint_t i = 0; i < dim; ++i )
            {
                u->get( fieldIter.cell(), i ) = uInData[ i ];
            }
        }
    }
}

template < typename FlagField_T >
void initFlagFieldBox( const shared_ptr< StructuredBlockStorage > & blocks, BlockDataID flagFieldID, FlagUID fluidFlagUID, FlagUID flagID1, FlagUID flagID2, cell_idx_t xMin, cell_idx_t yMin, cell_idx_t zMin, cell_idx_t xMax, cell_idx_t yMax, cell_idx_t zMax )
{
    using flag_t = typename FlagField_T::value_type;
    using FieldIterator_T = typename FlagField_T::iterator;

    cell_idx_t x, y, z;

    for( BlockIter_T iter = blocks->begin(); iter != blocks->end(); ++iter )
    {
        FlagField_T * f = iter->getData< FlagField_T >( flagFieldID );

        flag_t fluidFlag    = f->getFlag( fluidFlagUID );
        flag_t flag1        = f->registerFlag( flagID1 );
        flag_t flag2        = f->registerFlag( flagID2 );

        for( FieldIterator_T fieldIter = f->beginWithGhostLayerXYZ(); fieldIter != f->end(); ++ fieldIter )
        {
            x = fieldIter.cell().x();
            y = fieldIter.cell().y();
            z = fieldIter.cell().z();

            if( f->isPartOfMaskSet( x, y, z, fluidFlag ) )
            {
                if( ( x >= xMin ) && ( x <= xMax ) && ( y >= yMin ) && ( y <= yMax ) && ( z >= zMin ) && ( z <= zMax ) )
                {
                    f->addMask( x, y, z, flag2 );
                }
                else
                {
                    f->addMask( x, y, z, flag1 );
                }
            }
        }
    }
}

template < typename FlagField_T, typename Stencil_T >
void setMomentumExchangeFlag( const shared_ptr< StructuredBlockStorage > & blocks, BlockDataID flagFieldID, FlagUID fluidUID, FlagUID obstacleUID, FlagUID memUID )
{
    using flag_t = typename FlagField_T::value_type;
    using FieldIterator_T = typename FlagField_T::iterator;

    Cell cell;
    Cell neighbour;
    Cell offset;

    for( BlockIter_T iter = blocks->begin(); iter != blocks->end(); ++iter )
    {
        FlagField_T * f = iter->getData< FlagField_T >( flagFieldID );

        flag_t fluidFlag    = f->getFlag( fluidUID );
        flag_t obstacleFlag = f->getFlag( obstacleUID );
        flag_t memFlag      = f->registerFlag( memUID );

        for( FieldIterator_T fieldIter = f->begin(); fieldIter != f->end(); ++fieldIter )
        {
            cell = fieldIter.cell();

            if( f->isPartOfMaskSet( cell, obstacleFlag ) )
            {
                for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir )
                {
                    offset = Cell( dir.cx(), dir.cy(), dir.cz() );
                    neighbour = cell + offset;

                    if( f->isPartOfMaskSet( neighbour, fluidFlag ) && !( f->isPartOfMaskSet( neighbour, memFlag ) ) )
                    {
                        f->addMask( neighbour, memFlag );
                    }
                }
            }
        }
    }
}

cell_idx_t coordinateToGlobalIndex( real_t x, real_t xMin, real_t xMax, cell_idx_t Nx )
{
    return cell_idx_t( Nx * (x - xMin ) / (xMax - xMin ) );
}


auto inflowVelocityCallback = [] ( const Cell & cell, const shared_ptr< StructuredBlockForest > & forest, IBlock & block, Vector3< real_t> uMax )
{
    Cell globalCell;
    CellInterval domain = forest->getDomainCellBB();
    forest->transformBlockLocalToGlobalCell( globalCell, block, cell );

    real_t Ny = real_c( domain.ySize() );
    real_t y = real_c( globalCell[ 1 ] );

    real_t scaling = 4 * uMax[ 0 ] / ( Ny * Ny );

    real_t oscillation = real_c( 1.0 );
 
    real_t uIn = real_c( oscillation * scaling * y * ( Ny - y ) );

    return Vector3< real_t >( uIn, real_c( 0.0 ), real_c( 0.0 ) );
};

void timeIncrement( int64_t * t )
{
    *t += 1;
}

real_t latticeViscosityFromOmega( real_t omega )
{
    return real_c( ( 1.0 / omega  - 0.5 ) / 3.0 ); 
}

real_t omegaFromLatticeViscosity( real_t nu_l )
{
    return real_c( 2.0 / ( 6.0 * nu_l + 1 ) );
}

/**
 * @brief Computes the discrete square side length via Re = uMax_l * D_l / nu_l => D_l = Re * nu_l / uMax_l
 * 
 * @param Re Reynolds number defining flow regime
 * @param uMax_l maximum (lattice) inflow velocity
 * @param nu_l kinematic (lattice) viscosity
 * @return uint_t lattice side length of square (i.e. number of cells per side) (rounded UP to the next integer!)
 */
uint_t discreteObstacleDimensionFromRe( real_t Re, real_t uMax_l, real_t nu_l )
{
    uint_t D_l = uint_c( Re * nu_l / uMax_l ) + 1;
    
    if( ( D_l % 2 ) == 1 ) return D_l + 1;
    else return D_l;
}

/**
 * @brief corrects lattice inflow velocity s.t. Reynolds number is matched with obstacle dimension via Re = uMax_l * D_l / nu_l => uMax_l = Re * nu_l / D_l
 * 
 * @param Re Reynolds number
 * @param nu_l kinematic (lattice) viscosity
 * @param D_l discrete side length of obstacle
 * @return real_t corrected inflow velocity
 */
real_t velocityCorrection( real_t Re, real_t nu_l, uint_t D_l )
{
    return Re * nu_l / real_c( D_l );
}

/**
 * @brief sets discrete obstacle dimension and lattice velocity to match Reynolds number and lattice viscosity
 * 
 * @param Re Reynolds number (defining flow regime)
 * @param nu_l lattice viscosity (defined via simulation parameter omega)
 * @param uMax_l maximum lattice inflow velocity
 * @param D_l discrete obstacle size (zero-initialized)
 */
void matchReynoldsNumber( real_t Re, real_t nu_l, real_t & uMax_l, uint_t & D_l )
{
    D_l = discreteObstacleDimensionFromRe( Re, uMax_l, nu_l );
    uMax_l = velocityCorrection( Re, nu_l, D_l );
}

/**
 * @brief computes cells per spatial direction via L_l = L / C_l
 * 
 * @param L length in physical coordinates
 * @param C_l length conversion factor (obtained from physical and discrete obstacle dimensions)
 * @return uint_t number of cells in spatial direction
 */
uint_t cellsPerDirection( real_t L, real_t C_l )
{
    return uint_c( L / C_l ) + 1;
}

/**
 * @brief computes "bottom left" corner of square obstacle placed 1/3 downstream the channel at 1/2 channel height
 * 
 * @param Nx number of cells in x direction (discrete channel length)
 * @param Ny number of cells in y direction (discrete channel height)
 * @param D_l discrete obstacle side length#include <math.h>

 * @return Vector3< uint_t > < xMin, yMin, zMin > coordinate of bottom left corner
 */
Vector3< uint_t > obstacleBoxMin( uint_t Nx, uint_t Ny, uint_t D_l )
{
    uint_t xMin = uint_c( Nx / 3 );

    uint_t yMin = uint_c( Ny / 2 - D_l / 2 );

    return Vector3< uint_t >( xMin, yMin, uint_c( 0 ) );
}

/**
 * @brief computes "top right" corner of square obstacle placed 1/3 downstream the channel at 1/2 channel height
 * 
 * @param Nx number of cells in x direction (discrete channel length)
 * @param Ny number of cells in y direction (discrete channel height)
 * @param D_l discrete obstacle side length
 * @return Vector3< uint_t > < xMax, yMax, zMax > coordinate of bottom left corner
 */
Vector3< uint_t > obstacleBoxMax( uint_t Nx, uint_t Ny, uint_t D_l )
{
    uint_t xMax = uint_c( Nx / 3 + D_l );

    /**
     * @brief with same reasoning as in obstacleBoxMin, the upper edge of the obstacle must be at y = Ny/2 - 1 + D_l/2
     * 
     */
    uint_t yMax = uint_c( Ny / 2 + D_l / 2 - 1 );

    return Vector3< uint_t >( xMax, yMax, uint_c( 1 ) );
}



/**
 * @brief generates a parameter file to run a simulation with the physical quantities given as arguments 
 * 
 * @param pathToConfigFile path to file ending with file name, i.e. "/< path_to_directory >/config.prm"
 * @param pathToResultDir path to directory results are saved in
 * @param Re Reynolds number
 * @param omega LBM relaxation parameter
 * @param uMax_l maximum (lattice) inflow velocity
 * @param L physical length of channel [m]
 * @param H physical height of channel [m]
 * @param D physical obstacle dimension [m]
 */
void generateParameterFile( std::string pathToConfigFile, std::string pathToResultDir, real_t Re, real_t omega, real_t uMax_l, real_t timeScaling = 5.0, real_t L = 5.0, real_t H = 0.8, real_t D = 0.1 )
{
    // compute simulation parameters from given input:
    real_t nu_l = latticeViscosityFromOmega( omega );
    uint_t D_l( 0 );
    matchReynoldsNumber( Re, nu_l, uMax_l, D_l );

    // Compute length conversion factor C_l via D_l = D / C_l => C_l = D / D_l;
    real_t C_l = D / real_c( D_l );

    // Domain dimensions:
    uint_t Nx = cellsPerDirection( L, C_l );
    uint_t Ny = uint_c( 8 * D_l ); // blockage ratio 1/8 satisfied precisely!
    uint_t Nz( 1 );

    // Obstacle position:
    Vector3 obstacleMin = obstacleBoxMin( Nx, Ny, D_l );
    Vector3 obstacleMax = obstacleBoxMax( Nx, Ny, D_l );

    // Create file:
    std::fstream configFile;
    configFile.open( pathToConfigFile, std::ios::out | std::ios::trunc );
    
    // Write contents:
    configFile  << "DomainSetup\n"
                << "{\n"
                << "    blocks          < 1, 1, 1 >;\n"
                << "    cellsPerBlock   < " << Nx << ", " << Ny << ", " << Nz << " >;\n"
                << "    periodic        < 0, 0, 1 >;\n"
                << "}\n\n"
                << "Parameters\n"
                << "{\n"
                << "    Re          " << Re                 << ";\n"
                << "    omega       " << omega              << ";\n"
                << "    uMax_l      " << uMax_l             << ";\n"
                << "    D_l         " << D_l                << ";\n"
                << "    nu_l        " << nu_l               << ";\n"
                << "    timeScaling " << timeScaling        << ";\n"
                << "\n"
                << "    resultDir   " << pathToResultDir    << ";\n"
                << "\n"
                << "remainingTimeLoggerFrequency    " << 300 << "; // [s]\n"
                << "numberOfSnapshots               " << 400 << ";\n"
                << "}\n\n"
                << "Boundaries\n"
                << "{\n"
                << "    Border { direction N, S;    walldistance -1;    flag NoSlip; }\n"
                << "    Border { direction W;       walldistance -1;    flag InflowFunc; }\n"
                << "    Border { direction E;       walldistance -1;    flag SEO; }\n\n"
                << "    Body\n"
                << "    {\n"
                << "        shape box;\n\n"
                << "        min < " << obstacleMin[ 0 ] << ", " << obstacleMin[ 1 ] << ", " << obstacleMin[ 2 ] << " >;\n"
                << "        max < " << obstacleMax[ 0 ] << ", " << obstacleMax[ 1 ] << ", " << obstacleMax[ 2 ] << " >;\n\n"
                << "        flag Obstacle;\n"
                << "    }\n"
                << "}\n";

    configFile.close(); 

}

} // namespace lbm
} // namespace walberla