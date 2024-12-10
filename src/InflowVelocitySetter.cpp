#include "InflowVelocitySetter.hpp"

namespace walberla {
namespace lbm {

const Vector3 < real_t > InflowVelocitySetter::operator()( const Cell & cell, const shared_ptr< StructuredBlockForest > & forest, IBlock & block )
{
    Cell globalCell;
    CellInterval domain = forest->getDomainCellBB();
    forest->transformBlockLocalToGlobalCell( globalCell, block, cell );

    real_t Ny = real_c( domain.ySize() );
    real_t y = real_c( globalCell[ 1 ] );

    real_t scaling = 4 * uMax_[ 0 ] / ( Ny * Ny );

    real_t oscillation = unsteady_ ? real_c( std::sin( M_PI * real_c( t_ ) / real_c( t_max_ ) ) ) : real_c( 1.0 );

    real_t uIn = oscillation * scaling * y * ( Ny - y );

    t_++;

    return Vector3< real_t >( uIn, real_c( 0.0 ), real_c( 0.0 )  );
}
} // namespace lbm
} // namespace walberla