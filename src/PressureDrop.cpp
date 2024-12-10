#include "PressureDrop.hpp"

namespace walberla {
namespace lbm {

void PressureDrop::initialize( cell_idx_t x0, cell_idx_t y0, cell_idx_t x1, cell_idx_t y1, BlockDataID flagFieldId, FlagUID obstacleFlagUID )
{
    FlagField_T * f = blocks_->begin()->getData< FlagField_T >( flagFieldId );
    flag_t obstacle = f->getFlag( obstacleFlagUID );

    // Set x0_ to point right in front of cylinder:
    cell_idx_t x_temp = x0;
    while( f->isPartOfMaskSet( x_temp, y0, 0, obstacle ) )
    {
        --x_temp;
    }
    x0_ = x_temp;
    y0_ = y0;

    // Set x1_ to point right behind cylinder
    x_temp = x1;
    while( f-> isPartOfMaskSet( x_temp, y1, 0, obstacle) )
    {
        ++x_temp;
    }
    x1_ = x_temp;
    y1_ = y1;

    std::cout << "Pressure drop is evaluated between points:\n";
    std::cout << "P0 = ( "<< x0_ << ", " << y0_ << ");\tP1 = ( " << x1_ << ", " << y1_ << " )" << std::endl;
}

void PressureDrop::run( IBlock * block )
{
    ScalarField_T * rho = block->getData< ScalarField_T >( densityFieldId_ );

    real_t rho0 = rho->get( x0_, y0_, 0 );
    real_t rho1 = rho->get( x1_, y1_, 0 );

    // In LBM: p( x ) = (1 / C_s² ) * rho( x ) = rho( x ) / 3
    real_t dP_lattice = ( rho0 - rho1 ) / 3;

    real_t dP = C_p_ * dP_lattice;

    std::cout << "t = " << t_ << "\tdP( t ) = " << dP << std::endl;
}

} // namespace lbm
} // namespace walberla