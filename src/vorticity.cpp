#include "vorticity.hpp"

namespace walberla {
namespace lbm {

void VorticitySweep::run( IBlock * block )
{
    VectorField_T * velocity    = block->getData< VectorField_T >( velocityFieldId_ );
    ScalarField_T * vorticity   = block->getData< ScalarField_T >( vorticityFieldId_ );

    // Access data required:
    double * RESTRICT const _data_velocity  = velocity->dataAt( -1, -1, 0, 0 );
    double * RESTRICT _data_vorticity       = vorticity->dataAt( -1, -1, 0, 0 );

    // Get field dimensions:
    const int64_t _size_x = int64_t( int64_c( velocity->xSize() ) + 2 );
    const int64_t _size_y = int64_t( int64_c( velocity->ySize() ) + 2 );

    // Get strides in each dimension ( + components for velocity )
    const int64_t _stride_velocity_x    = int64_t( velocity->xStride() );
    const int64_t _stride_velocity_y    = int64_t( velocity->yStride() );
    const int64_t _stride_velocity_f    = int64_t( velocity->fStride() ); // required to access y-component of velocity

    const int64_t _stride_vorticity_x   = int64_t( vorticity->xStride() );
    const int64_t _stride_vorticity_y   = int64_t( vorticity->yStride() );

    for( int64_t y = 1; y < _size_y; y++ )
    {
        for( int64_t x = 1; x < _size_x; x++ )
        {
            _data_vorticity[ x * _stride_vorticity_x + y * _stride_vorticity_y ] = 
                .5 *    (  _data_velocity[ ( x + 1 ) * _stride_velocity_x + y * _stride_velocity_y + _stride_velocity_f ] - _data_velocity[ ( x - 1) * _stride_velocity_x + y * _stride_velocity_y + _stride_velocity_f ] 
                         - _data_velocity[ x * _stride_velocity_x + ( y + 1) * _stride_vorticity_y ] + _data_velocity[ x * _stride_velocity_x + ( y - 1 ) * _stride_velocity_y ] );
        }
    }
}

}
}