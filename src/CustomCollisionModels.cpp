#include "CustomCollisionModels.hpp"

namespace walberla {
namespace lbm {
namespace collision_model {

const real_t SRT_TRT::threeSixteenth = real_t(3) / real_t(16);

SRT_TRT SRT_TRT::constructWithMagicNumber( const real_t _omega, const real_t _magicNumber /* = threeSixteenth */, const uint_t _level /* = uint_t(0) */ )
{
    SRT_TRT srt_trt;
    srt_trt.initWithMagicNumber( _omega, _magicNumber, _level );
    return srt_trt;
}

void SRT_TRT::configure( IBlock & block, StructuredBlockStorage & sbs )
{
    const uint_t _level = level_;
    level_ = sbs.getLevel( block );
    omega_ = levelDependentRelaxationParameter( level_, omega_, _level );
    lambda_d_ = lambda_d( omega_, magicNumber_ );
    viscosity_ = viscosityFromOmega( omega_ );
}

void SRT_TRT::reset( const real_t _omega, const real_t _lambda_d, const uint_t _omega_level /* = uint_t(0) */ )
{
    omega_ = levelDependentRelaxationParameter( level_, _omega, _omega_level );
    magicNumber_ = magicNumber( _omega, _lambda_d );
    lambda_d_ = lambda_d( omega_, magicNumber_ );
    viscosity_ = viscosityFromOmega( omega_ );
}

void SRT_TRT::resetWithMagicNumber( const real_t _omega, const real_t _magicNumber /* = threeSixteenth */, const uint_t _omega_level /* = uint_t(0) */ )
{
    omega_ = levelDependentRelaxationParameter( level_, _omega, _omega_level );
    lambda_d_ = lambda_d( omega_, _magicNumber );
    magicNumber_ = _magicNumber;
    viscosity_ = viscosityFromOmega( omega_ );
}

void SRT_TRT::initWithMagicNumber( const real_t _omega, const real_t _magicNumber, const uint_t _level )
{
    omega_ = _omega;
    lambda_d_ = lambda_d( omega_, _magicNumber );
    magicNumber_ = _magicNumber;
    viscosity_ = viscosityFromOmega( omega_ );
    level_ = _level;
}

} // namespace collision_model
} // namespace lbm
} // namespace walberla