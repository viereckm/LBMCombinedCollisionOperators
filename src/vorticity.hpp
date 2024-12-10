#pragma once

#include "core/all.h"
#include "field/all.h"
#include "domain_decomposition/all.h"
#include "stencil/all.h"

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wunused-parameter"
#   pragma GCC diagnostic ignored "-Wreorder"
#endif

namespace walberla {
namespace lbm {

class VorticitySweep
{
    public:

    typedef field::GhostLayerField< double, 2 > VectorField_T;
    typedef field::GhostLayerField< double, 1 > ScalarField_T;

    VorticitySweep() = default;
    VorticitySweep( const BlockDataID velocityFieldId, const BlockDataID vorticityFieldId, const uint_t frequency )
        : velocityFieldId_( velocityFieldId ), vorticityFieldId_( vorticityFieldId ), frequency_( frequency ) { t_ = uint_t( 0 ); }

    ~VorticitySweep() = default;

    void operator()( IBlock * block )
    {
        if( ( t_ % frequency_ ) == 0 ) run( block );
        t_ += 1;
    }
    
    private:

    void run( IBlock * block );

    BlockDataID velocityFieldId_;
    BlockDataID vorticityFieldId_;

    uint_t frequency_, t_;
};

} // namespace lbm
} // namespace walberla