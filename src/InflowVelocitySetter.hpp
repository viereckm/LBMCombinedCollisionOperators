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
#include <numbers>
#include <functional>

namespace walberla {
namespace lbm {

class InflowVelocitySetter
{
    public:
    
    InflowVelocitySetter( int64_t t_max, real_t uMax, bool unsteady ) : t_( int64_t( 0 ) ), t_max_( t_max ), uMax_( uMax ), unsteady_( unsteady ) {}

    const Vector3< real_t > operator()( const Cell & cell, const shared_ptr< StructuredBlockForest > & forest, IBlock & block );
    
    private:
    int64_t t_, t_max_;
    const Vector3< real_t > uMax_;
    bool unsteady_;
};

} // namespace lbm
} // namespace walberla