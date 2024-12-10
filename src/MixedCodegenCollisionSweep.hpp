#pragma once

#include "core/DataTypes.h"
#include "core/logging/Logging.h"

#include "field/GhostLayerField.h"
#include "field/SwapableCompare.h"
#include "field/FlagField.h"
#include "field/EvaluationFilter.h"


#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include <set>
#include <functional>

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

class D2Q9MixedCodegenCollisionSweep
{
    public:
    typedef walberla::uint16_t flag_t;
    typedef FlagField< flag_t > FlagField_T;
    typedef typename walberla::field::FlagFieldEvaluationFilter< FlagField_T > Filter_T;
    typedef std::function< void ( 
        int64_t const /* x */, 
        int64_t const /* y */,
        double * RESTRICT const /* density */, 
        double * RESTRICT const /* src */, 
        double * RESTRICT const /* dst */, 
        double * RESTRICT const /* velocity */,
        int64_t const /* density x stride */,
        int64_t const /* density y stride */, 
        int64_t const /* src x stride */, 
        int64_t const /* src y stride */, 
        int64_t const /* src f stride */, 
        int64_t const /* dst x stride */, 
        int64_t const /* dst y stride */, 
        int64_t const /* dst f stride */, 
        int64_t const /* vel x stride */, 
        int64_t const /* vel y stride */, 
        int64_t const /* vel component stride */, 
        double /* omega */ ) > Operator_T; 

    D2Q9MixedCodegenCollisionSweep( BlockDataID densityID, BlockDataID pdfsID, BlockDataID velocityID, double omega, Filter_T fluidFilter, Filter_T filter1, Filter_T filter2, Operator_T operator1, Operator_T operator2 ) : 
        densityID_( densityID ), pdfsID_( pdfsID ), velocityID_( velocityID ), omega_( omega ), fluidFilter_( fluidFilter ), filter1_( filter1 ), filter2_( filter2 ), operator1_( operator1 ), operator2_( operator2 ) {}

    ~D2Q9MixedCodegenCollisionSweep()
    {
        for( auto p: cache_pdfs_ ) { delete p; }
    }

    void run( IBlock * block );

    void operator()( IBlock * block )
    {
        run( block );
    }

    static std::function< void ( IBlock * ) > getSweep( const shared_ptr< D2Q9MixedCodegenCollisionSweep > & kernel )
    {
        return [ kernel ] ( IBlock * b ) { kernel->run( b ); };
    }

    std::function< void( IBlock *) > getSweep()
    {
        return [ this ] ( IBlock * b ) { this->run( b ); };
    }

    // Same functionality but on a cell interval:

    void runOnCellInterval( const shared_ptr< StructuredBlockStorage > & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block );

    static std::function< void ( IBlock * ) > getSweepOnCellInterval( const shared_ptr< D2Q9MixedCodegenCollisionSweep > & kernel, const shared_ptr < StructuredBlockStorage > & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers = 1 )
    {
        return [ kernel, blocks, globalCellInterval, ghostLayers ] ( IBlock * b ) { kernel->runOnCellInterval( blocks, globalCellInterval, ghostLayers, b ); };
    }

    std::function< void ( IBlock * ) > getSweepOnCellInterval( const shared_ptr< StructuredBlockStorage > & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers = 1 )
    {
        return [ this, blocks, globalCellInterval, ghostLayers ] ( IBlock * b ) { this->runOnCellInterval( blocks, globalCellInterval, ghostLayers, b ); };
    }

    void configure( const shared_ptr< StructuredBlockStorage > & blocks, IBlock * block ) {}

    private:
    BlockDataID densityID_;
    BlockDataID pdfsID_;
    BlockDataID velocityID_;
    double omega_;

    Filter_T fluidFilter_;
    Filter_T filter1_;
    Filter_T filter2_;

    Operator_T operator1_;
    Operator_T operator2_;

    std::set< field::GhostLayerField< double, 9 > *, field::SwapableCompare< field::GhostLayerField< double, 9 > * > > cache_pdfs_;
}; // !class MixedCodegenCollisionSweep

} // !namespace lbm
} // !namespace walberla

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif