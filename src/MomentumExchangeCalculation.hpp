#pragma once

#include "core/all.h"
#include "field/all.h"
#include "domain_decomposition/all.h"
#include "stencil/all.h"

#include <vector>
#include <iostream>
#include <filesystem>
#include <string>
#include <fstream>

namespace walberla {
namespace lbm {

class D2Q9MomentumExchangeCalculation
{
    public:
    
    typedef walberla::uint16_t flag_t;
    typedef FlagField< flag_t > FlagField_T;
    typedef field::GhostLayerField< double, 9 > PdfField_T;

    typedef FlagField_T::iterator FieldIterator_T;
    typedef domain_decomposition::StructuredBlockStorage::iterator BlockIter_T;
    
    typedef stencil::D2Q9 Stencil_T;

    D2Q9MomentumExchangeCalculation() = default;
    
    D2Q9MomentumExchangeCalculation( const shared_ptr< StructuredBlockForest > & blocks, const BlockDataID pdfFieldID, const BlockDataID flagFieldID, const FlagUID obstacleUID, const FlagUID memUID, /* real_t rho, */ real_t uMean, real_t D, uint_t frequency, real_t dx, real_t dt, /* bool writeToFile ,*/ std::string path = "None" )
        : blocks_( blocks), pdfFieldId_( pdfFieldID ),  F_D_( 0.0 ), F_L_( 0.0 ), C_D_( 0.0 ), C_L_( 0.0 ), F_to_C_scaling_( 2.0 / ( /* rho * */ uMean * uMean * D ) ), Momentum_to_Force_scaling_( dx * dx / dt ), t_( 0 ), frequency_( frequency ), /* writeToFile_( writeToFile ) ,*/ path_( path )
    {
        error_ = initialize( flagFieldID, obstacleUID, memUID );
    }

    ~D2Q9MomentumExchangeCalculation() = default;
    
    uint_t initialize( const BlockDataID flagFieldID, const FlagUID obstacleUID, const FlagUID memUID );

    void run( IBlock *  block );

    void operator()( IBlock * block )
    {
        if( (t_ % frequency_ ) == 0 ) run( block );
        t_+= 1;
    }

    real_t F_D() { return F_D_; }
    real_t F_L() { return F_L_; }
    real_t C_D() { return C_D_; }
    real_t C_L() { return C_L_; }
    uint_t numMemCells() { return numMemCells_; }
    uint_t ERROR() { return error_; }

    private:
    const shared_ptr< StructuredBlockForest > & blocks_;
    std::vector< Cell > memCells_;
    std::vector< std::vector< Stencil_T::iterator > > boundaryLinks_;
    BlockDataID pdfFieldId_;
    real_t F_D_, F_L_, C_D_, C_L_, F_to_C_scaling_, Momentum_to_Force_scaling_;
    uint_t numMemCells_;
    uint_t error_, t_, frequency_;

//    bool writeToFile_;

    std::string path_;    
    std::string C_D_path_;
    std::string C_L_path_;
    std::string time_path_; 
};

} // namespace lbm
} // namespace walberla