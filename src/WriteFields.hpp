#pragma once

#include "core/all.h"
#include "field/all.h"
#include "domain_decomposition/all.h"
#include "stencil/all.h"

#include <filesystem>
#include <string>
#include <fstream>

namespace walberla {
namespace lbm {

template< typename PdfField_T, typename VectorField_T, typename ScalarField_T >
class WriteFields
{
    public:

    WriteFields( shared_ptr< StructuredBlockForest > & blocks, const BlockDataID pdfFieldId, const BlockDataID velocityFieldId, const BlockDataID densityFieldId, std::string pdfFieldString, std::string velocityFieldString, std::string densityFieldString, uint_t timesteps )
        : blocks_( blocks ), pdfFieldId_( pdfFieldId ), velocityFieldId_( velocityFieldId ), densityFieldId_( densityFieldId ), pdfFieldString_( pdfFieldString ), velocityFieldString_( velocityFieldString ), densityFieldString_( densityFieldString ), timesteps_( timesteps ), t_( uint_c( 0 )  ) {}

    ~WriteFields() = default;

    void operator()( IBlock * /* block */ )
    {
        if( ++t_ == ( timesteps_ - 10 ) ) { execute_(); }
    }

    private:

    void execute_();

    shared_ptr< StructuredBlockForest > & blocks_;
    BlockDataID pdfFieldId_, velocityFieldId_, densityFieldId_;
    std::string pdfFieldString_, velocityFieldString_, densityFieldString_;
    uint_t timesteps_, t_;
};

template< typename PdfField_T, typename VectorField_T, typename ScalarField_T >
void WriteFields< PdfField_T, VectorField_T, ScalarField_T >::execute_()
{
    BlockStorage & blockStorage = blocks_->getBlockStorage();

    field::writeToFile< PdfField_T >( pdfFieldString_, blockStorage, pdfFieldId_ );
    field::writeToFile< VectorField_T >( velocityFieldString_, blockStorage, velocityFieldId_ );
    field::writeToFile< ScalarField_T >( densityFieldString_, blockStorage, densityFieldId_ );
}

}
}

 