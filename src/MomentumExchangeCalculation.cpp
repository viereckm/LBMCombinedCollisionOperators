#include "MomentumExchangeCalculation.hpp"

namespace walberla {
namespace lbm {

uint_t D2Q9MomentumExchangeCalculation::initialize( const BlockDataID flagFieldID, const FlagUID obstacleUID, const FlagUID memUID )
{
    Cell cell;
    Cell globalCell;
    Cell offset;
    Cell neighbour;

    for( BlockIter_T blockIter = blocks_->begin(); blockIter != blocks_->end(); ++blockIter )
    {
        IBlock & block = *( blockIter );
        FlagField_T * f = blockIter->getData< FlagField_T >( flagFieldID );
        
        flag_t obstacleFlag = f->getFlag( obstacleUID );
        flag_t memFlag      = f->getFlag( memUID );

        for( FieldIterator_T fieldIter = f->beginWithGhostLayerXYZ(); fieldIter != f->end(); ++fieldIter )
        {
            cell = fieldIter.cell();
            if( f->isPartOfMaskSet( cell, memFlag ) )
            {
                blocks_->transformBlockLocalToGlobalCell( globalCell, block, cell );
                memCells_.push_back( globalCell );

                std::vector< Stencil_T::iterator > cellBoundaryLinks;

                for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir )
                {
                    offset = Cell( dir.cx(), dir.cy(), dir.cz() );
                    neighbour = cell + offset;

                    if( f->isPartOfMaskSet( neighbour, obstacleFlag ) )
                    {
                        cellBoundaryLinks.push_back( dir );
                    }
                }

                boundaryLinks_.push_back( cellBoundaryLinks );
            }
        }
    }

    if( /* writeToFile_ */ path_ != "None" )
    {
//        path_ = "./" + path_;
        C_D_path_ = path_ + "/C_D.csv";
        C_L_path_ = path_ + "/C_L.csv";
        time_path_ = path_ + "/time.csv";

        // Clear files
        std::fstream file;
        file.open( C_D_path_, std::ios::out | std::ios::trunc );
        file.close();

        file.open( C_L_path_, std::ios::out | std::ios::trunc );
        file.close();
        
        file.open( time_path_, std::ios::out | std::ios::trunc );
        file.close();
    }

    if( memCells_.size() == boundaryLinks_.size() )
    {
        numMemCells_ = memCells_.size();
        return 0;
    }
    
    return 1;
}

void D2Q9MomentumExchangeCalculation::run( IBlock * block )
{
    Cell globalCell;
    Cell fluidCell;
    Cell offset;
    Cell obstacleCell;

    real_t P_x_temp( 0.0 ), P_y_temp( 0.0 ), f_diff( 0.0 ) , rho( 1.0 ), rho_inv( 0.0 );

    PdfField_T * f = block->getData< PdfField_T >( pdfFieldId_ );

    for( uint_t i = 0; i < numMemCells_; ++i )
    {
        globalCell = memCells_.at( i );
        blocks_->transformGlobalToBlockLocalCell( fluidCell, *block, globalCell );

        rho = 1.0;

        for( Stencil_T::iterator dir = Stencil_T::begin(); dir != Stencil_T::end(); ++dir )
        {
            rho += f->get( fluidCell, Stencil_T::idx[ dir.direction() ] );
        }

        rho_inv = 1.0 / rho;
        
        for( Stencil_T::iterator dir : boundaryLinks_.at( i ) )
        {
            offset = Cell( dir.cx(), dir.cy(), dir.cz() );
            obstacleCell = fluidCell + offset;
            
            real_t f_fluid      = f->get( fluidCell, Stencil_T::idx[ dir.direction() ] );
            real_t f_obstacle   = f->get( obstacleCell, Stencil_T::idx[ dir.inverseDir() ] );
                                   
            f_diff              = 2 * ( f_fluid ); // - f_obstacle );

            P_x_temp            += rho_inv * f_diff * dir.cx();
            P_y_temp            += rho_inv * f_diff * dir.cy();
        }
    }

    F_D_ =  P_x_temp;
    F_L_ =  P_y_temp;

    C_D_ = F_to_C_scaling_ * F_D_;
    C_L_ = F_to_C_scaling_ * F_L_;

//    std::cout << "t = " << t_ << "\tC_D = " << C_D_ <<"\tC_L = " << C_L_ << std::endl;

    if( /* writeToFile_ */ path_ != "None" )
    {
        std::fstream cdfile;
        cdfile.open( C_D_path_, std::ios::out | std::ios::app );
        cdfile << C_D_ << ",";
        cdfile.close();
        
        std::fstream clfile;
        clfile.open( C_L_path_, std::ios::out | std::ios::app );
        clfile << C_L_ << ",";
        clfile.close();

        std::fstream tfile;
        tfile.open( time_path_, std::ios::out | std::ios::app );
        tfile << t_ << ",";
        tfile.close();
    }
}

} // namespace lbm
} // namespace walberla