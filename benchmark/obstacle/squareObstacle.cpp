#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "geometry/all.h"
#include "lbm/all.h"
#include "timeloop/all.h"
#include "stencil/D2Q9.h"

#include "../../src/all.h"

#include <functional>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>

#include <cstdio>

namespace walberla {

// Typedefs
typedef pystencils::D2Q9CodegenPackInfo                             PackInfo_T;
typedef stencil::D2Q9                                               Stencil_T;

typedef field::GhostLayerField< real_t, Stencil_T::Size >           PdfField_T;
typedef field::GhostLayerField< real_t, Stencil_T::D >              VectorField_T;
typedef field::GhostLayerField< real_t, 1 >                         ScalarField_T;

typedef uint16_t                                                    flag_t;
typedef FlagField< flag_t >                                         FlagField_T;

typedef lbm::D2Q9CodegenNoSlip                                      NoSlip_T;
typedef lbm::D2Q9CodegenUBB                                         UBB_T;
typedef lbm::D2Q9CodegenSEO_E                                       SEO_T;
typedef lbm::D2Q9InflowFunc                                         Inflow_T;

typedef typename lbm::D2Q9MixedCodegenCollisionSweep::Filter_T      Filter_T;
typedef typename lbm::D2Q9MixedCodegenCollisionSweep::Operator_T    Operator_T;
typedef typename lbm::D2Q9CodegenInitialPdfsSetter::Setter_T        Setter_T;

typedef domain_decomposition::StructuredBlockStorage::iterator      BlockIter_T;

int main( int argc, char** argv )
{
    walberla::Environment walberlaEnv( argc, argv );

    // Create domain from .prm file:
    shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGridFromConfig( walberlaEnv.config() );

    // Domain dimensions in lattice units (required to derive unit conversion factors):
    real_t Nx = real_c( blocks->getNumberOfXCellsPerBlock() ); // Only one block in this simulation
    real_t Ny = real_c( blocks->getNumberOfYCellsPerBlock() ); // Only one block in this simulation

    // Get Reynolds number (defines flow regime and couples D, u and nu):
    auto parameters = walberlaEnv.config()->getOneBlock( "Parameters" );
    real_t Re = parameters.getParameter< real_t >( "Re", real_c( -1.0 ) );

    if( Re < 0 )
    {
        std::cout << "Reynolds number Re must be specified in section Parameters of .prm file" << std::endl;
        return -1;
    }

    // Max inflow velocity in lattice units:
    real_t uMax_l = parameters.getParameter< real_t >( "uMax_l", real_c( 0.1 ) );
    real_t uMean_l = real_c( 2.0 / 3.0 ) * uMax_l;

    // Inflow velocity as vector:
    Vector3< real_t > uIn_l( uMax_l, real_c( 0.0 ), real_c( 0.0 ) );

    // (Zero) Initial velocity:
    Vector3<real_t > u0_l( real_c( 0.0 ), real_c( 0.0 ), real_c( 0.0 ) ); 

    // Physical domain size and simulation parameters (required for unit conversions):
    real_t L = 5.0;     // [m] - length of channel
    real_t H = 0.8;     // [m] - height of channel
    real_t D = 0.1;     // [m] - side length of square obstacle
    real_t rho = 1.0;   // [kg / m³] - fluid density
    real_t nu = 0.001;  // [m² / s ] - kinematic viscosity

    // Physical max inflow velocity from Re = uMax * D / nu <-> uMax = Re * nu / D
    real_t uMax = Re * nu / D;
    real_t uMean = real_c( 2.0 / 3.0 ) * uMax;

    // Derive unit conversion factors via: X_l = X / C_X <-> C_X = X / X_l
    // From C_l (length conversion factor) ne can derive D_l, the obstacle side length in lattice units.
    // Then, nu_l is available via Re_l = uMax_l * D_l / nu_l <-> nu_l = uMax_l * D_l / Re_l AND Re_l = Re. 
    real_t C_l  = H / Ny;
    real_t D_l_  = D / C_l;
    real_t D_l = parameters.getParameter< real_t >( "D_l", D_l_ );
    real_t nu_l = parameters.getParameter< real_t >( "nu_l", real_c( uMax_l * D_l / Re ) );
    
    // Simulation parameter:
    real_t omega_ = lbm::omegaFromLatticeViscosity( nu_l );
    real_t omega = parameters.getParameter< real_t >( "omega", omega_);

    real_t C_u  = uMax / uMax_l;
    real_t C_nu = nu / nu_l;

    real_t C_t = C_l / C_u;

    // Compute grid Reynolds number (should not exceed O(10) [Krueger, Timm])
    real_t Re_g = uMax_l / nu_l;

    // How long should the simulation run?
    // Without obstacle, an average particle takes L / uMean [s] to travel the channel.
    // This can be scaled using timeScaling in the Parameters section of the .prm file.
    real_t timeScaling = parameters.getParameter< real_t >( "timeScaling", real_c( 2.0 ) );
    
    real_t t_end = timeScaling * L / uMean;

    // Check: obstacle position:
    cell_idx_t xMinObstacle = cell_idx_c( Nx / 3 );
    cell_idx_t yMinObstacle = cell_idx_c( Ny / 2 - D_l / 2 );

    cell_idx_t xMaxObstacle = cell_idx_c( Nx / 3 ) + cell_idx_c( D_l );
    cell_idx_t yMaxObstacle = cell_idx_c( Ny / 2 + D_l / 2 - 1 );

    // Subdomain Cumulants are used in (same as refinemen region in Guo's paper):
    cell_idx_t xMinCumulants = cell_idx_c( Nx / 3 - D_l );
    cell_idx_t yMinCumulants = cell_idx_c( Ny / 2 -  2 * D_l );

    cell_idx_t xMaxCumulants = cell_idx_c( Nx / 3 + 3 * D_l );
    cell_idx_t yMaxCumulants = cell_idx_c( Ny / 2 +  2 * D_l - 1 );

    std::string strRe( "Re" + std::to_string( uint_c( Re ) ) );
    std::string strNx( "Nx" + std::to_string( uint_c( Nx ) ) );
    std::string strNy( "Ny" + std::to_string( uint_c( Ny ) ) );

    // Path to directory where results are saved:
    std::filesystem::path buildDirPath( std::filesystem::current_path() );
    std::string buildDirString( buildDirPath.string() );
    std::string resultPath_ = "/results/BreuerGuo" + strRe + strNx + strNy + "/";
    std::string resultPath__ = parameters.getParameter< std::string >( "resultDir", resultPath_ );
    std::string resultPath = buildDirString + resultPath__;

    // Path to directory where last snapshot of PDF, velocity and density field are saved:
    // Saving the corect snapshots somehow does not work...
    std::string snapshotLoadPath    = "examples/obstacle/BreuerGuo/" + strRe + "/" + strNx + strNy + "/";
    
    std::string loadPdfField        = snapshotLoadPath + "pdfField";
    std::string loadVelocityField   = snapshotLoadPath + "velocityField";
    std::string loadDensityField    = snapshotLoadPath + "densityField";

    const char * CLoadPdfField      = loadPdfField.c_str();
    const char * CLoadVelocityField = loadVelocityField.c_str();
    const char * CLoadDensityField  = loadDensityField.c_str();
    
    std::string snapshotSafePath    = "../" + snapshotLoadPath;
    
    std::string safePdfField        = snapshotSafePath + "pdfField";
    std::string safeVelocityField   = snapshotSafePath + "velocityField";
    std::string safeDensityField    = snapshotSafePath + "densityField";

    // decide whether or not to load fields from saved snapshots:
    bool loadSnapshots = false;

    std::FILE * pdfSnapshot         = std::fopen( CLoadPdfField, "r" );
    std::FILE * velocitySnapshot    = std::fopen( CLoadVelocityField, "r" );
    std::FILE * densitySnapshot     = std::fopen( CLoadDensityField, "r" );

    if( !pdfSnapshot || !velocitySnapshot || !densitySnapshot )
    {
        std::cout << "Snapshots of fields not found; starting with standard initial values" << std::endl;
    }
    else 
    {
        loadSnapshots = true;
        std::cout << "Starting with last snapshot from previous simulation!" << std::endl;

        t_end = 2 * L / uMean;

        std::fclose( pdfSnapshot );
        std::fclose( velocitySnapshot );
        std::fclose( densitySnapshot );       
    }

    // Set number of timesteps in simulation:
    uint_t timesteps = uint_c( t_end / C_t );

    // Set Frequency to write .vtk files:
    uint_t numberOfSnapshots = parameters.getParameter< uint_t >( "numberOfSnapshots", uint_c( 100 ) );
    uint_t VTKwriteFrequency = std::max( uint_c( 1 ), uint_c( timesteps / numberOfSnapshots ) );

    // Time logging frequency:
    real_t remainingTimeLoggerFrequency = parameters.getParameter< real_t >( "remainingTimeLoggerFrequency", real_c( 10 ) );    

    // Log parameters to console:
    std::cout   << "Simulating flow around a square obstacle in the following setup:"   << std::endl
                << "------------------------------------------------------------------" << std::endl
                << "Geometry:"                                                          << std::endl
                << "L = " << L << "\tH = " << H << "\tD = " << D << "\tC_l = " << C_l   << std::endl
                << "Nx = " << Nx << "\tNy = " << Ny << "\tD_l = " << D_l << "#Cells = " << Nx*Ny << std::endl
                << "BottomLeft = ( " << xMinObstacle << ", " << yMinObstacle << " ),\tTopRight = ( " << xMaxObstacle << ", " << yMaxObstacle << " )" << std::endl
                << "------------------------------------------------------------------" << std::endl
                << "Simulation parameters:"                                             << std::endl
                << "uMax = " << uMax << "\tuMax_l = " << uMax_l << "C_u = " << C_u      << std::endl
                << "nu = " << nu << "\tnu_l = " << nu_l << "\tC_nu = " << C_nu          << std::endl
                << "t_end = " << t_end << "\tNt = " << timesteps << "C_t = " << C_t     << std::endl
                << "omega = " << omega << "\tRe_grid = " << Re_g << "\tRe = " << Re     << std::endl
                << "------------------------------------------------------------------" << std::endl
                << "Cumulants are evaluated in this subdomain:"                         << std::endl
                << "BottomLeft = ( " << xMinCumulants << ", " << yMinCumulants << " ),\tTopRight = ( " << xMaxCumulants << ", " << yMaxCumulants << " )" << std::endl
                << "------------------------------------------------------------------" << std::endl
                << "Results are written to: " << resultPath                             << std::endl
                << "------------------------------------------------------------------" << std::endl;

/*    // Log simulation parameters to file for postprocessing:
    std::fstream file;
    file.open( std::string( resultPath + "/parameters.csv"), std::ios::out | std::ios::trunc );
    file    << "L,"             << L                << "\n"
            << "H,"             << H                << "\n"
            << "D,"             << D                << "\n"
            << "C_l,"           << C_l              << "\n"
            << "Nx,"            << Nx               << "\n"
            << "Ny,"            << Ny               << "\n"
            << "D_l,"           << D_l              << "\n"
            << "xMinObstacle,"  << xMinObstacle     << "\n"
            << "yMinObstacle,"  << yMinObstacle     << "\n"
            << "xMaxObstacle,"  << xMaxObstacle     << "\n"
            << "yMaxObstacle,"  << yMaxObstacle     << "\n"
            << "uMax,"          << uMax             << "\n"
            << "uMax_l,"        << uMax_l           << "\n"
            << "C_u,"           << C_u              << "\n"
            << "nu,"            << nu               << "\n"
            << "nu_l,"          << nu_l             << "\n"
            << "C_nu,"          << C_nu             << "\n"
            << "t_end,"         << t_end            << "\n"
            << "Nt,"            << timesteps        << "\n"
            << "C_t,"           << C_t              << "\n"
            << "omega,"         << omega            << "\n"
            << "Re_grid,"       << Re_g             << "\n"
            << "xMinCumulants," << xMinCumulants    << "\n"
            << "yMinCumulants," << yMinCumulants    << "\n"
            << "xMaxCumulants," << xMaxCumulants    << "\n"
            << "yMaxCumulants," << yMaxCumulants    << "\n";
    file.close();
*/
    // Actual simulation setup starts here
    // Operators and PDF Setter:
    Operator_T  srtKernel = lbm::codegenKernels::D2Q9SRTKernel;
    Setter_T    srtSetter = lbm::codegenKernels::D2Q9SRTPdfsSetter;

    Operator_T  cumulantKernel = lbm::codegenKernels::D2Q9CumulantKernel;
    Setter_T    cumulantSetter = lbm::codegenKernels::D2Q9CumulantPdfsSetter;

    // Register fields:
    BlockDataID velocityFieldId     = field::addToStorage< VectorField_T >( blocks, "velocity", real_c( 0.0 ), field::fzyx );
    BlockDataID flagFieldId         = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field");
    BlockDataID pdfFieldId          = field::addToStorage< PdfField_T >( blocks, "pdf field", real_c( 0.0 ), field::fzyx );
    BlockDataID densityFieldId      = field::addToStorage< ScalarField_T >( blocks, "density", real_c( 1.0 ), field::fzyx );
    BlockDataID vorticityFieldId    = field::addToStorage< ScalarField_T >( blocks, "vorticity", real_c( 0.0 ), field::fzyx );

    BlockStorage & blockStorage = blocks->getBlockStorage();

    if( loadSnapshots )
    {
        field::readFromFile< PdfField_T >( loadPdfField, blockStorage, pdfFieldId );
        field::readFromFile< VectorField_T >( loadVelocityField, blockStorage, velocityFieldId );
        field::readFromFile< ScalarField_T >( loadDensityField, blockStorage, densityFieldId );
    }

    // Initialize flags:
    // Domain:
    const FlagUID fluidFlagUID( "Fluid flag");
    const FlagUID srtFlagUID( "SRT operator flag" );
    const FlagUID cumulantFlagUID( "Cumulant operator flag" );
    const FlagUID memFlagUID( "Momentum exchange cells" );

    // Boundaries:
    const FlagUID noSlipFlagUID( "NoSlip" );
    const FlagUID obstacleFlagUID( "Obstacle" );
    const FlagUID inflowFuncFlagUID( "InflowFunc" );
    const FlagUID SEOFlagUID( "SEO" );

    // Generate bundary handling:
    auto boundariesConfig = walberlaEnv.config()->getOneBlock( "Boundaries" );

    // Inflow function:
    std::function < Vector3< real_t >( const Cell &, const shared_ptr< StructuredBlockForest > &, IBlock & ) >
        uInProfile = std::bind( lbm::inflowVelocityCallback, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, uIn_l );

    // Generate boundary instances:
    NoSlip_T    noSlip( blocks, pdfFieldId );
    NoSlip_T    obstacle( blocks, pdfFieldId );
    Inflow_T    inflow( blocks, pdfFieldId, uInProfile ); 
    SEO_T       seo( blocks, pdfFieldId );
    
    geometry::initBoundaryHandling< FlagField_T >( *blocks, flagFieldId, boundariesConfig );
    geometry::setNonBoundaryCellsToDomain< FlagField_T >( *blocks, flagFieldId, fluidFlagUID );

    noSlip.fillFromFlagField< FlagField_T >( blocks, flagFieldId, noSlipFlagUID, fluidFlagUID );
    obstacle.fillFromFlagField< FlagField_T >( blocks, flagFieldId, obstacleFlagUID, fluidFlagUID );
    inflow.fillFromFlagField< FlagField_T >( blocks, flagFieldId, inflowFuncFlagUID, fluidFlagUID );
    seo.fillFromFlagField< FlagField_T >( blocks, flagFieldId, SEOFlagUID, fluidFlagUID );

    // Set flags for SRT / Cumulant operator
    lbm::initFlagFieldBox< FlagField_T >( blocks, flagFieldId, fluidFlagUID, srtFlagUID, cumulantFlagUID, xMinCumulants, yMinCumulants, 0, xMaxCumulants, yMaxCumulants, 0 );

    // Mark Momentum Exchange cells
    lbm::setMomentumExchangeFlag< FlagField_T, Stencil_T >( blocks, flagFieldId, fluidFlagUID, obstacleFlagUID, memFlagUID );

    // Initialize operator kernel specific filters 
    walberla::Set< FlagUID > fluidCells( fluidFlagUID );
    walberla::Set< FlagUID > srtCells( srtFlagUID );
    walberla::Set< FlagUID > cumulantCells( cumulantFlagUID );

    Filter_T fluidFilter( flagFieldId, fluidCells );
    Filter_T srtFilter( flagFieldId, srtCells );
    Filter_T cumulantFilter( flagFieldId, cumulantCells );

    // Initialize velocity field:
    lbm::initConstantVelocityField< VectorField_T >( blocks, velocityFieldId, u0_l );

    // Initialize PDFs accordingly;
    lbm::D2Q9CodegenInitialPdfsSetter pdfSetter(
        pdfFieldId,
        velocityFieldId,
        rho,
        fluidFilter,
        srtFilter,
        cumulantFilter,
        srtSetter,
        cumulantSetter
    );

    for( BlockIter_T blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
    {
        pdfSetter( &( *blockIt ) );
    }

    lbm::D2Q9MixedCodegenCollisionSweep simulationSweep(
        densityFieldId,
        pdfFieldId,
        velocityFieldId,
        omega,
        fluidFilter,
        srtFilter,
        cumulantFilter,
        srtKernel,
        cumulantKernel
    );

    // 
    lbm::D2Q9MomentumExchangeCalculation memSweep( 
        blocks, 
        pdfFieldId,
        flagFieldId,
        obstacleFlagUID,
        memFlagUID, 
        // rho,
        uMean_l,
        D_l,
        VTKwriteFrequency,
        C_l,
        C_t,
//        true,
        resultPath
    );

    lbm::VorticitySweep vorticitySweep(
        velocityFieldId,
        vorticityFieldId,
        VTKwriteFrequency
    );

    if( memSweep.ERROR() )
    {
        std::cout << "Error initializing MEM cells..." << std::endl;
        return -1;
    }

/*
    lbm::WriteFields< PdfField_T, VectorField_T, ScalarField_T > fieldWriteSweep(
        blocks,
        pdfFieldId,
        velocityFieldId,
        densityFieldId,
        safePdfField,
        safeVelocityField,
        safeDensityField,
        timesteps
    );
*/
    // Initialize communication scheme:
    blockforest::communication::UniformBufferedScheme< Stencil_T > communication( blocks );
    communication.addPackInfo( make_shared< PackInfo_T >( pdfFieldId ) );

    // Set up timeloop:
    SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

    timeloop.addFuncBeforeTimeStep( communication, "communication" );
    timeloop.add() << Sweep( noSlip );
    timeloop.add() << Sweep( obstacle );
    timeloop.add() << Sweep( inflow );
    timeloop.add() << Sweep( seo );
    timeloop.add() << Sweep( memSweep );
    timeloop.add() << Sweep( simulationSweep );
    timeloop.add() << Sweep( vorticitySweep );
//    timeloop.add() << Sweep( fieldWriteSweep );

    // remaining time logger:
    timeloop.addFuncAfterTimeStep( timing::RemainingTimeLogger( timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency ), "remaining time logger" );

    if( VTKwriteFrequency > 0 )
    {
        auto vtkOutput = vtk::createVTKOutput_BlockData( *blocks, "MixedCollisionVelocityField", VTKwriteFrequency, 0, false, resultPath, "simulation_step", false, true, true, false, 0 );
    
        auto velWriter = make_shared< field::VTKWriter< VectorField_T > >( velocityFieldId, "Velocity" );
        vtkOutput->addCellDataWriter( velWriter );
        auto rhoWriter = make_shared< field::VTKWriter< ScalarField_T > >( densityFieldId, "Density" );
        vtkOutput->addCellDataWriter( rhoWriter );
        auto vortWriter = make_shared< field::VTKWriter< ScalarField_T > >( vorticityFieldId, "Vorticity" );
        vtkOutput->addCellDataWriter( vortWriter );
        
        timeloop.addFuncBeforeTimeStep( vtk::writeFiles( vtkOutput), "VTK Output" );
    }

    timeloop.run();

    std::fstream file;
    file.open( std::string( resultPath + "/parameters.csv"), std::ios::out | std::ios::trunc );
    file    << "L,"             << L                << "\n"
            << "H,"             << H                << "\n"
            << "D,"             << D                << "\n"
            << "C_l,"           << C_l              << "\n"
            << "Nx,"            << Nx               << "\n"
            << "Ny,"            << Ny               << "\n"
            << "D_l,"           << D_l              << "\n"
            << "xMinObstacle,"  << xMinObstacle     << "\n"
            << "yMinObstacle,"  << yMinObstacle     << "\n"
            << "xMaxObstacle,"  << xMaxObstacle     << "\n"
            << "yMaxObstacle,"  << yMaxObstacle     << "\n"
            << "uMax,"          << uMax             << "\n"
            << "uMax_l,"        << uMax_l           << "\n"
            << "C_u,"           << C_u              << "\n"
            << "nu,"            << nu               << "\n"
            << "nu_l,"          << nu_l             << "\n"
            << "C_nu,"          << C_nu             << "\n"
            << "t_end,"         << t_end            << "\n"
            << "Nt,"            << timesteps        << "\n"
            << "C_t,"           << C_t              << "\n"
            << "omega,"         << omega            << "\n"
            << "Re_grid,"       << Re_g             << "\n"
            << "xMinCumulants," << xMinCumulants    << "\n"
            << "yMinCumulants," << yMinCumulants    << "\n"
            << "xMaxCumulants," << xMaxCumulants    << "\n"
            << "yMaxCumulants," << yMaxCumulants    << "\n";
    file.close();


    return EXIT_SUCCESS;
}   

} // namespace walberla

int main( int argc, char** argv ) { return walberla::main( argc, argv ); }