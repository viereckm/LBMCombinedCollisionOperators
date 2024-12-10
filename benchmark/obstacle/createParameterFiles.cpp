/**
 * @file  createParameterFiles.cpp
 * @author Markus Viereck
 * @brief function automatically creating .prm files + directory tree for grid study on omega; MUST BE RUN FROM build/ to correctly implement directory structure
 * @version 0.1
 * @date 2024-10-14
 * 
 * @copyright Copyright (c) 2024
 * 
**/

#include "../../src/all.h"

#include <vector>

#include <filesystem>
#include <fstream>
#include <string>

namespace walberla {


int main( int argc, char** argv )
{
    // simulation parameters files are generated for:
    std::vector< real_t > Res{ real_c( 10 ), real_c( 30 ), real_c( 60 ), real_c( 100 ), real_c( 133 ), real_c( 150 ) };
    std::vector< real_t > omegas{ real_c( 1.2 ), real_c( 1.4 ), real_c( 1.6 ) };
    std::vector< real_t > uMax_ls{ real_c( 0.2 ), real_c( 0.15 ), real_c( .1 ) };

    // file stream variable to write to files:
    std::fstream file;
    std::string filePath;

    // Strings used repeatedly in this application:
    std::string cmlString( "CMakeLists.txt" );
    std::string confString( "config.prm" );

    // IF run from build/ this command should return path to build directory:
//    std::filesystem::path buildDirPath( std::filesystem::current_path() );

    // Path to directory containing BreuerGuo benchmark files:
//    std::filesystem::path breuerGuoPath( buildDirPath /= std::filesystem::path( ".." ) /= std::filesystem::path( "examples" ) /= std::filesystem::path( "obstacle" ) /= std::filesystem::path( "BreuerGuo" ) );
//    std::filesystem::path paramStudyPath( breuerGuoPath /= std::filesystem::path( "parameterStudy" ) );

    std::string breuerGuoPath = "../examples/obstacle/BreuerGuo";
    std::string paramStudyPath = breuerGuoPath + "/parameterStudy";

    // add parameterStudy directory and add to CMakeLists:
    std::filesystem::create_directory( paramStudyPath );
    filePath = breuerGuoPath + "/" + cmlString;
    file.open( filePath, std::ios::out | std::ios::app );
    file << "\nadd_subdirectory( parameterStudy )\n";
    file.close();

    for( real_t Re : Res )
    {
        std::string reStr( "Re" + std::to_string( uint_c( Re ) ) );
        std::string rePath( paramStudyPath + "/" + reStr );
        std::filesystem::create_directory( rePath );
        
        filePath = paramStudyPath + "/" + cmlString;
        file.open( filePath, std::ios::out | std::ios::app );
        file << "add_subdirectory( " << reStr << " )\n";
        file.close();

        for( real_t omega : omegas )
        {
            std::string omegaStr( "omega" + std::to_string( uint_c( 10 * omega ) ) );
            std::string omegaPath( rePath + "/" + omegaStr );
            std::filesystem::create_directory( omegaPath );

            filePath = rePath + "/" + cmlString;
            file.open( filePath, std::ios::out | std::ios::app );
            file << "add_subdirectory( " << omegaStr << " )\n";
            file.close();

            for( real_t uMax_l : uMax_ls )
            {
                std::string uStr( "uMax" + std::to_string( uint_c( 100 * uMax_l ) ) );
                std::string uPath( omegaPath + "/" + uStr );
                std::filesystem::create_directory( uPath );

                filePath = omegaPath + "/" + cmlString;
                file.open( filePath, std::ios::out | std::ios::app );
                file << "add_subdirectory( " << uStr << " )\n";
                file.close();

                filePath = uPath + "/" + cmlString;
                file.open( filePath, std::ios::out | std::ios::app );
                file << "waLBerla_link_files_to_builddir( *.prm )\n";
                file.close();

                std::string confFilePath( uPath + "/" +  "config.prm" );
                std::string resultString( "BreuerGuo" + reStr + omegaStr + uStr );
                std::string resultPath( "/results/" +  resultString );

                lbm::generateParameterFile( confFilePath, resultPath, Re, omega, uMax_l );
            }
        }
    }
    return EXIT_SUCCESS;
}

} // namespace walberla

int main( int argc, char** argv ) { return walberla::main( argc, argv ); }