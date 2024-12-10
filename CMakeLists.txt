CMAKE_MINIMUM_REQUIRED( VERSION 3.14 )

PROJECT( THESIS_PROJECT )

# Set C++ Standard
set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS OFF )

# Set build options for waLBerla
option( WALBERLA_BUILD_WITH_PYTHON   ON )
option( WALBERLA_BUILD_WITH_CODEGEN  ON )
# accelerate build process:
option( WALBERLA_BUILD_TESTS        OFF )
option( WALBERLA_BUILD_BENCHMARKS   OFF )
option( WALBERLA_BUILD_TUTORIALS    OFF )
option( WALBERLA_BUILD_SHOWCASES    OFF )

# third party libraries:
option( WALBERLA_BUILD_WITH_MPI      OFF )
option( WALBERLA_BUILD_WITH_FFTW     ON )

# Add third party source directories to prefix path:
set( CMAKE_PREFIX_PATH "/home/viem/linux/local/fftw-3.3.10/install; /home/viem/linux/local/pfft/install;" )

# Set walberla source directory:
set( WALBERLA_DIR "${CURRENT_DIR}/../walberla/" )

# Add local directories:
include_directories( ${THESIS_PROJECT_BINARY_DIR}/src )

# Add local directory to module path so cmake finds waLBerla via FindwaLBerla.cmake
list( APPEND CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${THESIS_PROJECT_SOURCE_DIR} "/home/viem/linux/local/OpenMesh-10.0.0/cmake/")
find_package( waLBerla )

# Add subdirectories:
add_subdirectory( benchmark )
add_subdirectory( src )
add_subdirectory( PythonCodegenSweeps )

waLBerla_link_files_to_builddir( *.sh )
waLBerla_link_files_to_builddir( *.py )
