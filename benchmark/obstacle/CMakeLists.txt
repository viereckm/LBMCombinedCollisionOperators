waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_executable( NAME squareObstacle
                         FILES squareObstacle.cpp
                         DEPENDS blockforest core field lbm geometry timeloop )

waLBerla_add_executable( NAME createParameterFiles
                         FILES createParameterFiles.cpp
                         DEPENDS lbm )

add_subdirectory( BreuerGuo )
