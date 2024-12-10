#!/bin/bash

# Create python environment and install dependencies
echo "Setting up python environment"

python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt

echo "python environment ready"

# Setup waLBerla installation
echo "Cloning waLBerla repository"
git clone https://i10git.cs.fau.de/walberla/walberla.git

echo "Setting waLBerla build options"

WALBERLA_ARGS=( "-DBUILD_TESTING=OFF" \
                "-DWALBERLA_BUILD_BENCHMARKS=OFF" \
                "-DWALBERLA_BUILD_DOC=OFF" \
                "-DWALBERLA_BUILD_SHOWCASES=OFF" \
                "-DWALBERLA_BUILD_TESTS=OFF" \
                "-DWALBERLA_BUILD_TOOLS=OFF" \
                "-DWALBERLA_BUILD_TUTORIALS=OFF" \
                "-DWALBERLA_BUILD_WITH_CODEGEN=ON" \
                "-DWALBERLA_BUILD_WITH_CUDA=OFF" \
                "-DWALBERLA_BUILD_WITH_FASTMATH=OFF" \
                "-DWALBERLA_BUILD_WITH_FFTW=OFF" \
                "-DWALBERLA_BUILD_WITH_GCOV=OFF" \
                "-DWALBERLA_BUILD_WITH_GPROF=OFF" \
                "-DWALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT=OFF" \
                "-DWALBERLA_BUILD_WITH_HIP=OFF" \
                "-DWALBERLA_BUILD_WITH_LIKWID_MARKERS=OFF" \
                "-DWALBERLA_BUILD_WITH_LTO=OFF" \
                "-DWALBERLA_BUILD_WITH_METIS=OFF" \
                "-DWALBERLA_BUILD_WITH_MPI=OFF" \
                "-DWALBERLA_BUILD_WITH_OPENMESH=OFF" \
                "-DWALBERLA_BUILD_WITH_OPENMP=OFF" \
                "-DWALBERLA_BUILD_WITH_PARMETIS=OFF" \
                "-DWALBERLA_BUILD_WITH_PYTHON=OFF" \
                "-DWALBERLA_DEPS_ERROR=OFF" \
                "-DWALBERLA_DOUBLE_ACCURACY=ON" \
                "-DWALBERLA_ENABLE_GUI=OFF" \
                "-DWALBERLA_GUI_USE_3D_BLOCKVIEW=OFF" \
                "-DWARNING_ERROR=OFF" \
                "-DWARNING_PEDANTIC=ON" \
                "-DWALBERLA_DIR=../walberla/" )

# Setup build directory
mkdir build
cd build

cmake "${WALBERLA_ARGS[@]}" ..
make -j16

cd ..