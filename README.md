# LBMCombinedCollisionOperators
My TUM CSE Master's thesis' code project

## Installation

### Automated procedure

After cloning this repository, there are two options: on a Linux machine, the script `install.sh` can be used to build everything, given `CMAKE` and a suitable `C++`compiler are installed. To run the installation (and afterwards simulations), make shell scripts executable via

```
chmod +x *.sh
```

Aferwards, the project (with dependencies) can be built via 

```
./install.sh
```

### Custom procedure

Refer to waLBerla's [manual page](https://walberla.net/doxygen/setup_instructions.html) for custom installation. Dependencies such as `OpenMesh`, `MPI`, or `OpenMP` are not necessary to run this project. waLBerla's compilation process can be accellerated by not building tutorials, showcases, and benchmarks. This can be done by adjusting build options in waLBerla's `CMakeLists.txt`. For further detail, please refer to waLBerla's documentation. To install all python dependencies, it is recommended to set up and activate a virtual python environment first, e.g. via 

```
python3 -m venv .venv
source .venv/bin/activate
```

In case your python environment is named differently, you should add it to the [.gitignore](./.gitignore). All python dependencies (with or without a virtual environment) can be installed via

```
pip install -r requirements.txt
```

With all requirements insalled, create a build directory and step into it via

```
mkdir build
cd build
```

From he build directory, call

```
cmake -DWALBERLA_DIR=<path to waLBerla sources> ..
make -j16
```

Note that this procedure takes quite a while when executed for the first time. To run simulations using the scripts provided, step back into the projects root directory via

```
cd ..
```

## Running the benchmark

In case you did not run the automated installation routine, give execution rights to the project's shell scripts now via

```
chmod +x *.sh
```

Otherwise, you have to call all simulations manually, which is possible but cumbersome and typo-prone.

### Coarse grid validation (recommended)
To run the simulation on the two grids that are similar to the ones used in [reference literature](./referenceLiterature/), just type 

```
./runCoarseGrids.sh
```

from the project's root directory. Plots will be generated and saved [here](./figures/).

### Parameter study (if you are very patient)
To run the entire parameter study, type

```
./runParamStudy.sh
```

from the project's root directory. Plots are saved [here](./figures/)

### Alternative way to run parameter study

instead of running the entire parameter study in one process, you can split it up and run the separate `runRe<Reynolds Number>.sh` scripts. To do so for e.g. Re = 100, step into the build directory and call the script from there (all scripts are linked to the build directory) via

```
cd build
./runRe100.sh
```

To reproduce plots from the thesis, run the respective postprocessing script (still from the build directory) afterwards via

```
python3 parameterStudy.py
```

Again, plots will be saved [here](./figures/).