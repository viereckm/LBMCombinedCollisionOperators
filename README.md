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

Refer to waLBerla's [manual page](https://walberla.net/doxygen/setup_instructions.html) for custom installation. Dependencies such as `OpenMesh`, `MPI`, or `OpenMP` are not necessary to run this project. For codegeneration, install `lbmpy` and dependencies via 

```
pip install lbmpy
```

For installation options, refer to the `lbmpy` [repository](https://i10git.cs.fau.de/pycodegen/lbmpy). For postprocessing, [`matplotlib`](https://matplotlib.org) is required. It can be installed via

```
pip install matplotlib
```

## Running the benchmark

### Coarse grid validation (reccommended)
To run the simulation on the two grids that are similar to the ones used in [reference literature](./referenceLiterature/), just type 

```
./runCoarseGrids.sh
```

from the project's base directory. Plots will be generated and saved [here](./figures/).

### Parameter study (if you are very patient)
To run the entire parameter study (NOT reccommended due to its long duration), type

```
./runParamStudy.sh
```

from the project's base directory. Plots are saved [here](./figures/)
