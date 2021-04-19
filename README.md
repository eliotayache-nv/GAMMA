# GAMMA

GAMMA is a code for modeling relativistic hydrodynamics and non-thermal emission on a moving mesh.

Parts of the code are currently still under development. A full documentation will be made available shortly.
For now we provide information on how to compile and run the test cases included with the code.

## Requirements

- gnu scientific library
- OpenMP
- MPI

## Installation

Clone this repository in a directory of your choice. Move to the directory and build the code.

```bash
git clone https://github.com/eliotayache/GAMMA/ dirname
cd dirname
make -B
```

## Usage

### Choice of module to use

The choice of geometry, time-integration, solver, dimensions are specified in ./Makefile
```Makefile
INITIAL    = Tests/BM1D      # Initial setup: .cpp file in src/Initial (see test examples)
TIMESTEP   = rk3             # euler / rk3
GEOMETRY   = spherical1D     # cartesian / spherical / spherical1D
HYDRO      = rel_sph         # rel_cart / rel_sph
RADIATION  = radiation_sph   # only one option for now
SOLVER     = hllc            # only one option for now
DIMENSIONS = 1d              # 1d / 2d
IO         = text1d          # text1d / text2d
```

### environment variables

A range of self-explanatory environement variables are specified in src/environment.h and should be set before running the code.

### initial setup

New initial setups can be created as new .cpp files in src/Initial. 
These files should follow the same architecture as the example tests in src/Initial/Tests.

In `initialValues()` the fluid state should be specified in terms of primitive variables.


### running a simulation

Output is stored in ./results/Last. Before you can run calculations you need to create this directory.

```bash
mkdir -p results/Last
```

Use mpirun to launch computation.

```bash
mpirun -n N_nodes ./bin/GAMMA -w  # to overwrite files in results/Last
mpirun -n N_nodes ./bin/GAMMA -r  # to resume from the last file in results/Last
```


## License
[MIT](https://choosealicense.com/licenses/mit/)

## Contact
Eliot Ayache: e.h.r.ayache[at]bath.ac.uk
