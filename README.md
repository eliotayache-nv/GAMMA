# GAMMA

GAMMA is a code for modeling relativistic hydrodynamics and non-thermal emission on a moving.
Parts of the code are currently still under development. A full documentation will be made available shortly.
For now we provide a information on how to compile and run the test cases included with the code.

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

Output is stored in dirname/results/Last. Before you can run calculations you need to create this directory.

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
