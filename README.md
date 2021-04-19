# GAMMA

GAMMA is a code for modeling relativistic hydrodynamics and non-thermal emission on a moving.

## Installation

Clone this repository in a directory of your choice.

```bash
git clone https://github.com/eliotayache/GAMMA/ dirname
```

## Usage

Output is stored in ./results/Last. you can run calculations you need to create this directory.

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
