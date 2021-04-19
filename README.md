# GAMMA

GAMMA is a code for modeling relativistic hydrodynamics and non-thermal emission on a moving.

## Installation

Clone this repository in a directory of your choice.

```bash
git clone https://github.com/eliotayache/GAMMA/ dirname
```

## Usage

Output is stored in ./results/Last

Before a new computation, create this directory.

```bash
mkdir -p results/Last
```

Use mpirun to launch computation.

```bash
mpirun -n N_nodes ./bin/GAMMA -w  # to overwrite files in results/Last
mpirun -n N_nodes ./bin/GAMMA -r  # to resume from the last file in results/Last
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)
