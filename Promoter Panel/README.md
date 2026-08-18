# Induced endosymbiosis between Rhizopus microsporus and Ralstonia pickettii indicates a shift from antagonism to commensalism 
### Accompanying analyses scripts

## Citation
This repository contains a collection of Python scripts for generating figure panels for the following publication: 

```
Gassler, T., et al.(2025). Induced endosymbiosis between Rhizopus microsporus and Ralstonia pickettii indicates a shift from antagonism to commensalism.
```

## Requirements

All required packages are specified in the `environment.yml` file. The scripts were developed and tested with Python 3.13.


## Setup

1. Clone this repository:
   ```bash
   git clone https://github.com/MicrobiologyETHZ/tgassler_rhizopus.git
   cd tgassler_rhizopus
   ```

2. Create and activate the conda environment:
   ```bash
   conda env create -f environment.yml
   conda activate rhizopus-env
   ```
   
   Alternatively, you can use mamba for faster dependency resolution:
   ```bash
   mamba env create -f environment.yml
   mamba activate rhizopus-env
   ```

## Usage

You can run individual scripts to generate specific figure panels:

```bash
cd code/Figure_1

python Figure_1_i_j.py

```

## License

This project is licensed under the GNU GENERAL PUBLIC License - see the [LICENSE](LICENSE) file for details.
