# xtbdft

This is the official repository of the `xtbdft` program, a wrapper script for multi-level molecular modelling powered by CREST/GFN2-XTB and NWChem (DFT)

## Pre-Requisites
- Python3.x
- [XTB](https://github.com/grimme-lab/xtb/releases) (version 6.3.2)
- [CREST](https://github.com/grimme-lab/crest/releases) (version 2.11)
- [NWChem](https://github.com/nwchemgit/nwchem/releases) (version 6.8)

XTBDFT has not been tested with, but may happen to work with, newer versions of XTB, CREST, and NWChem.

## Installation
[latest release page](https://github.com/sibo/xtbdft/releases/latest)
- 

## Usage
```bash
nohup python3 ~/path/to/script/xtbdft.py guess.xyz [-chrg int] [-uhf int] [-xc str,str,str,str] [-bs str,str,str,str] [-mode autoConf|autoTS] > outputFile.out &
```
For determining the lowest energy conformation of a neutral, singlet molecule (as in reference 1 below, the input can be simplified to:
```bash
nohup python3 ~/path/to/script/xtbdft.py guess.xyz > outputFile.out &
```
If instead you'd like to scan for a transition state of a monocationic, doublet species, in which the distance between atoms X and Y is adjusted to Z Angstroms over 100 steps, the input is:
```bash
nohup python3 ~/path/to/script/xtbdft.py guess.xyz -chrg 1 -uhf 1 -mode autoTS X Y Z > outputFile.out &
```

## Citations

1. Lin, S.; Fromer, J. C.; Ghosh, Y.; Hanna, B.; Elanany, M.; Xu, Wei "Computer-assisted catalyst development via automated modelling of conformationally complex molecules: application to diphosphinoamine ligands" <i>Sci. Rep.</i> <b>2021</b>, <i>11</i>, 4534. DOI: <a href="https://doi.org/10.1038/s41598-021-82816-x">10.1038/s41598-021-82816-x</a>

## License

`xtbdft` is free software: you can redistribute it and/or modify it under
the terms of the MIT License.


