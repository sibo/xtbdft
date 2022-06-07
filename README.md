# xtbdft

This is the official repository of the `xtbdft` program, a wrapper script for multi-level molecular modelling powered by CREST/GFN2-XTB and NWChem (DFT)

## Pre-Requisites
- Python3.x
- [XTB](https://github.com/grimme-lab/xtb/releases) (version 6.3.2)
- [CREST](https://github.com/grimme-lab/crest/releases) (version 2.11)
- [NWChem](https://github.com/nwchemgit/nwchem/releases) (version 6.8)
- [Goodvibes](https://github.com/patonlab/GoodVibes/tree/54d0750b0ba7aa9121c284519271a9a0bd0764a9) (version 3.0.1, commit 54d0750)

XTBDFT has not been tested with, but may happen to work with, newer versions of XTB, CREST, and NWChem.

## Standard Installation
[latest release page](https://github.com/sibo/xtbdft/releases/latest)
1. Clone this repository in your desired directory and navigate into it:
```bash
git clone https://https://github.com/sibo/xtbdft.git
cd xtbdft
```
2a. If performing a complete clean installation, just run complete_install.py
```bash
./complete_install.sh
```
2b. Or if there are pre-existing installations of xtb, crest, nwchem, and/or goodvibes, comment out the relevant section(s) from .complete_install.sh and then execute it.
3. Modify the first 25 lines of bin/xtbdft.py to fit your computing environment, if necessary.

## Usage
```bash
nohup xtbdft.py guess.xyz [-chrg int] [-uhf int] [-xc str,str,str,str] [-bs str,str,str,str] [-mode autoConf|autoTS] [-other=["skipCrest"|crestParameters] &
```
For determining the lowest energy conformation of a neutral, singlet-spin molecule (as in reference 1 below, the input can be simplified to:
```bash
nohup xtbdft.py guess.xyz &
```
For determining the lowest energy conformation of a cationic, quartet-spin molecule, containing a PhCl bond that undergoes undesired oxidative addition under normal CREST parameters:
```bash
nohup xtbdft.py guess.xyz -chrg 1 -uhf 3 -other="-cbonds 0.1" &
```
If instead you'd like to scan for a transition state of a monocationic, doublet species, in which the distance between atoms X and Y is adjusted to Z Angstroms over 100 steps, the input is:
```bash
nohup xtbdft.py guess.xyz -chrg 1 -uhf 1 -mode autoTS X Y Z &
```
XTBDFT is currently configured to run on a compute cluster with MSUB scheduler, however, it can be run locally for debugging purposes:
```bash
xtbdft.py guess.xyz -local true
```


## Citations

1. Lin, S.; Fromer, J. C.; Ghosh, Y.; Hanna, B.; Elanany, M.; Xu, Wei "Computer-assisted catalyst development via automated modelling of conformationally complex molecules: application to diphosphinoamine ligands" <i>Sci. Rep.</i> <b>2021</b>, <i>11</i>, 4534. DOI: <a href="https://doi.org/10.1038/s41598-021-82816-x">10.1038/s41598-021-82816-x</a>

## License

`xtbdft` is free software: you can redistribute it and/or modify it under
the terms of the MIT License.


