# xtbdft

This is the official repository of the `xtbdft` program, a wrapper script for multi-level molecular modelling powered by CREST/GFN2-XTB and NWChem (DFT)

## Pre-requisites
- Python3.x
- [XTB](https://github.com/grimme-lab/xtb/releases) (version 6.3.2)
- [CREST](https://github.com/grimme-lab/crest/releases) (version 2.11)
- [NWChem](https://github.com/nwchemgit/nwchem/releases) (version 6.8)

XTBDFT has not been tested with, but may happen to work with, newer versions of XTB, CREST, and NWChem.

## Installation
[latest release page](https://github.com/sibo/xtbdft/releases/latest)

## Usage
```bash
nohup python3 ~/path/to/script/xtbdft.py coords.xyz > autoConf.out &
```

## Citations

Lin, S.; Fromer, J. C.; Ghosh, Y.; Hanna, B.; Elanany, M.; Xu, Wei "Computer-assisted catalyst development via automated modelling of conformationally complex molecules: application to diphosphinoamine ligands" <i>Sci. Rep.</i> <b>2021</b>, <i>11</i>, 4534. DOI: <a href="https://doi.org/10.1038/s41598-021-82816-x">10.1038/s41598-021-82816-x</a>

## License

`xtbdft` is free software: you can redistribute it and/or modify it under
the terms of the MIT License.


