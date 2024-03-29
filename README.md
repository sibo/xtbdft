# xtbdft

This is the official repository of the `xtbdft` program, a wrapper script for multi-level molecular modelling powered by CREST/GFN2-xTB, NWChem (DFT), and Goodvibes (anharmonic frequency corrections). XTBDFT and all underlying softwares are open-source.

## Pre-Requisites
- Python3.x
- [GFN-xTB](https://github.com/grimme-lab/xtb/releases) (version 6.3.2, GNU Lesser General Public License v3)
- [CREST](https://github.com/grimme-lab/crest/releases) (version 2.11, GNU Lesser General Public License v3)
- [NWChem](https://github.com/nwchemgit/nwchem/releases) (version 6.8, Educational Community License v2)
- [Goodvibes](https://github.com/patonlab/GoodVibes/tree/54d0750b0ba7aa9121c284519271a9a0bd0764a9) (version 3.0.1, commit 54d0750, MIT License)
- [Moab](https://adaptivecomputing.com/moab-hpc-suite/) or [Slurm](https://github.com/mattthias/slurm) (Mattthias version 0.4.3) workload scheduler

XTBDFT has not been tested with newer versions of XTB, CREST, NWChem, or Slurm. 

## Installation
0. The installation of a workload scheduler such as Moab or Slurm is beyond the scope of these instructions, but may be accomplished [following instructions elsewhere](http://web.archive.org/web/20210708045411/https://blog.llandsmeer.com/tech/2020/03/02/slurm-single-instance.html).
2. Clone this repository (or unzip a [static release](https://github.com/sibo/xtbdft/releases/latest)) and navigate into the new directory in a terminal:
```bash
git clone https://https://github.com/sibo/xtbdft.git
cd xtbdft
```
2a. If performing a complete clean installation, just run complete_install.py
```bash
./complete_install.sh
```
2b. Or if there are pre-existing installations of xtb, crest, nwchem, and/or goodvibes, comment out the relevant section(s) from .complete_install.sh and then execute it. This option is particularly useful if you choose to compile NWChem from source to optimize DFT performance.

3. Modify the first 25 lines of bin/xtbdft.py to fit your computing environment, if necessary.
- goodvibes.py location (default ~/xtbdft/GoodVibes/goodvibes/GoodVibes.py)
- path for xtb and crest binaries (default ~/xtbdft/xtb-6.3.2)
- nodes per calculation (default 1)
- CPUs per calculation (default 24)
- maximum wall time per calculation (default 14 days)
- memory per calculation (default 100GB)

## Usage
Note: Commonly, XTBDFT is run on a remote computing cluster via SSH. It is convenient to bracket the XTBDFT command with "nohup" and "&", to run the job in the background and allow it to keep running even after closing the SSH terminal or powering off your computer. If running XTBDFT locally, it is okay to omit "nohup" and "&".

```bash
nohup xtbdft.py guess.xyz [-chrg int] [-uhf int] [-xc str,str,str,str] [-bs str,str,str,str] [-mode autoConf|autoTS] [-other=["skipCrest"|crestParameters] &
```
For determining the lowest energy conformation of a neutral, singlet-spin molecule (as in reference 1 below, the input can be simplified to:
```bash
nohup xtbdft.py guess.xyz &
```
For determining the lowest energy conformation of a cationic, quartet-spin transition metal complex, containing a PhCl ligand that undergoes undesired oxidative addition under normal CREST parameters:
```bash
nohup xtbdft.py guess.xyz -chrg 1 -uhf 3 -other="-cbonds 0.02" &
```
For generating a guess transition state of a monocationic, doublet species, in which the distance between atoms X and Y is adjusted to Z Angstroms over 100 steps, but skipping CREST conformation searching, the input is:
```bash
nohup xtbdft.py guess.xyz -chrg 1 -uhf 1 -mode autoTS X Y Z -other="skipCrest" &
```
To scan for a transition state of a monocationic, doublet species, in which the distance between atoms X and Y is adjusted to Z Angstroms over 100 steps, and then constrain all bonds during CREST conformational searching with spring constant 1.0 Hartree/au, the input is:
```bash
nohup xtbdft.py guess.xyz -chrg 1 -uhf 1 -mode autoTS X Y Z -other="-cbonds 1.0" &
```
** New in XTBDFT V1.2 ** To do a conformer serach of a reasonable guess structure for a transition state of a monocationic, doublet species, in which atoms 1,2,[...],N are held constant, the input is:
```bash
nohup xtbdft.py guess.xyz -chrg 1 -uhf 1 -mode TSconf atomNo1 atomNo2 [...] atomNoN &
```



XTBDFT is currently configured to run on a compute cluster with MSUB scheduler, however, it can be run locally for debugging purposes:
```bash
xtbdft.py guess.xyz -local true
```


## Citations

1. Lin, S.; Fromer, J. C.; Ghosh, Y.; Hanna, B.; Elanany, M.; Xu, W. "Computer-assisted catalyst development via automated modelling of conformationally complex molecules: application to diphosphinoamine ligands" <i>Sci. Rep.</i> <b>2021</b>, <i>11</i>, 4534. DOI: <a href="https://doi.org/10.1038/s41598-021-82816-x">10.1038/s41598-021-82816-x</a>
2. Lin, S.; Elanany, M.; Khawji, M. "XTBDFT: Automated workflow for conformer searching of minima and transition states powered by extended tight binding and density functional theory" <i>SoftwareX</i> <b>2022</b>, <i>20</i>, 101242. DOI: <a href="https://www.sciencedirect.com/science/article/pii/S2352711022001601?via%3Dihub">10.1016/j.softx.2022.101242</a>

## License

`xtbdft` is free software: you can redistribute it and/or modify it under
the terms of the MIT License.


