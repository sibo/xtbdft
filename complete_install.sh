#!/bin/bash

#install XTB
git clone https://github.com/grimme-lab/xtb.git
git checkout -b xtb_6.3.2 452493c

#install CREST
git clone https://github.com/grimme-lab/crest.git
git checkout -b crest_2.11 348d53e

#install NWChem
sudo apt-get install nwchem

#install GoodVibes 3.0.1 committ 54d0750 (newer versions of GoodVibes have dropped Orca and NWChem support)
git clone https://github.com/patonlab/GoodVibes.git
git checkout -b goodVibes_54d0750 54d0750
