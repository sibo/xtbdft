#!/bin/bash

#install XTB v6.3.2 ("oldkernel" link may need to be uncommented, if running an older Linux version)
#newer versions of XTB may or may not work with XTBDFT v1.0
link=https://github.com/grimme-lab/xtb/releases/download/v6.3.2/xtb-200702.tar.xz
#link=https://github.com/grimme-lab/xtb/releases/download/v6.3.2/xtb-200709-old-kernel.tar.xz
curl $link > xtb.tar.xz
tar -xvf xtb.tar.xz
chmod +x xtb/bin/xtb
#need to edit config_env.bash file???
echo -e "#XTB bashrc entries below: \n export PATH=\$PATH:$PWD/xtb/bin/ \n ulimit -s unlimited \n OMP_STACKSIZE=4G \n source $PWD/xtb/share/xtb/config_env.bash \n" >> ~/.bashrc

#install CREST v2.11 ("oldkernel" link may need to be uncommented, if running an older Linux version)
#newer versions of CREST may or may not work with XTBDFT v1.0
link=https://github.com/grimme-lab/crest/releases/download/v2.11/crest.tgz
#link=https://github.com/grimme-lab/crest/releases/download/v2.11/crest-oldkernel.tgz
curl $link > crest.tgz
tar -xvzf crest.tgz
chmod +x crest/bin/crest
echo -e "#CREST bashrc entries below: \n export PATH=\$PATH:$PWD/crest/bin/ \n" >> ~/.bashrc

#install NWChem
sudo apt-get install nwchem

#install GoodVibes 3.0.1 committ 54d0750 
#newer versions of GoodVibes have dropped Orca and NWChem support, and will NOT work with XTBDFT
git clone https://github.com/patonlab/GoodVibes.git
git checkout -b goodVibes_54d0750 54d0750
echo -e "#Goodvibes bashrc entries below: \n export PYTHONPATH=\$PATH:$PWD/goodVibes_54d0750/goodvibes/" >> ~/.bashrc 

source ~/.bashrc
