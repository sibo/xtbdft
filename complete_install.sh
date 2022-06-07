#!/bin/bash

#adding XTBDFT to path
echo -e "export PATH=\$PATH:$PWD/bin/ \n" >> ~/.bashrc

#install XTB v6.3.2 ("oldkernel" link may need to be uncommented, if running an older Linux version)
#newer versions of XTB may or may not work with XTBDFT v1.0
echo "Installing XTB v6.3.2..."
link=https://github.com/grimme-lab/xtb/releases/download/v6.3.2/xtb-200702.tar.xz
#link=https://github.com/grimme-lab/xtb/releases/download/v6.3.2/xtb-200709-old-kernel.tar.xz
curl -LJ0 $link > xtb.tar.xz
tar -xvf xtb.tar.xz
chmod +x xtb-6.3.2/bin/xtb
#need to edit config_env.bash file???vi
echo -e "#XTB bashrc entries below: \n export PATH=\$PATH:$PWD/xtb-6.3.2/bin/ \n ulimit -s unlimited \n OMP_STACKSIZE=4G \n source $PWD/xtb-6.3.2/share/xtb/config_env.bash \n" >> ~/.bashrc


#install CREST v2.11 ("oldkernel" link may need to be uncommented, if running an older Linux version)
#newer versions of CREST may or may not work with XTBDFT v1.0
echo "Installing CREST v2.11"
link=https://github.com/grimme-lab/crest/releases/download/v2.11/crest.tgz
#link=https://github.com/grimme-lab/crest/releases/download/v2.11/crest-oldkernel.tgz
curl -LJ0 $link > crest.tgz
tar -xvzf crest.tgz --directory xtb-6.3.2/bin/
chmod +x xtb-6.3.2/bin/crest

#install NWChem
echo "Installing NWChem..."
sudo apt-get install nwchem

#install GoodVibes 3.0.1 committ 54d0750 
#newer versions of GoodVibes have dropped Orca and NWChem support, and will NOT work with XTBDFT
echo "Installing GoodVibes commit 54d0750..."
git clone https://github.com/patonlab/GoodVibes.git
cd GoodVibes
git checkout -b goodVibes_54d0750 54d0750
echo -e "#Goodvibes bashrc entries below: \n export PYTHONPATH=\$PATH:$PWD/goodvibes/" >> ~/.bashrc 

#clean up installation
cd ..
rm xtb.tar.xz crest.tgz
echo "Sourcing .bashrc file..."
source ~/.bashrc
echo "Installation complete"
