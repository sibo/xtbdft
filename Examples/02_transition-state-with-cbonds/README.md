Example 2: Transition State Conformational Search with Constrained Bonds
========================================
by Sibo Lin

### Command line:

Enter this example code in the command line, in this directory, to explore the conformational space of ethylene insertion into a chromacycloheptane: 

`nohup xtbdft.py start.xyz -chrg 1 -uhf 3 -mode autoTS 28 65 1.6 -other="-cbonds 1.0" > autoTS.out &`

- **start.xyz**: a structure of the reactant geometry
- **-chrg and -uhf**: this species is cationic and quartet-spin
- **-mode**: this flag toggles between autoConf (default) and autoTS modes
- **autoTS 28 65 1.6**: XTBDFT will scan from start.xyz, adjusting the distance between carbon atoms 28 and 65 toward 1.6 Å (in 100 steps)
- **-other="cbonds 1.0"**: Apply a harmonic constraining potential on all bonds of 1.0 Hartree/Bohr<sup>2</sup>. Default CREST conformational searching parameters (when omitting this flag) may result in undesired Cr-P bond dissociation. 
