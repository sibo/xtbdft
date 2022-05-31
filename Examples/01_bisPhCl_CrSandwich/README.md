Example 1: [Cr(η<sup>6</sup>-PhCl)<sub>2</sub>]<sup>+</sup>
========================================
by Sibo Lin

### Command line:

Enter this example code in the command line, in this directory, to explore the conformational space of [Cr(η<sup>6</sup>-PhCl)<sub>2</sub>]<sup>+</sup>:

`nohup xtbdft.py start.xyz -chrg 1 -uhf 1 -other="-cbonds 0.02" &`

Or, to perform conformer ensemble DFT re-optimization at the B3LYP+D3/def2-TZVP level of theory, enter:

`nohup xtbdft.py start.xyz -chrg 1 -uhf 1 -xc "b3lyp,,b3lyp,b3lyp" -bs "def2-tzvp,,def2-tzvp,def2-tzvp" -other="-cbonds 0.02" &`
