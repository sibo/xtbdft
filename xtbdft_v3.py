#!/usr/bin/python
#
#running in headless mode: nohup python3 ~/path/to/this/file/xtb-nwchem.py geom.xyz [-chrg int] [-uhf int] [-xc str,str,str,str] [-bs str,str,str,str] > autoConf.out &
#
import os, sys, subprocess, time, math
from datetime import datetime
goodvibesPy = "~/_programs/goodvibes/goodvibes-3.0.1/goodvibes/GoodVibes.py"
xtbPath="~/_programs/xtb" # xtb and crest binaries must be located in xtbPath/bin/
default_cutoff = 2.0
default_xc = "b3lyp,,b3lyp,b3lyp"
default_bs = "def2-sv(p),,def2-svp,def2-tzvp"
default_chrg = 0
default_uhf = 0

### edit below to match your computing environment
NP=24
msubHeader="""#!/bin/bash
#MSUB -l nodes=1:ppn={0}
#MSUB -l walltime=99:00:00:00
#MSUB -l mem=100gb
""".format(NP)

### do not change below code unless you know what you're doing!
realpath = os.path.realpath(__file__)
calcName=os.path.basename(os.getcwd())
startPath=os.getcwd()

def run_crest(file,chrg,uhf):
    #if(os.path.isfile(os.getcwd()+"/crest.out")):
    #    go=input("crest.out already exists. Are you sure you want to proceed? (y/n)")
    #    if(go != 'y' and go != 'Y'):
    #        print("Okay, aborting")
    #        exit(1)
    script=open("submit_crest_{0}.sh".format(calcName), "w")
    script.write(msubHeader+"""
export OMP_STACKSIZE=1G
#source {0}/share/xtb/config_env.bash
ulimit -s unlimited
cd {1}
{2}/bin/crest {3} --chrg {4} --uhf {5} -T {6} > crest.out""".format(xtbPath,os.getcwd(),xtbPath,file,chrg,uhf,NP))
    script.close()
    #calcID = subprocess.run(['msub','submit_crest_{0}.sh'.format(calcName)], stdout=subprocess.PIPE).stdout.decode('utf-8')
    calcID = subprocess.check_output(['msub','submit_crest_{0}.sh'.format(calcName)]).decode('utf-8').replace("\n","")
    os.system("echo \"calcID is {0}\" | tee {0}.calcID".format(calcID))
    return calcID

def run_nwchem(chrg,uhf,calcType,xc,bs,cutoff,**kwargs):
    calcName2 = calcName +"_"+calcType
    if calcType == "crude":
        xc1,xc2 = xc[0],xc[1]
        bs1,bs2 = bs[0],bs[1]
    elif calcType == "refine" or calcType=="autoTS":
        xc1,xc2 = xc[2],xc[3]
        bs1,bs2 = bs[2],bs[3]     
        atomNo1 = kwargs.get('atom1')
        atomNo2 = kwargs.get('atom2')
        if calcType == "refine":
            xyz="../minimum_lowest.xyz"
            optType=""
        else:
            xyz="../TSguess.xyz"
            optType="saddle"
    input=open("{0}.nw".format(calcName2),"w")
    input.write("""memory heap 200 mb stack 1000 mb global 2800 mb
start calc
title "{3}-D3/{4}: {0}"
echo
charge {1}
dft
    mult {2}
    xc {3}
    disp vdw 3
    iterations 200
    print low
    direct
end
basis bs1 spherical
  * library {4}
end
basis "cd basis" spherical
  * library "Weigend Coulomb Fitting"
end
""".format(calcType,chrg,uhf+1,xc1,bs1))
    if bs2 != "":
        input.write("""basis bs2 spherical
  * library {0}
end
""".format(bs2))
    if calcType == "crude":
        input.write("""dft
    grid medium
    convergence energy 1d-6
end
driver
    gmax 0.002 ; grms 0.0003 ; xrms 1 ; xmax 1
    maxiter 400         
end
set "ao basis" bs1
""")
        numConfs=readConfs(cutoff)
        os.system("echo There are {0} conformers with deltaE_xtb < {1} kcal/mol".format(numConfs,cutoff))
        for i in range(1,numConfs+1):
            input.write("""###conf {0}
geometry units angstroms noautosym
    load frame {0} ../crest_conformers.xyz
end
dft
    vectors input atomic output bs1.mos
end
task dft optimize
task shell "echo conf {0} geometry opt complete"
""".format(i))
            if(xc2 != ""):
                input.write("""### single point at {0}/{1}//{2}/{3]}
set "ao basis" bs2
dft
  xc {0}
  vectors input atomic output bs2.mos
end
""".format(xc2,bs2,xc1,bs1))
    elif calcType == "refine" or calcType == "autoTS":
        input.write("""dft
    grid fine
    convergence energy 1d-8
end
driver #tightopt criteria from Orca
    gmax 0.0001 ; grms 0.00003 ; xrms 0.0006 ; xmax 0.001
    maxiter 400  
    sadstp 0.03
end
geometry units angstroms noautosym
    load {0}
end
set "ao basis" bs1
""".format(xyz))
        if calcType == "autoTS":
            os.system("mkdir -p opt")
            input.write("""geometry adjust #fix reaction coordinate (bond)
  zcoord
    bond {0} {1} constant
  end
end
task shell "echo @starting constrained opt"
task dft optimize
geometry adjust noautoz #unfix reaction coordinate, switch to cartesian opt for saddle pt search
  zcoord
    bond {0} {1}
  end
end
driver 
  xyz opt/optSaddle
end
task shell "echo @starting saddle optimization" 
""".format(atomNo1,atomNo2))        
        input.write("""task dft optimize {0}
task shell "echo @starting vibrational calculation"
task dft freq numerical
task shell "echo @starting single point calculation"
set "ao basis" bs2
dft
  xc {1}
  vectors input atomic output bs2.mos
end
task dft energy
""".format(optType,xc2))    
    else:
        print("NWChem must have a calctype: crude, refine, or autoTS")
        exit(1)
    
    #submit nwchem job to msub system
    script=open("submit_nwchem_{0}.sh".format(calcName2), "w")
    script.write(msubHeader+"""
#export ARMCI_DEFAULT_SHMMAX=8192
module load mpich-x86_64
cd {0}
mpirun -np {1} nwchem_mpich {2}.nw > nwchem.out
echo "{3} (nwchem) submitted" """.format(os.getcwd(),NP,calcName2,calcName2))
    script.close()
    calcID = subprocess.check_output(['msub','submit_nwchem_{0}.sh'.format(calcName2)]).decode('utf-8').replace("\n","")
    os.system("echo \"calcID is {0}\" | tee {0}.calcID".format(calcID))
    return calcID

def readConfs(cutoff):
	f = open("../crest.energies",'r')
	lines = f.readlines()
	numConfs=0
	for line in lines:
		l = line.split()
		if float(l[1]) < cutoff:
			numConfs = numConfs + 1
	return numConfs

def track_crest(file,calcID):
    #check every 5min if crest.out exists before proceeding
    while True:
        if (os.path.isfile("crest.out")):
            print(os.getcwd()+"\crest.out exists")
            break
        time.sleep(5*60)
    #check if crest.out is done
    while True:
        if ("CREST terminated" in subprocess.check_output(['tail','-1','crest.out']).decode('utf-8')):
            print("CREST terminated normally")
            break
        time.sleep(5*60)
    
def track_nwchem(calcID):
    #check every 5min if nwchem.out exists before proceeding
    while True:
        if (os.path.isfile("nwchem.out")):
            print(os.getcwd()+"\nnwchem.out exists")
            break
        time.sleep(5*60)
    #check if nwchem.out every 30min
    while True:
        if ("Total times" in subprocess.check_output(['tail','-1','nwchem.out']).decode('utf-8')):
            print("NWChem terminated normally")
            break
        time.sleep(30*60)

# readrawdata reads NWChem output containing multiple geometries and returns an array of geometries (xyz, Angstrom), energies (Hartrees), and boolean flag denoting if this geometry marks a geometric minimum
def getrawdata(infile):
    f=open(infile,'r')
    geo=0
    energies=[]
    xyzs=[]
    structure=[]
    isMin=[]
    itsOn = 0
    findEnergy = 0
    dumpData = 0
    recordGeo = 0
    confNum = []
    confCounter = 1
    for line in f:
        if "Step" in line and line.split()[0] == "Step":
            try: 
                int(line.split()[1])
                recordGeo = 1
            except ValueError:
                pass
        if 'Geometry \"' in line and recordGeo == 1:
            geo = 1
            recordGeo = 0
        elif ' ---- ' in line and geo == 1:
            geo = 0
            itsOn = 1
        elif itsOn == 1:
            #print line
            data = line.split()
            if data != []:
                structure += [[data[1],data[3],data[4],data[5]]]
                #print structure
            elif data == []:
                itsOn = 0
                findEnergy = 1
        elif findEnergy == 1 and 'Total DFT' in line:
            energy = float(line.split()[4])
            findEnergy = 2
        elif findEnergy == 2 and 'Optimization converged' in line:
            isMin += [True]
            dumpData = 1
        elif findEnergy == 2 and 'Line search:' in line:
            isMin += [False]
            dumpData = 1
        if dumpData == 1:
            dumpData = 0
            xyzs+=[structure]
            structure=[]
            energies+=[energy]
            confNum += [confCounter]
            if isMin[-1] == True:
                confCounter += 1
    #print len(xyzs)
    #print len(energies)
    #print energies
    #print len(isMin)
    #print isMin
    return xyzs,energies,isMin,confNum

        
def pes_parse(infile):
    xyzs,energies,isMin,confNum=getrawdata(infile)
    f=open('all.xyz','w')
    g=open('energies.dat','w')
    h=open('optHist.dat','w')
    i=open('minima.xyz','w')
    j=open('minima_basis2.xyz','w')
    k=open('minimum_lowest.xyz','w')
    minIndex = 0
    hasAltBasis = False # True
    altBasis = False
    altBasisOptHist = ""
    for n in range(len(energies)):
        if energies[n] < energies[minIndex]:
            minIndex = n
    k.write(str(len(xyzs[minIndex])) + "\n Lowest energy conformer: " + str(minIndex) + " \n")
    for n in range(len(energies)):
        energykcalmol = str((energies[n]-energies[minIndex])*627.509)
        if energykcalmol == "0.0":
            energykcalmol = "0.0000000000"
        energyString = str(confNum[n]) + ' \t' + energykcalmol +'\t rel kcal/mol \t '+str(energies[n])+' \t Hartrees'
        g.write(str(n) + "\t" + energyString + "\t" + str(isMin[n]) + "\n")
        f.write(str(len(xyzs[n])) + "\n" + energyString + "\n")
        if isMin[n] == True:
            if altBasis == True:
                j.write(str(len(xyzs[n])) + "\n" + energyString + "\n")
            elif altBasis == False:
                i.write(str(len(xyzs[n])) + "\n" + energyString + "\n")
        for atom in xyzs[n]:
            line=" ".join(atom) + "\n"
            f.write(line)
            if(isMin[n] == True and altBasis == True):
                j.write(line)
            elif isMin[n] == True and altBasis == False:
                i.write(line)
            if(n == minIndex):
                k.write(line)
        if altBasis == False:
            h.write(energykcalmol + '\t')
            if isMin[n] == True:
                h.write('\n')
                if hasAltBasis == True:
                    altBasis = True
        else:
            altBasisOptHist += (energykcalmol + '\t')
            if isMin[n] == True:
                altBasisOptHist += '\n'
                altBasis = False    
    h.write("\n\n" + altBasisOptHist)
    g.close()
    f.close()
    h.close()
    i.close()
    j.close()
    k.close()
    if hasAltBasis == False:
        os.remove("minima_basis2.xyz")

def clean_nwchem():
    os.system("rm -vf *.2ceri* *.b *.b^-1 *.c *.calcID *.cdfit *.db *.drv.hess *.gridpts.* *.p *.zmat submit*.sh*")

def goodvibes(file):
    os.system("echo starting goodvibes...")
    os.system("python3 {0} --invert -100 {1}".format(goodvibesPy,file))

def parseArgs():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("-chrg", action="store",dest="chrg",default=default_chrg,type=int)
    parser.add_argument("-uhf", action="store",dest="uhf",default=default_uhf,type=int)
    parser.add_argument("-xc",action="store",dest="xc",default=default_xc,type=str)
    parser.add_argument("-bs",action="store",dest="bs",default=default_bs,type=str)
    parser.add_argument("-cutoff",action="store",dest="cutoff",default=default_cutoff,type=float)
    parser.add_argument("-mode", action="store",dest="mode",default="autoConf",type=str,nargs='+')
    parser.add_argument("filename",action="store")
    results=parser.parse_args()
    if (not os.path.isfile(results.filename)):
        print("Error: {0} is not a file!".format(results.filename))
        sys.exit(1)
    elif results.filename[-4:] != ".xyz":
        print ("Error: Cartesian coord file must be .xyz")
        sys.exit(1)
    xc = results.xc.split(",")
    bs = results.bs.split(",")
    if (len(xc) != 4 or len(bs) != 4):
        print ("Error: -bs and -xc must be comma-delimited strings of exactly four items. Omitting them will default to them to: \n-xc b3lyp,,b3lyp,b3lyp \n-bs def2-sv(p),,def2-svp,def2-tzvp")
        sys.exit(1)
    mode = results.mode[0]
    if mode == "autoConf":
        pass
    elif mode == "autoTS" and results.mode[1].isdigit() and results.mode[2].isdigit():
        try:
            float(results.mode[3])
        except:
            print("could not convert {0} to a float".format(results.mode[3])) 
    else:
        print("Error: valid calls are '-mode autoConf' or '-mode autoTS atom1 atom2 finalDistance'")
        sys.exit(1)
    return results.filename,results.chrg,results.uhf,xc,bs,results.cutoff,mode,results.mode[1:]
    
def checkEnv():
    #check that required shell scripts are executable
    #subprocess.check_output(['chmod','+x',realpath+"/calcInp_svP.sh"])
    #subprocess.check_output(['chmod','+x',realpath+"/calcInp_svp_final.sh"])
    #subprocess.check_output(['chmod','+x',realpath+"/pes_parse.py"])
    pass

def autoConf(file,chrg,uhf,xc,bs,cutoff):
    if (os.path.isfile("crest.out")):
        if ("CREST terminated" in subprocess.check_output(['tail','-1','crest.out']).decode('utf-8')):
            print("CREST already complete, skipping to NWChem Crude Refinement")
    #run crest on input .xyz file
    else:
        os.system("echo {0}: running CREST".format(datetime.now()))
        calcID_crest = run_crest(file,chrg,uhf)
        os.system("echo {0}: tracking CREST".format(datetime.now()))
        track_crest(file, calcID_crest)
    
    #run crude NWChem re-optimization
    os.system("mkdir nwchem")
    os.chdir("nwchem")
    os.system("echo {0}: running NWChem Crude".format(datetime.now()))
    calcID_nwchem = run_nwchem(chrg,uhf,"crude",xc,bs,cutoff)
    os.system("echo {0}: tracking NWChem Crude".format(datetime.now()))
    track_nwchem(calcID_nwchem)
    clean_nwchem()
    
    #determine the lowest energy conformer and refine using NWChem
    os.system("echo {0}: running pes_parse subroutine".format(datetime.now()))
    pes_parse("nwchem.out")
    os.system("mkdir refine")
    os.chdir("refine")
    os.system("echo {0}: running NWChem refine".format(datetime.now()))
    calcID_nwchem = run_nwchem(chrg,uhf,"refine",xc,bs,cutoff)
    os.system("echo {0}: tracking NWChem refine".format(datetime.now()))
    track_nwchem(calcID_nwchem)
    clean_nwchem()
    
    #optional: script goodvibes correction
    os.system("echo {0}: running GoodVibes.py".format(datetime.now()))
    goodvibes("nwchem.out")    

def getBondDistance(xyzFile,atomNo1,atomNo2):
    f=open(xyzFile,"r")
    f1 = f.readlines()
    atomCoords=[]
    for l in f1[2:]:
        print(l)
        coordinates = l.split()
        atomCoords.append([float(coordinates[1]),float(coordinates[2]),float(coordinates[3])])
    atom1 = ["Atom1"] + atomCoords[int(atomNo1)-1]
    atom2 = ["Atom2"] + atomCoords[int(atomNo2)-1]
    return distance(atom1,atom2)             
    #deltaX = atomCoords[atomNo1-1][0] - atomCoords[atomNo2-1][0]
    #deltaY = atomCoords[atomNo1-1][1] - atomCoords[atomNo2-1][1]
    #deltaZ = atomCoords[atomNo1-1][2] - atomCoords[atomNo2-1][2]
    #print deltaX
    #return math.pow(math.pow(deltaX,2)+math.pow(deltaY,2)+math.pow(deltaZ,2),0.5)

def distance(atom1, atom2):
    #atom1 and atom2 are arrays of strings
    #atomX[0] = chemical symbol
    #atomX[1] = x coord
    #atomX[2] = y coord
    #atomX[3] = z coord
    deltaX = float(atom1[1]) - float(atom2[1])
    deltaY = float(atom1[2]) - float(atom2[2])
    deltaZ = float(atom1[3]) - float(atom2[3])
    return math.pow(math.pow(deltaX,2)+math.pow(deltaY,2)+math.pow(deltaZ,2),0.5)

def autoTS(file,chrg,uhf,xc,bs,atomNo1,atomNo2,finalScanDistance,direction):
    #import matplotlib.pyplot as plt
    
    #set up XTB scan
    os.system("echo {0}: setting up XTB xconstrains file".format(datetime.now()))
    steps=100
    startScanDistance=getBondDistance(file,atomNo1,atomNo2)
    dir = "scan_" + str(atomNo1) + "_" + str(atomNo2) + "_" + direction
    os.system("mkdir -p " + dir)
    os.system("cp " + file + " " + dir + "/")
    f=open(dir + "/xcontrol","w+")
    command = "$constrain \n  force constant = 1 \n  distance: " + str(atomNo1) + "," + str(atomNo2) + "," + str(startScanDistance) + "\n$scan\n  1: " + str(startScanDistance) + "," + str(finalScanDistance) + "," + str(steps) + " \n$end\n"
    f.write(command)
    f.close()
    os.chdir(dir)
    os.system("echo {0}: running XTB".format(datetime.now()))
    os.system("xtb {0} --chrg {1} --uhf {2} --opt veryfine --input xcontrol | tee opt.out".format(os.path.basename(file),chrg,uhf))
    
    #parse XTB scan results
    f=open("xtbscan.log","r")
    f1=f.readlines()
    totalAtoms=int(f1[0])
    print("total atoms = " + str(totalAtoms))
    energies=[]
    structures=[]
    bondLengths=[]
    energy0 = float(f1[1].split()[1])*627.509
    for i in range(0,steps):
        energy = float(f1[i*(totalAtoms+2)+1].split()[1])*627.509
        energies.append(energy-energy0)
        print(str(energy-energy0))
        coords = []
        for line in f1[i*(totalAtoms+2)+2:(i+1)*(totalAtoms+2)]:
            coords.append(line.split())
        structures.append(coords)
        bondLengths.append(distance(coords[int(atomNo1)-1],coords[int(atomNo2)-1]))
    f=open("scan.csv","w+")
    for i in range(0,len(energies)):
        f.write(str(bondLengths[i]) + "," + str(energies[i]) + "\n")
    f.close()
    #plt.scatter(bondLengths,energies)
    #plt.xlabel(r"Bond Distance ($\AA$)")
    #plt.ylabel("relative E (kcal/mol)")
    #plt.savefig("scanPES.png")
    maxEnergy = energies[0]
    maxStruct = []
    for i in range(1,len(energies)):
        if energies[i] > maxEnergy:
            maxEnergy = energies[i]
            maxStruct = structures[i]
    f=open("TSguess.xyz","w+")
    f.write(str(totalAtoms) + "\n TS guess, E(xtb) = " + str(maxEnergy) + " \n")
    for coord in maxStruct:
        f.write(" ".join(coord) + "\n")
    f.close()
    
    #run NWChem constrained optimization, TS opt, frequency calc, and single-point energy evaluation
    os.system("mkdir nwchem")
    os.chdir("nwchem")
    os.system("echo {0}: running NWChem Constrained Opt, then TS opt, and then Freq at {1}/{2} and Single-Point evaluation at {3}/{4}".format(datetime.now(),xc[2],bs[2],xc[3],bs[3]))
    calcID_nwchem = run_nwchem(chrg,uhf,"autoTS",xc,bs,cutoff,atom1=atomNo1,atom2=atomNo2)
    os.system("echo {0}: tracking NWChem...".format(datetime.now()))
    track_nwchem(calcID_nwchem)
    clean_nwchem()
    os.system("echo {0}: finished NWChem!".format(datetime.now()))

    #optional: script goodvibes correction
    os.system("echo {0}: running GoodVibes.py".format(datetime.now()))
    goodvibes("nwchem.out")   


if __name__ == "__main__":
    print("Reminder on how to run xtb-nwchem.py in headless mode:")
    print("nohup python3 ~/path/to/this/file/xtbdft.py geom.xyz [-chrg int] [-uhf int] > yourOutputFile.out &")
    pid = os.getpid()
    print("Process id = " + str(pid))
    file,chrg,uhf,xc,bs,cutoff,mode,params=parseArgs()
    checkEnv()
    if (mode == "autoConf"):
        autoConf(file,chrg,uhf,xc,bs,cutoff)
    elif (mode == "autoTS"):
        autoTS(file,chrg,uhf,xc,bs,params[0],params[1],params[2],"forward")
    