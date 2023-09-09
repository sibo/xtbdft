#!/usr/bin/env python3
#
#running in headless mode: nohup python3 ~/path/to/this/file/xtb-nwchem.py geom.xyz [-chrg int] [-uhf int] [-xc str,str,str,str] [-bs str,str,str,str] > autoConf.out &
#
import os, sys, subprocess, time, math
from datetime import datetime
goodvibesPy = "~/_programs/goodvibes/goodvibes-3.0.1/goodvibes/GoodVibes.py" 
# goodvibesPy = "~/xtbdft/GoodVibes/goodvibes/GoodVibes.py"
xtbPath="~/_programs/xtb" # xtb and crest binaries must be located in xtbPath/bin/ 
# xtbPath = "~/xtbdft/xtb-6.3.2" #xtb and crest binaries must be located in xtbPath/bin/
default_cutoff = 3.0
default_xc = "b3lyp,,b3lyp,b3lyp"
default_bs = "def2-svp,,def2-svp,def2-tzvp" #xtbdft_v4 default: "def2-sv(p),,def2-svp,def2-tzvp"
default_chrg = 0
default_uhf = 0

### edit below to match your computing environment
nodes=1
NP=24
max_walltime_days=14
mem_gb=125

if os.system("which squeue") == 0:
	msubHeader="""#!/bin/bash
	#SBATCH --nodes={1}
	#SBATCH --ntasks=1
	#SBATCH --cpus-per-task={0}
	#SBATCH --time={2}-00:00:00
	#SBATCH --mem={3}G
	""".format(NP,nodes,max_walltime_days,mem_gb)
	qsub="sbatch -c {0}".format(NP)
else:
	msubHeader="""#!/bin/bash
	#PBS -l nodes={1}:ppn={0}
	#PBS -l walltime={2}:00:00:00
	#PBS -l mem={3}gb
	""".format(NP,nodes,max_walltime_days,mem_gb)
	qsub="qsub"
	
### do not change below code unless you know what you're doing!
if not (os.path.exists(os.path.expanduser(goodvibesPy))):
    print("Error: GoodVibes.py does not exist at " + goodvibesPy + "\n Please modify line 7 of xtbdft.py to reflect where GoodVibes.py is located.\nExiting...")
    exit()
if not (os.path.exists(os.path.expanduser(xtbPath + "/bin/xtb"))):
    print("Error: xtb does not exist at " + xtbPath + "/bin/xtb \n Please modify line 8 of xtbdft.py to reflect where xtb is located.\nExiting...")
    exit()
if not (os.path.exists(os.path.expanduser(xtbPath + "/bin/crest"))):
    print("Error: crest does not exist at " + xtbPath + "/bin/crest \n Please modify line 8 of xtbdft.py to reflect where crest is located.\nExiting...")
    exit()

realpath = os.path.realpath(__file__)
calcName=os.path.basename(os.getcwd())
startPath=os.getcwd()

def run_crest(file,chrg,uhf,**kwargs):
    TS=False
    TS=kwargs.get('TS')
    otherParams=kwargs.get('otherParams',[""])
    paramString=""
    for p in otherParams:
      paramString += p
    tsString = ""
    if TS == True:
        tsString = "-nozs -cinp .xcontrol"
        #need to ignore -cbonds here, or else crest will freeze
        params = paramString.split('-')
        paramString=""
        for p in params:
            if p.find("cbonds") == 0:
                print("CBonds for TS search already implemented manually to combine constraint files. Combining TS constraints and cbond constraints can result in frozen CREST calcs")
                pass
            elif p != '':
                paramString += "-" + p
                
    print("otherParams = " + paramString)
    script=open("submit_crest_{0}.sh".format(calcName), "w")
    script.write(msubHeader+"""
export OMP_STACKSIZE=1G
#source {0}/share/xtb/config_env.bash
ulimit -s unlimited
cd {1}
{2}/bin/crest {3} --chrg {4} --uhf {5} -T {6} {7} {8}> crest.out""".format(xtbPath,os.getcwd(),xtbPath,file,chrg,uhf,NP,tsString,paramString))
    script.close()
    #calcID = subprocess.run(['msub','submit_crest_{0}.sh'.format(calcName)], stdout=subprocess.PIPE).stdout.decode('utf-8')
    calcID = subprocess.check_output((qsub + ' submit_crest_{0}.sh'.format(calcName)).split()).decode('utf-8').replace("\n","")
    os.system("echo \"calcID is {0}\" | tee {0}.calcID".format(calcID))
    return calcID

def run_nwchem(chrg,uhf,calcType,xc,bs,cutoff,**kwargs):
    solv=kwargs.get("solv","decane") #decane as default solvent--solvent must have a COSMO shortname programmed into NWChem
    calcName2 = calcName +"_"+calcType
    if calcType == "crude" or calcType == "crudeTS":
        xc1,xc2 = xc[0],xc[1]
        bs1,bs2 = bs[0],bs[1]
    elif calcType == "refine" or calcType == "refineTS":
        xc1,xc2 = xc[2],xc[3]
        bs1,bs2 = bs[2],bs[3]     
        xyz="../minimum_lowest.xyz"
        if calcType == "refine":
            optType="optimize"
        else:
            optType="saddle"
    if calcType == "crudeTS" or calcType == "refineTS":
        #atomNo1 = kwargs.get('atom1')
        #atomNo2 = kwargs.get('atom2')
        atomNos = kwargs.get('atomNos')
        zcoordTxt = ""
        bonds = []
        atomPairs = list(combinations(atomNos,2))
        for atomPair in atomPairs:
            bonds = bonds + genBondsData(file,atomPair[0],atomPair[1],"wbo")
        os.system("echo constraining bonds: " + str(bonds))
        for bond in bonds:
            zcoordTxt += "  bond {0} {1} constant\n".format(bond[0],bond[1])
   
    input=open("{0}.nw".format(calcName2),"w")
    input.write("""
memory total 5200 mb
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
basis bs-fitting spherical
  * library "Weigend Coulomb Fitting"
end
""".format(calcType,chrg,uhf+1,xc1,bs1))
    if bs2 != "":
        input.write("""basis bs2 spherical
  * library {0}
end
""".format(bs2))
    if calcType == "crude" or calcType == "crudeTS":
        input.write("""dft
    grid fine
    convergence energy 1d-6
end
driver
    gmax 0.002 ; grms 0.0003 ; xrms 1 ; xmax 1
    maxiter 400         
end
""")
        numConfs=readConfs(cutoff)
        os.system("echo There are {0} conformers with deltaE_xtb < {1} kcal/mol".format(numConfs,cutoff))
        for i in range(1,numConfs+1):
            input.write("""###conf {0}
geometry units angstroms noautosym
    load frame {0} ../crest_conformers.xyz
end
set "ao basis" bs1
dft
    vectors input atomic output bs1.mos
end""".format(i))
            if calcType == "crudeTS":
                input.write("""
geometry adjust # fix reaction coordinate (bond length)
    zcoord
        {0}
    end
end""".format(zcoordTxt))
            input.write("""
task dft optimize
task shell "echo @ conf {0} geometry opt complete"
""".format(i))
            if(xc2 != ""):
                input.write("""
### single point at {0}/{1}//{2}/{3}
set "ao basis" bs2
dft
  xc {0}
  vectors input atomic output bs2.mos
end
task dft energy
""".format(xc2,bs2,xc1,bs1))
    elif calcType == "refine" or calcType == "refineTS":
        input.write("""dft
    grid fine
    convergence energy 1d-8
end
driver #tightopt criteria from Orca
    gmax 0.0001 ; grms 0.00003 ; xrms 0.0006 ; xmax 0.001
    maxiter 1000  
    sadstp 0.03
end
task shell "mkdir -p opt"
geometry units angstroms noautosym
    load {0}
end
set "ao basis" bs1
set "cd basis" bs-fitting
""".format(xyz))
        if calcType == "refineTS":
            input.write("""geometry adjust #fix reaction coordinate (bond)
  zcoord
    {0}
  end
end

task shell "echo @starting constrained opt"
task dft optimize

driver
  maxiter 0
  xyz opt/optC
end
task dft optimize ignore

geometry noautoz #unfix reaction coordinate, switch to cartesian opt for saddle pt search
  load opt/optC-000.xyz
end
driver
  clear
  maxiter 600
  sadstp 0.02
  xyz opt/optS
end
task shell "echo @starting saddle optimization" 
""".format(zcoordTxt))        
        input.write("""
task dft {0}

driver
  maxiter 0
  xyz opt/optF
end
task dft optimize ignore

unset "cd basis" #remove cd basis to enable anal. freq. calc.
task shell "echo @starting vibrational calculation"
task dft freq
set "cd basis" bs-fitting

task shell "echo @starting {2} solvated calc"
cosmo
  do_cosmo_smd true
  do_gasphase false
  solvent {2}
end
task dft energy
unset cosmo:* #need this line between every solvent in a multi-solvent calculation
cosmo
  off
end

task shell "echo @starting single point calculation"
set "ao basis" bs2
dft
  xc {1}
  vectors input atomic output bs2.mos
end
task dft energy
""".format(optType,xc2,solv))    
    else:
        print("NWChem must have a calctype: crude, refine, crudeTS, or refineTS")
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
    calcID = subprocess.check_output((qsub + ' submit_nwchem_{0}.sh'.format(calcName2)).split()).decode('utf-8').replace("\n","")
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
        if ("CREST terminated" in subprocess.check_output(['tail','-2','crest.out']).decode('utf-8')):
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
        if ("Total times" in subprocess.check_output(['tail','-2','nwchem.out']).decode('utf-8')):
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

def clean_crest():
    os.system("rm -vrf NORM* MRMSD METADYN*")

def clean_nwchem():
    os.system("rm -vf *.2ceri* *.b *.b^-1 *.c *.calcID *.cdfit *.db *.drv.hess *.gridpts.* *.p *.zmat submit*.sh*")

def goodvibes(file):
    os.system("echo starting goodvibes...")
    os.system("python3 {0} --invert -100 {1}".format(goodvibesPy,file))

def parseArgs():
    import argparse
    import sys
    print("Argument list:",str(sys.argv))
    parser=argparse.ArgumentParser()
    parser.add_argument("-chrg", action="store",dest="chrg",default=default_chrg,type=int)
    parser.add_argument("-uhf", action="store",dest="uhf",default=default_uhf,type=int)
    parser.add_argument("-xc",action="store",dest="xc",default=default_xc,type=str)
    parser.add_argument("-bs",action="store",dest="bs",default=default_bs,type=str)
    parser.add_argument("-cutoff",action="store",dest="cutoff",default=default_cutoff,type=float)
    parser.add_argument("-mode", action="store",dest="mode",default=["autoConf"],type=str,nargs='+')
    parser.add_argument("-other", action="store",dest="otherParams",default="",type=str,nargs='+')
    parser.add_argument("-solv", action="store",dest="solv",default="decane",type=str)
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
        print ("Error: -bs and -xc must be comma-delimited strings of exactly four items. Omitting them will default to them to: \n-xc b3lyp,,b3lyp,b3lyp \n-bs def2-svp,,def2-svp,def2-tzvp")
        sys.exit(1)
    print("mode = {}".format(results.mode))
    mode = results.mode[0]
    if mode == "autoConf":
        pass
    elif mode == "autoTS" and results.mode[1].isdigit() and results.mode[2].isdigit():
        try:
            float(results.mode[3])
        except:
            print("could not convert {0} to a float".format(results.mode[3]))
    elif mode == "TSconf" and len(results.mode[1:]) >= 2:
        for a in results.mode[1:]:
            try:
                int(a)
            except:
                print("could not convert {} to an int for atom index".format(a))
    else:
        print("Error: valid calls are '-mode autoConf' or '-mode autoTS atom1 atom2 finalDistance' or '-mode TSconf atom1 [...] atomN'")
        sys.exit(1)
    return results.filename,results.chrg,results.uhf,xc,bs,results.cutoff,mode,results.mode[1:],results.otherParams,results.solv
    
def checkEnv():
    #check that required shell scripts are executable
    #subprocess.check_output(['chmod','+x',realpath+"/calcInp_svP.sh"])
    #subprocess.check_output(['chmod','+x',realpath+"/calcInp_svp_final.sh"])
    #subprocess.check_output(['chmod','+x',realpath+"/pes_parse.py"])
    pass

def autoConf(file,chrg,uhf,xc,bs,cutoff,otherParams,solvent):
    print("otherParams = " + str(otherParams))
    if otherParams == ['skipCrest']:
        os.system("cp {0} minimum_lowest.xyz".format(file))
    else:
        if (os.path.isfile("crest.out")) and ("CREST terminated" in subprocess.check_output(['tail','-1','crest.out']).decode('utf-8')):
            print("CREST already complete, skipping to NWChem Crude Refinement")
         #run crest on input .xyz file
        else:
            os.system("echo {0}: running XTB pre-opt".format(datetime.now()))
            os.system("xtb {0} --chrg {1} --uhf {2} --opt | tee xtb_sp.out".format(file,chrg,uhf))
            os.system("echo {0}: running CREST".format(datetime.now()))
            calcID_crest = run_crest("xtbopt.xyz",chrg,uhf,otherParams=otherParams)
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
    calcID_nwchem = run_nwchem(chrg,uhf,"refine",xc,bs,cutoff,solv=solvent)
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

def autoTS(file,chrg,uhf,xc,bs,mode,TSparams,otherParams,solvent):
    from shutil import copyfile
    atomNos=[]
    if mode == "autoTS":
        [atomNo1,atomNo2,finalScanDistance] = TSparams
        atomNos=[atomNo1,atomNo2]
        #set up XTB scan
        os.system("echo {0}: setting up XTB xconstrains file".format(datetime.now()))
        startScanDistance=getBondDistance(file,atomNo1,atomNo2)
        barrier1=xtbScan(file,chrg,uhf,atomNo1,atomNo2,startScanDistance,finalScanDistance,"forward")
        xtbRxnPath(file,"end.xyz",chrg,uhf)

    elif mode == "TSconf":
        atomNos = TSparams
        os.system("echo TSparams: " + str(TSparams))
        os.system("cp {} TSguess.xyz".format(file))

    #run xtb single-point on TS guess, in order to determine if any hydrogen atoms are involved in the transition state, what the bonded neighbors of the hydrogen atom are
    os.system("mkdir -p xtb_sp")
    os.chdir("xtb_sp")
    os.system("xtb {0} --chrg {1} --uhf {2} | tee xtb_sp.out".format("../TSguess.xyz",chrg,uhf))
    os.chdir("../")
    
    if otherParams == ['skipCrest']:
        os.system("cp TSguess.xyz minimum_lowest.xyz")
    else:
        if (os.path.isfile("crest.out")) and ("CREST terminated" in subprocess.check_output(['tail','-2','crest.out']).decode('utf-8')):
            print("CREST already complete, skipping to NWChem Crude Refinement")
        else:
            #run crest on TSguess.xyz, constraining bonds/atoms of interest
            os.system("echo {0}: setting up CREST for TS guess files".format(datetime.now()))
            crestTS("TSguess.xyz",chrg,uhf,mode,atomNos,otherParams)

        #run NWChem constrained optimization on Crest ensemble
        os.system("mkdir nwchem")
        os.chdir("nwchem")
        os.system("echo {0}: running NWChem crude Constrained Opt on crest ensemble".format(datetime.now()))
        calcID_nwchem = run_nwchem(chrg,uhf,"crudeTS",xc,bs,cutoff,atom1=atomNo1,atom2=atomNo2)
        os.system("echo {0}: tracking NWChem...".format(datetime.now()))
        track_nwchem(calcID_nwchem)
        clean_nwchem()
    
        #determine lowest energy conformer and refine using NWChem
        os.system("echo {0}: running pes_parse subroutine".format(datetime.now()))
        pes_parse("nwchem.out")


    os.system("mkdir refine")
    os.chdir("refine")
    os.system("echo {0}: running NWChem refine Constrained Opt, TS opt, frequency calc, and single-point energy evaluation".format(datetime.now()))
    calcID_nwchem = run_nwchem(chrg,uhf,"refineTS",xc,bs,cutoff,atom1=atomNo1,atom2=atomNo2,solv=solvent)
    os.system("echo {0}: tracking NWChem refine".format(datetime.now()))
    track_nwchem(calcID_nwchem)
    clean_nwchem()

    #optional: script goodvibes correction
    os.system("echo {0}: running GoodVibes.py".format(datetime.now()))
    goodvibes("nwchem.out")   
    
    
def xtbRxnPath(file1,file2,chrg,uhf):
    dir = "xtbRxnPath"
    if (os.path.isfile(dir + "/xtb_rxnpath.out")):
        last_line = ""
        with open(dir + "/xtb_rxnpath.out") as f:
            for line in f:
                if line != "\n":
                    last_line = line
            print("***" + last_line)
        if last_line.split()[0] == "*" and last_line.split()[-1] == "speedup":
            print("xtb reaction path already completed...reading existing results...")
            return
  
    os.system("mkdir -p " + dir)
    os.chdir(dir)
    os.system("echo {0}: running XTB reaction path".format(datetime.now()))
    input=open("path.inp","w")
    input.write("""$path
   nrun=1
   npoint=25
   anopt=10
   kpush=0.003
   kpull=-0.015
   ppull=0.05
   alp=1.2
$end""")
    os.system("xtb ../{0} --chrg {1} --uhf {2} --opt".format(os.path.basename(file1),chrg,uhf))
    os.system("xtb ./{0} --path ../{3} --chrg {1} --uhf {2} -I path.inp | tee xtb_rxnpath.out".format("xtbopt.xyz",chrg,uhf,file2))
    os.chdir("../")
    os.system("cp {0}/xtbpath_ts.xyz ./TSguess.xyz".format(dir))
    
    
def xtbScan(file,chrg,uhf,atomNo1,atomNo2,startScanDistance,finalScanDistance,direction): 
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        pass
    steps=100
    dir = "scan_" + str(atomNo1) + "_" + str(atomNo2)
    #test if this scan has already been done
    if (os.path.isfile(dir + "/xtbopt.xyz")):
        print("scan already completed...reading existing results...")
        os.chdir(dir)
    else:    
        os.system("mkdir -p " + dir)
        os.system("cp " + file + " " + dir + "/")
        f=open(dir + "/xcontrol","w+")
        command = "$constrain \n  force constant = 1 \n  distance: " + str(atomNo1) + "," + str(atomNo2) + "," + str(startScanDistance) + "\n$scan\n  1: " + str(startScanDistance) + "," + str(finalScanDistance) + "," + str(steps) + " \n$opt \n  maxcycle=10 \n$end\n"
        f.write(command)
        f.close()
        os.chdir(dir)
        os.system("echo {0}: running XTB".format(datetime.now()))
        os.system("xtb {0} --chrg {1} --uhf {2} --opt veryfine --input xcontrol | tee opt.out".format(os.path.basename(file),chrg,uhf))
        #os.system("mkdir -p product")
        #os.chdir("product")
        #os.system("xtb ../xtbopt.xyz --chrg {0} --uhf {1} --opt veryfine | tee opt.out".format(chrg,uhf))
        #os.system("cp xtbopt.xyz ../../end.xyz")
        #os.chdir("../")
        os.system("cp xtbopt.xyz ../end.xyz")
            
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
    try:
        plt.scatter(bondLengths,energies)
        plt.xlabel(r"Bond Distance ($\AA$)")
        plt.ylabel("relative E (kcal/mol)")
        plt.savefig("scanPES.png")
    except:
        pass
    maxEnergy = energies[0]
    maxStruct = []
    maxIndex = 0
    rising=True
    for i in range(1,len(energies)):
        #if the scan was rising, but then experiences a sudden >2.5 kcal/mol drop, then there was a large geometric rearrangement
        if energies[i] > maxEnergy or (energies[i] - energies[i-1] < -2.5 and rising == True):
            maxEnergy = energies[i]
            maxStruct = structures[i]
            maxIndex = i
        if energies[i] - energies[i-1] > 0:
            rising = True
        else: 
            rising = False
    f=open("./TSguess.xyz","w+")
    f.write(str(totalAtoms) + "\n TS guess, E(xtb) = " + str(maxEnergy) + " \n")
    for coord in maxStruct:
        f.write(" ".join(coord) + "\n")
    f.close()
    os.system("cp TSguess.xyz ../")
    os.chdir("..")
    return maxEnergy

def crestTS(file,chrg,uhf,mode,atomNos,otherParams): 
    os.system("echo @@@ atomNos: {}".format(atomNos))
    #otherParams can be other commandline arguments for crest such as "-cbonds 0.1"
    from itertools import combinations
    dir = "crestTS"
    os.system("mkdir -p " + dir)
    os.system("cp " + file + " " + dir + "/")
    os.system("cp xtb_sp/wbo " + dir + "/") 
    os.chdir(dir)
    f = open(".xcontrol","w+")
    f.write("$constrain\n  force constant=1.0\n")
    bonds = []
    atomPairs = list(combinations(atomNos,2))
    for atomPair in atomPairs:
        bonds = bonds + genBondsData(file,atomPair[0],atomPair[1],"wbo")
    os.system("echo constraining bonds: " + str(bonds))
    for bond in bonds:
        f.write("  distance: {0}, {1}, {2}\n".format(bond[0],bond[1],bond[2]))
        #if -cbonds is used, for some reason the CREST calculation fails unless the .xcontrains file is streamlined manually
    paramString=""
    for p in otherParams:
      paramString += p
    params = paramString.split('-')
    paramString=""
    for p in params:
        if p.find("cbonds") == 0:
            print("CBonds for TS search being implemented manually to combine constraint files. Combining TS constraints and cbond constraints can result in frozen CREST calcs")
            if p.split()[1].replace('.','',1).isdigit():
                f.write("   force constant=" + p.split()[1] + "\n")
            f2 = open("wbo","r")
            lines = f2.readlines()
            for line in lines:
                m = line.split()[0]
                n = line.split()[1]
                f.write("  distance: {}, {}, {}\n".format(m, n, getBondDistance(file,m,n)))
    f.write("$end")
    f.close()
    
    os.system("echo {0}: running CREST".format(datetime.now()))
    run_crest(file,chrg,uhf,TS=True,otherParams=otherParams)
    track_crest("place","holder")
    clean_crest()
    return True

# Error: constraining atomNo1 and atomNo2 and -cbonds would result in crashing CREST calculations
def crestTS_old(file,chrg,uhf,atomNo1,atomNo2,otherParams): 
    #otherParams can be other commandline arguments for crest such as "-cbonds 0.1"
    dir = "crestTS"
    os.system("mkdir -p " + dir)
    os.system("cp " + file + " " + dir + "/")
    os.system("cp xtb_sp/wbo " + dir + "/") 
    os.chdir(dir)
    
    #if either of the contrained atoms is a hydrogen atom, then all the neighbors of that atom must also be fixed
    os.system("crest {} --constrain {},{}{}".format(file,atomNo1,atomNo2,hydrogenNeighbors(file,atomNo1,atomNo2,"wbo")))
    #os.system("crest {} --constrain {},{}".format(file,atomNo1,atomNo2))
    #check every 1 min if .xcontrol.sample was made
    while True:
        if (os.path.isfile(".xcontrol.sample")):
            print(os.getcwd()+"\.xcontrol.sample generated and moving to .xcontrol")
            os.system("mv .xcontrol.sample .xcontrol")
            break
        time.sleep(1*60)
        print("waiting for .xcontrol.sample to be generated...")
    os.system("echo {0}: running CREST".format(datetime.now()))
    run_crest(file,chrg,uhf,TS=True,otherParams=otherParams)
    track_crest("place","holder")
    return True

def genBondsData(xyz,atomNo1,atomNo2,wbo):
    bonds = []
    bonds.append([atomNo1,atomNo2,getBondDistance(xyz,atomNo1,atomNo2)])
    f = open(xyz,"r")
    lines = f.readlines()
    atom1=lines[1+int(atomNo1)]
    atom2=lines[1+int(atomNo2)]
    f.close()
    f = open(wbo,"r")
    lines = f.readlines()

    #if atom1 or atom2 are hydrogen atoms, looks for other bonds that may need to be constrained
    if atom1.split()[0].lower() == 'h':
        for line in lines:
            if str(atomNo1) in line.split() and str(atomNo2) not in line.split():
                bond=line.split()
                #only add constraint if the bond is weak
                if float(bond[-1]) < 0.4:
                    bonds.append([bond[0],bond[1],getBondDistance(xyz,bond[0],bond[1])])
    if atom2.split()[0].lower() == 'h':
        for line in lines:
            if str(atomNo2) in line.split() and str(atomNo1) not in line.split():
                bond=line.split()
                #only add constraint if the bond is weak
                if float(bond[-1]) < 0.4:
                    bonds.append([bond[0],bond[1],getBondDistance(xyz,bond[0],bond[1])])
    print("constraining the following bonds: " + str(bonds))
    return bonds

def hydrogenNeighbors(xyz,atomNo1,atomNo2,wbo):
    f = open(xyz,"r")
    lines = f.readlines()
    atom1=lines[1+int(atomNo1)]
    atom2=lines[1+int(atomNo2)]
    print(atom1 + " " + atom2)
    f.close()
    f = open(wbo,"r")
    lines = f.readlines()
    hydrogenNeighbors=[]
    if atom1.split()[0].lower() == 'h':
        for line in lines:
            if str(atomNo1) in line.split():
                neighbor=line.split()[:-1]
                neighbor.remove(str(atomNo1))
                hydrogenNeighbors.append(neighbor[0])
    if atom2.split()[0].lower() == 'h':
        for line in lines:
            if str(atomNo2) in line.split():
                neighbor=line.split()[:-1]
                neighbor.remove(str(atomNo2))
                hydrogenNeighbors.append(neighbor[0])    
    if hydrogenNeighbors == []:
        print("found no hydrogen atom neighbors")
        return ""
    else:
        print("found hydrogenNeighbors: " + str(hydrogenNeighbors))
        return ","+str(hydrogenNeighbors)[1:-1].replace(" ","")
    

def sharedNeighbor(file,atomNo1,atomNo2):
#reads a wbo output file from GFN2-XTB, and returns index of shared neighbor atom (metal center) between atomNo1 and atomNo2
    f=open(file,"r")
    lines = f.readlines()
    atom1Neighbors=[]
    atom2Neighbors=[]
    overlap=0
    for line in lines:
        atoms=line.split()[:2]
        if atoms[0] == str(atomNo1) and atoms[1] != str(atomNo2):
            atom1Neighbors.append(atoms[1])
        elif atoms[1] == str(atomNo1)  and atoms[0] != str(atomNo2):
            atom1Neighbors.append(atoms[0])
        elif atoms[0] == str(atomNo2) and atoms[1] != str(atomNo1):
            atom2Neighbors.append(atoms[1])
        elif atoms[1] == str(atomNo2)  and atoms[0] != str(atomNo1):
            atom2Neighbors.append(atoms[0])
    for i in atom1Neighbors:
        if i in atom2Neighbors:
            overlap = int(i)
    print("Atom 1 neighbors =" + str(atom1Neighbors))
    print("Atom 2 neighbors =" + str(atom2Neighbors))
    print("overlap =" + str(overlap))
    return overlap
            

def atomList(file,atomNo1,atomNo2):
    with open(file) as f:
        firstline = f.readlines()[0].rstrip()
    totalAtoms=int(firstline)
    print(totalAtoms)
    a=""
    inRange=False
    for n in range(1,totalAtoms+1):
        if n != atomNo1 and n != atomNo2:
            if inRange == False:
                a += "," + str(n)
                inRange = True
            elif inRange == True and n == totalAtoms:
                a += "-" + str(n)
        elif n == atomNo1 or n == atomNo2:
            if inRange == True:
                a += "-" + str(n-1)
                inRange = False
    return a 

if __name__ == "__main__":
    pid = os.getpid()
    os.system("touch " + str(pid))
    print("Reminder on how to run xtb-nwchem.py in headless mode:")
    print("nohup python3 ~/path/to/this/file/xtbdft.py geom.xyz [-chrg int] [-uhf int] > yourOutputFile.out &")
    pid = os.getpid()
    print("Process id = " + str(pid))
    file,chrg,uhf,xc,bs,cutoff,mode,TSparams,otherParams,solv=parseArgs()
    checkEnv()
    if (mode == "autoConf"):
        autoConf(file,chrg,uhf,xc,bs,cutoff,otherParams,solv)
    elif (mode == "autoTS" or mode == "TSconf"):
        autoTS(file,chrg,uhf,xc,bs,mode,TSparams,otherParams,solv)
    
