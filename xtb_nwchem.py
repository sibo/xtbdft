#!/usr/bin/python
#
#running in headless mode: nohup python3 ~/path/to/this/file/xtb-nwchem.py geom.xyz [-chrg int] [-uhf int] [-xc str,str,str,str] [-bs str,str,str,str] > autoConf.out &
#
import os, sys, subprocess, time
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
#MSUB -l nodes=1:ppn={}
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
    script=open("submit_crest_{}.sh".format(calcName), "w")
    script.write(msubHeader+"""
export OMP_STACKSIZE=1G
source {}/Config_xtb_env.bash
ulimit -s unlimited
cd {}
{}/bin/crest {} --chrg {} --uhf {} -T {} > crest.out""".format(xtbPath,os.getcwd(),xtbPath,file,chrg,uhf,NP))
    script.close()
    #calcID = subprocess.run(['msub','submit_crest_{}.sh'.format(calcName)], stdout=subprocess.PIPE).stdout.decode('utf-8')
    calcID = subprocess.check_output(['msub','submit_crest_{}.sh'.format(calcName)]).decode('utf-8').replace("\n","")
    os.system("echo \"calcID is {}\" | tee {}.calcID".format(calcID,calcID))
    return calcID

def run_nwchem(chrg,uhf,calcType,xc,bs,cutoff):
    calcName2 = calcName +"_"+calcType
    if calcType == "crude":
        xc1,xc2 = xc[0],xc[1]
        bs1,bs2 = bs[0],bs[1]
    elif calcType == "refine":
        xc1,xc2 = xc[2],xc[3]
        bs1,bs2 = bs[2],bs[3]        
    input=open("{}.nw".format(calcName2),"w")
    input.write("""memory heap 200 mb stack 1000 mb global 2800 mb
start calc
title "B3LYP-D3/{}"
echo
charge {}
dft
    mult {}
    xc {}
    disp vdw 3
    iterations 200
    print low
    direct
end
basis bs1 spherical
  * library {}
end

basis "cd basis" spherical
  * library "Weigend Coulomb Fitting"
end
""".format(calcType,chrg,uhf+1,xc1,bs1))
    if bs2 != "":
        input.write("""basis bs2 spherical
  * library {}
end""".format(bs2))
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
        os.system("echo There are {} conformers with deltaE_xtb < {} kcal/mol".format(numConfs,cutoff))
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
    elif calcType == "refine":
        #run pes_parse.py on parent directory
        
        input.write("""dft
    grid fine
    convergence energy 1d-8
end
driver #tightopt criteria from Orca
    gmax 0.0001 ; grms 0.00003 ; xrms 0.0006 ; xmax 0.001
    maxiter 400         
end
geometry units angstroms noautosym
    load ../minimum_lowest.xyz
end
set "ao basis" bs1
task dft optimize
task shell "echo @starting vibrational calculation"
task dft freq numerical
task shell "echo @starting single point calculation"
set "ao basis" bs2
dft
  xc {0}
  vectors input atomic output bs2.mos
end
task dft energy
""".format(xc2))    
    else:
        print("NWChem must have a calctype: crude or refine")
        exit(1)
    
    #submit nwchem job to msub system
    script=open("submit_nwchem_{}.sh".format(calcName2), "w")
    script.write(msubHeader+"""
#export ARMCI_DEFAULT_SHMMAX=8192
module load mpich-x86_64
cd {}
mpirun -np {} nwchem_mpich {}.nw > nwchem.out
echo "{} (nwchem) submitted" """.format(os.getcwd(),NP,calcName2,calcName2))
    script.close()
    calcID = subprocess.check_output(['msub','submit_nwchem_{}.sh'.format(calcName2)]).decode('utf-8').replace("\n","")
    os.system("echo \"calcID is {}\" | tee {}.calcID".format(calcID,calcID))
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
            print(os.getcwd()+"/crest.out exists")
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
	opt=0
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

def autoSearch(f):
    #confirm existence of file f
    if(not os.path.isfile(f)):
        print("The file does not exist: " + f)
        sys.exit(1)
    
    #start CREST conf generation and crude ranking
    print("* starting CREST/GFN2-XTB conformational search on " + f)
    run_crest(f)
    calcMonitor = track_crest(f)
    if (calcMonitor != 1):
        print("!!! CREST failed! Printing end of crest.out:")
        os.system("tail -n 25 " + f)
        sys.exit(1)
    print("* finished CREST \n* starting NWChem crude re-optimization of crest conformers")
    run_nwchem(f)
    calcMonitor = track_nwchem(f)
    print("* finished NWChem crude re-optimization")

def goodvibes(file):
    os.system("echo starting goodvibes...")
    os.system("python3 {} --invert -100 {}".format(goodvibesPy,file))

def parseArgs():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("-chrg", action="store",dest="chrg",default=default_chrg,type=int)
    parser.add_argument("-uhf", action="store",dest="uhf",default=default_uhf,type=int)
    parser.add_argument("-xc",action="store",dest="xc",default=default_xc,type=str)
    parser.add_argument("-bs",action="store",dest="bs",default=default_bs,type=str)
    parser.add_argument("-cutoff",action="store",dest="cutoff",default=default_cutoff,type=float)
    parser.add_argument("filename",action="store")
    results=parser.parse_args()
    if (not os.path.isfile(results.filename)):
        print("Error: {} is not a file!".format(results.filename))
        sys.exit(1)
    elif results.filename[-4:] != ".xyz":
        print ("Error: Cartesian coord file must be .xyz")
        sys.exit(1)
    xc = results.xc.split(",")
    bs = results.bs.split(",")
    if (len(xc) != 4 or len(bs) != 4):
        print ("Error: -bs and -xc must be comma-delimited strings of exactly four items. Omitting them will default to them to: \n-xc b3lyp,,b3lyp,b3lyp \n-bs def2-sv(p),,def2-svp,def2-tzvp")
        sys.exit(1)
    return results.filename,results.chrg,results.uhf,xc,bs,results.cutoff
    
def checkEnv():
    #check that required shell scripts are executable
    #subprocess.check_output(['chmod','+x',realpath+"/calcInp_svP.sh"])
    #subprocess.check_output(['chmod','+x',realpath+"/calcInp_svp_final.sh"])
    #subprocess.check_output(['chmod','+x',realpath+"/pes_parse.py"])
    pass

if __name__ == "__main__":
    print("Reminder on how to run autoConf.py in headless mode:")
    print("nohup python3 ~/path/to/this/file/autoConf.py geom.xyz [-chrg int] [-uhf int] > autoConf.out &")
    pid = print(os.getpid())
    file,chrg,uhf,xc,bs,cutoff=parseArgs()
    checkEnv()
    
    #run crest on input .xyz file
    os.system("echo {}: running CREST".format(datetime.now()))
    calcID_crest = run_crest(file,chrg,uhf)
    os.system("echo {}: tracking CREST".format(datetime.now()))
    track_crest(file, calcID_crest)
    
    #run crude NWChem re-optimization
    os.system("mkdir nwchem")
    os.chdir("nwchem")
    os.system("echo {}: running NWChem Crude".format(datetime.now()))
    calcID_nwchem = run_nwchem(chrg,uhf,"crude",xc,bs,cutoff)
    os.system("echo {}: tracking NWChem Crude".format(datetime.now()))
    track_nwchem(calcID_nwchem)
    clean_nwchem()
    
    #determine the lowest energy conformer and refine using NWChem
    os.system("echo {}: running pes_parse subroutine".format(datetime.now()))
    pes_parse("nwchem.out")
    os.system("mkdir refine")
    os.chdir("refine")
    os.system("echo {}: running NWChem refine".format(datetime.now()))
    calcID_nwchem = run_nwchem(chrg,uhf,"refine",xc,bs,cutoff)
    os.system("echo {}: tracking NWChem refine".format(datetime.now()))
    track_nwchem(calcID_nwchem)
    clean_nwchem()
    
    #optional: script goodvibes correction
    os.system("echo {}: running GoodVibes.py".format(datetime.now()))
    goodvibes("nwchem.out")
    
