#!/usr/bin/python
#call: xtbTS(xyzFile.atomNo1,atomNo2,finalScanDistance)
#does: (1) sets up potential energy scan for distance atomNo1--atomNo2
#      (2) finds maximum structure along scan and runs Hessian calcluation (ohess)
#      (3) finds nearby TS (optts)
#
import os,sys,math
import matplotlib.pyplot as plt

print "Usage: xtbTS[.py] xyzFile atomNo1 atomNo2 finalScanDistance [xtbParams:--chrg x --uhf y]"
xyzFile = sys.argv[1]
atomNo1 = int(sys.argv[2])
atomNo2 = int(sys.argv[3])
finalScanDistance = float(sys.argv[4])
params = " "
for s in sys.argv[5:]:
	params += s + " "
#print params
#print "atomNo1 = " + str(atomNo1)
#print "atomNo2 = " + str(atomNo2)

def getBondDistance(xyzFile,atomNo1,atomNo2):
	f=open(xyzFile,"r")
	f1 = f.readlines()
	atomCoords=[]
	for l in f1[2:]:
		print l
		coordinates = l.split()
		atomCoords.append([float(coordinates[1]),float(coordinates[2]),float(coordinates[3])])
	deltaX = atomCoords[atomNo1-1][0] - atomCoords[atomNo2-1][0]
	deltaY = atomCoords[atomNo1-1][1] - atomCoords[atomNo2-1][1]
	deltaZ = atomCoords[atomNo1-1][2] - atomCoords[atomNo2-1][2]
	#print deltaX
	return math.pow(math.pow(deltaX,2)+math.pow(deltaY,2)+math.pow(deltaZ,2),0.5)

def parseXYZFile(xyzFile):
  f.open(xyzFile,"r")

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


steps=100
startScanDistance = getBondDistance(xyzFile,atomNo1,atomNo2)

def PESscan(direction, xyzFile, atomNo1, atomNo2, finalScanDistance, startScanDistance, steps):
	dir = "scan_" + str(atomNo1) + "_" + str(atomNo2) + "_" + direction
	os.system("mkdir -p " + dir)
	os.system("cp " + xyzFile + " " + dir + "/")
	f=open(dir + "/xcontrol","w+")
	command = "$constrain \n  force constant = 1 \n  distance: " + str(atomNo1) + "," + str(atomNo2) + "," + str(startScanDistance) + "\n$scan\n  1: " + str(startScanDistance) + "," + str(finalScanDistance) + "," + str(steps) + " \n$end\n"
	f.write(command)
	f.close()
	os.chdir(dir)
	os.system("xtb " + os.path.basename(xyzFile) + params + "--opt veryfine --input xcontrol | tee opt.out")
	f=open("xtbscan.log","r")
	f1=f.readlines()
	totalAtoms=int(f1[0])
	print "total atoms = " + str(totalAtoms)
	energies=[]
	structures=[]
	bondLengths=[]
	energy0 = float(f1[1].split()[1])*627.509
	for i in range(0,steps):
		energy = float(f1[i*(totalAtoms+2)+1].split()[1])*627.509
		energies.append(energy-energy0)
		print energy
		coords = []
		for line in f1[i*(totalAtoms+2)+2:(i+1)*(totalAtoms+2)]:
			coords.append(line.split())
		structures.append(coords)
		bondLengths.append(distance(coords[atomNo1-1],coords[atomNo2-1]))

	f=open("scan.csv","w+")
	for i in range(0,len(energies)):
		f.write(str(bondLengths[i]) + "," + str(energies[i]) + "\n")
	f.close()

	plt.scatter(bondLengths,energies)
	plt.xlabel(r"Bond Distance ($\AA$)")
	plt.ylabel("relative E (kcal/mol)")
	plt.savefig("scanPES.png")
	

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
	os.chdir("../")
	return maxEnergy,maxStruct,dir

(maxEnergy1,maxStruct1,dir1)=PESscan("forward", xyzFile, atomNo1, atomNo2, finalScanDistance, startScanDistance,steps)
#(maxEnergy2,maxStruct2,dir2)=PESscan("backward", dir1 + "/xtbopt.xyz", atomNo1, atomNo2, startScanDistance, finalScanDistance,steps)

#
#os.system("mkdir -p sp_ohess")
#if maxEnergy1 < maxEnergy2:
#	os.system("cp scan_forward/TSguess.xyz sp_ohess/")
#else:
#	os.system("cp scan_backward/TSguess.xyz sp_ohess/")
#os.chdir("./sp_ohess")
#f=open("xcontrol","w+")
#command = "$constrain \n  force constant = 1 \n  distance: " + str(atomNo1) + "," + str(atomNo2) + "," + str(getBondDistance("TSguess.xyz",atomNo1,atomNo2)) + "\n"
#f.write(command)
#f.close()
#os.system("xtb TSguess.xyz " + params + "--ohess veryfine --input xcontrol | tee ohess.out")
#os.system("xtb xtbopt.xyz " + params + "--optts veryfine | tee optts.out")


#scan for energy vs bond length
#plot energy vs bond length
#if  maxEnergy, then store xyz
#write maxEnergy XYZ file
#set up fixed ohess and optts calculations
