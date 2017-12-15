import sys
import numpy as np
from StringIO import StringIO
import argparse
import math
import os
import re

#convert slim output files (sim*.out) to phased bgl files for iHS input
#create ancestral allele, gmap files for iHS input
#write .sh files to perform iHS prep & scan for each sim replicate

#put directory of out files to use into command line
if len(sys.argv) == 3:
	headdir = sys.argv[1]
	rate = float(sys.argv[2])
else:
    sys.exit("input: python slim_to_bgl_for_iHS.py /path/to/slim/outputFiles/ recRate(r)")

# headdir = os.getcwd() + "/"

#make output directory
files = os.listdir(headdir)
outpath = headdir + "iHS/"
if not os.path.exists(outpath):
	os.mkdir(outpath)

#prep commands file
prepPath = outpath + "iHS_prep_commands.sh"
if os.path.isfile(prepPath):
	#append to existing file
	prepFile = open(prepPath,'a')
else:
	prepFile = open(prepPath,'w')

#scan commands file
scanPath = outpath + "iHS_scan_commands.sh"
if os.path.isfile(scanPath):
	#append to existing file
	scanFile = open(scanPath,'a')
else:
	scanFile = open(scanPath,'w')


#print all the files and directories
for file in files:
	#sim files
	if re.search('out$',file):
		#get sim # (filename format: simXX.out)
		simNum = file[3:len(file)-4]
		
		base = file[0:len(file)-4]
		
		#skip if bgl already exists
		if os.path.exists(outpath + base + ".phased"):
			continue
		
		#
		inFile = open(headdir + file, 'r')
		section = "o"

		#store mutations in a dictionary of lists
		mutList = []
		#store genomes in dictionary of lists
		inds = {}

		for line in inFile:
			line = line.rstrip('\n')
			if re.match('Mutations:',line):
				section = "m"
				continue
			elif re.match('Genomes:',line):
				section = "g"
				continue

			#if in mutations section, store mut info
			if (section == "m"):
				data = line.split()
				#make sure this is a mutation line
				if (data[0].isdigit()):
					#array [id,pos]
					mutList.append([data[0],int(data[3])])

			#if in genomes section, store genome data
			if (section == "g"):
				data = line.split()
				#make sure this is a genome line
				if (re.match('p1:',line)):
					#genome lines: p1:ID A 2 3 5 ... (mut IDs)
					indID = data[0].split(':')[1]
					mutsIn = data[2:len(data)]
					inds[indID] = mutsIn
					
		#if no muts or inds were stored, then skip this sim
		if not dict(inds):
			continue
				
		#sort mutlist by position
		mutList = sorted(mutList,key=lambda l:l[1])
		
		#print anc states file
		ancFile = open(outpath + base + ".ancStates",'w')
		ancFile.write("ID\tCHR\tPOS\tDER\tANC\n")

		#write gmap file using variant positions & rec rate
		gmap = open(outpath + base + ".gmap",'w')
		gmap.write("CHR\tID\tRHO\tPOS\tDER\tANC\n")
		
		#all variants are anc=A, der=C
		for v in range(0,len(mutList)):
			#skip writing variant if 2nd in a row w/ same pos (recurrent mutation)
			if v>0:
				if mutList[v][1]==mutList[v-1][1]:
					continue
			#write rsid chr(1) pos C A
			ancFile.write("rs"+mutList[v][0]+"\t1\t"+str(mutList[v][1])+"\tC\tA\n")
			rho = str(round(4*14041.61*rate*mutList[v][1],4))
			gmap.write("1\trs"+mutList[v][0]+"\t"+rho+"\t"+str(mutList[v][1])+"\tC\tA\n")
		
		ancFile.close()
		gmap.close()
		
		#open output bgl file
		outFile = open(outpath + base + ".phased",'w')
		#print header
		outFile.write("I id ")
		for i in range(0,len(inds)):
			outFile.write(str(i)+" ")
		outFile.write("\n")
		
		#for each variant, line starts w/ M rsID
		v=0
		while v<len(mutList):
			#if recurrent mutation, use same line for 2 variants
			if (v<len(mutList)-1) and (mutList[v][1]==mutList[v+1][1]):
				outFile.write("M rs"+mutList[v][0]+" ")
				for i in range(0,len(inds)):
					if (mutList[v][0] in inds[str(i)]) or (mutList[v+1][0] in inds[str(i)]):
						outFile.write("C ")
					else:
						outFile.write("A ")
				outFile.write("\n")
				#skip next mutation bc incorporated here
				v+=2
			
			else:
				#o.w. non recurrent	
				outFile.write("M rs"+mutList[v][0]+" ")
				#for each ind, ancestral allele (A) if not in ind's genome
				#else derived allele (C)
				for i in inds:
					if mutList[v][0] in inds[str(i)]:
						outFile.write("C ")
					else:
						outFile.write("A ")
				outFile.write("\n")
				v+=1
		outFile.close()
				
		
		#write pload file
		pload = open(outpath + base + ".pload",'w')
		pload.write("1 "+ outpath + base + ".phased\n")
		pload.close()
		
		#write prep commands
		prepFile.write("bsub -o out%J -e error%J WHAMM.pl --samplefile CHB_sims.sinfo --ploaderfile "+outpath + base + ".pload --iHS-prep "+outpath + base + ".ancStates "+outpath + base + ".gmap --out "+base+"\n")
		
		#write scan commands
		scanFile.write("bsub -o out%J -e error%J WHAMM.pl --iHS-scan "+base+"_iHS_chr1.info "+base+"_iHS_chr1.data --out "+base+"\n")
	
prepFile.close()
scanFile.close()









