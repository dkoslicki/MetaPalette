#This script will generate the tree plots from the inferred profile
import numpy as np
import h5py
from Bio import Phylo
from Bio.Phylo import NewickIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from ete2 import Tree, faces, TreeStyle, COLOR_SCHEMES, TextFace, BarChartFace, CircleFace, AttrFace, NodeStyle, RectFace
import math
import os
import ClassifyPackage
import PlotPackage


#Input
# Location of CKM files
# Location of Taxonomy file
# Location of y-file and .profile file (call it the ProfileDir)

#Options
# Either generate the tree for a specific taxa, or else do all the species that have inferred abundance

#Output
# The plots of the individual trees, in a specific directory, with file names given by taxa of interest

outgroup = "Halobacterium_sp_DL1"

try:
	opts, args = getopt.getopt(sys.argv[1:],"hd:o:p:i:t:g",["Help=", "DataDir=", "OutputFolder=", "ProfileDir=", "InputFileName=", "Taxon=", "Outgroup="])
except getopt.GetoptError:
	print 'Unknown option, call using: python Plot.py -d <DataDir> -o <OutputFolder> -p <ProfileDir> -i <InputFileName> -t <Taxon> -g <Outgroup>'
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print 'python Plot.py -d <DataDir> -o <OutputFolder> -p <ProfileDir> -i <InputFileName> -t <Taxon> -g <Outgroup>'
		sys.exit(2)
	elif opt in ("-d", "--DataDir"):
		data_dir = arg
	elif opt in ("-o", "--OutputFolder"):
		output_folder = arg
	elif opt in ("-p","--ProfileDir"):
		profile_folder = arg
	elif opt in ("-i", "--InputFileName"):
		input_file_name = arg
	elif opt in ("-t", "--Taxon"):
		taxon = arg
	elif opt in ("-g", "--Outgroup"):
		outgroup = arg

#Get name of the file of interest
input_file_basename = os.path.basename(input_file_name)

#Parse the taxonomy
if not os.path.isfile(os.path.join(data_dir,"Taxonomy.txt")):
	print("Error: Missing taxonomy file: %s" % os.path.join(data_dir,"Taxonomy.txt"))
	sys.exit(2)

#Next, read in the taxonomy file
fid = open(os.path.join(data_dir,"Taxonomy.txt"),"r")
taxonomy = []
for line in fid:
	taxonomy.append(line.strip().split()[0])

fid.close()

#Read in Y_norms
Y_norms = list()
for kmer_size in kmer_sizes:
	fid = open(os.path.join(profile_folder,input_file_basename+"-y"+str(kmer_size)+".txt"),'r')
	Y = fid.readlines()
	fid.close()
	Y = np.array(map(lambda y: float(y), Y), dtype=np.float64)
	Y_norms.append(Y)


#Read in CKMs
CKM_matrices = list()
for kmer_size in kmer_sizes:
	fid = h5py.File(os.path.join(data_dir,"CommonKmerMatrix-"+str(kmer_size)+"mers.h5"),'r')
	CKM_matrices.append(np.array(fid["common_kmers"][:,:], dtype = np.float64))

#Do the classification
x = ClassifyPackage.Classify(taxonomy, CKM_matrices, Y_norms)

#Make the tree and export it###############
outfile = os.path.join(output_folder, input_file_basename, "testout.png")
PlotPackage.MakePlot(x, taxonomy, CKM_matrices[0], CKM_matrices[1], outgroup, outfile)



#if taxon == "species":
#	#Read in the y30 file, find the basis, split into species chunks, do the plot for each species
#elif: taxon == "genus":
#	#Read in the y30 file, find the basis, split into genus chunks, do the plot for each species
#else:
#	#Select only the organisms of interest, get the basis, reduce the CKM matrices, then do plot

	
	
	
	
	
	
	
	

