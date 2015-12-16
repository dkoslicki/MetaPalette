#This script will generate the tree plots from the inferred profile
import numpy as np
import h5py
import os, sys, shutil, subprocess, getopt
from Bio import Phylo
from Bio.Phylo import NewickIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from ete2 import Tree, faces, TreeStyle, COLOR_SCHEMES, TextFace, BarChartFace, CircleFace, AttrFace, NodeStyle, RectFace
import math
import os
import ClassifyPackage
import PlotPackage

outgroup = "Halobacterium_sp_DL1"

try:
	opts, args = getopt.getopt(sys.argv[1:],"hd:o:p:i:t:g:",["Help=", "DataDir=", "OutputFolder=", "ProfileDir=", "InputFileName=", "Taxon=", "Outgroup="])
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

kmer_sizes=[30,50]

#Check input taxon
if taxon!="genus" and taxon!="species":
	print("Error: taxon (-t) must be either species or genus. Value of " + taxon + " given.")
	sys.exit(2)
if not os.path.isdir(data_dir):
	print("Error: Data directory " + data_dir + " does not exist.")
	sys.exit(2)
if not os.path.isfile(os.path.join(data_dir,"Taxonomy.txt")):
	print("Error: Missing taxonomy file: %s" % os.path.join(data_dir,"Taxonomy.txt"))
	sys.exit(2)
if not os.path.isfile(os.path.join(profile_folder, input_file_basename+".profile")):
	print("Error: Missing profile file: " + os.path.join(profile_folder, input_file_basename+".profile"))
	print("Please run the Classify.py script and try again")
	sys.exit(2)
for kmer_size in kmer_sizes:
	if not os.path.isfile(os.path.join(profile_folder,input_file_basename+"-y"+str(kmer_size)+".txt")):
		print("Error: Missing file " + os.path.join(profile_folder,input_file_basename+"-y"+str(kmer_size)+".txt"))
		print("Please run the Classify.py script and try again.")
		sys.exit(2)
for kmer_size in kmer_sizes:
	if not os.path.isfile(os.path.join(data_dir,"CommonKmerMatrix-"+str(kmer_size)+"mers.h5")):
		print("Error: Missing file " + os.path.join(data_dir,"CommonKmerMatrix-"+str(kmer_size)+"mers.h5"))
		print("Please run the Train.py script (or download the pre-trained data) and try again.")
		sys.exit(2)

#Get name of the file of interest
input_file_basename = os.path.basename(input_file_name)

#Next, read in the taxonomy file
fid = open(os.path.join(data_dir,"Taxonomy.txt"),"r")
organism_names = []
tax_paths = []
for line in fid:
	temp = line.strip().split()[0]
	organism_names.append("_".join(temp.split("_")[1:]))
	tax_paths.append(line.strip().split()[2])

fid.close()

#Check if outgroup is in the taxonomy
if outgroup not in organism_names:
	print("Error: outgroup " + outgroup + " is not one of the organisms in " + os.path.join(data_dir,"Taxonomy.txt"))
	sys.exit(2)
else:
	outgroup_index = organism_names.index(outgroup)

#Read in the .profile file, concentrate on organisms that were actually classified to
fid = open(os.path.join(profile_folder, input_file_basename+".profile"))
species = list()
genera = list()
for line in fid:
	if line[0]!="@" and line[0]!="#" and line!="\n":
		name = line.strip().split()[1]
		tax_path = line.strip().split()[3]
		if name=="species":
			species.append(tax_path.split("|")[-1])
		elif name=="genus":
			genera.append(tax_path.split("|")[-1])

species = list(set(species))
genera = list(set(genera))

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

#Read find the basis corresponding to the particular taxa, split into taxa chunks, do the plot for each taxa
if taxon == "species":
	for specie in species:
		#select all the training organisms of this same species
		select_indicies = list()
		for tax_path in tax_paths:
			tax_names = tax_path.split("|")
			for tax_name in tax_names:
				if tax_name[0]=="s":
					temp = "_".join(tax_name.split("_")[3:])
					if temp == specie:
						select_indicies.append(tax_paths.index(tax_path))
		if len(select_indicies)>1: #Have to have more than two organisms to do the plot
			select_indicies.append(outgroup_index)
			CKM_matrices_reduced = list()
			CKM_matrices_reduced.append(CKM_matrices[0][select_indicies,:][:,select_indicies])
			CKM_matrices_reduced.append(CKM_matrices[1][select_indicies,:][:,select_indicies])
			organism_names_reduced = [organism_names[i] for i in select_indicies]
			Y_norms_reduced = list()
			Y_norms_reduced.append(Y_norms[0][select_indicies])
			Y_norms_reduced.append(Y_norms[1][select_indicies])
			print("Creating NJ tree and plot for " + specie)
			x = ClassifyPackage.Classify(organism_names_reduced, CKM_matrices_reduced, Y_norms_reduced)
			sum_x = sum(x)
			#Normalize the x vector
			x = map(lambda y: y/sum(x),x)
			outfile = os.path.join(output_folder, input_file_basename+"-"+specie+".png")
			outfilexml = os.path.join(output_folder, input_file_basename+"-"+specie+".xml")
			PlotPackage.MakePlot(x, organism_names_reduced, CKM_matrices_reduced[0], CKM_matrices_reduced[1], outgroup, outfile, outfilexml, sum_x)
elif taxon == "genus":
	for genus in genera:
		#select all the training organisms of this same species
		select_indicies = list()
		for tax_path in tax_paths:
			tax_names = tax_path.split("|")
			for tax_name in tax_names:
				if tax_name[0]=="g":
					temp = "_".join(tax_name.split("_")[3:])
					if temp == genus:
						select_indicies.append(tax_paths.index(tax_path))
		if len(select_indicies)>1: #Have to have more than two organisms to do the plot
			select_indicies.append(outgroup_index)
			CKM_matrices_reduced = list()
			CKM_matrices_reduced.append(CKM_matrices[0][select_indicies,:][:,select_indicies])
			CKM_matrices_reduced.append(CKM_matrices[1][select_indicies,:][:,select_indicies])
			organism_names_reduced = [organism_names[i] for i in select_indicies]
			Y_norms_reduced = list()
			Y_norms_reduced.append(Y_norms[0][select_indicies])
			Y_norms_reduced.append(Y_norms[1][select_indicies])
			print("Creating NJ tree and plot for " + genus)
			x = ClassifyPackage.Classify(organism_names_reduced, CKM_matrices_reduced, Y_norms_reduced)
			sum_x = sum(x)
			#Normalize the x vector
			x = map(lambda y: y/sum(x),x)
			outfile = os.path.join(output_folder, input_file_basename+"-"+genus+".png")
			outfilexml = os.path.join(output_folder, input_file_basename+"-"+genus+".xml")
			PlotPackage.MakePlot(x, organism_names_reduced, CKM_matrices_reduced[0], CKM_matrices_reduced[1], outgroup, outfile, outfilexml, sum_x)

#Do the classification
#x = ClassifyPackage.Classify(organism_names, CKM_matrices, Y_norms)

#Make the tree and export it###############
#outfile = os.path.join(output_folder, input_file_basename+"-testout.png")
#PlotPackage.MakePlot(x, organism_names, CKM_matrices[0], CKM_matrices[1], outgroup, outfile)
	
	
	
	

