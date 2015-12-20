#This is a python clone of Classify.jl
import os, sys, shutil, subprocess, getopt, math, copy
from itertools import *
from multiprocessing import Pool, freeze_support
import numpy as np
import h5py
import numpy.matlib
import scipy.optimize
import ClassifyPackage

#python /nfs1/Koslicki_Lab/koslickd/CommonKmers/Python/Scripts/Classify.py -d /nfs1/Koslicki_Lab/koslickd/CommonKmers/Python/Output -o /nfs1/Koslicki_Lab/koslickd/CommonKmers/Python/Output -i /nfs1/Koslicki_Lab/koslickd/CommonKmers/Python/Data/test-reads.fa -k default -j /home/pi/koslickd/jellyfish-2.2.3/bin/./jellyfish -q /nfs1/Koslicki_Lab/Backup/koslickd/CAMI/CommonKmers/src/QueryPerSeq/./query_per_sequence -Q C -y -x -t 48

#save_y = False
save_x = False
quality = 'C'
re_run = False
normalize = False

try:
	opts, args = getopt.getopt(sys.argv[1:],"hd:o:i:k:j:q:Q:xrnt:",["Help=", "DataDir=","OutputFolder=", "InputFile=", "Kind=", "JellyfishBinary=", "QueryPerSequenceBinary=", "Quality=", "SaveX=", "ReRun=", "Normalize=", "Threads="])
except getopt.GetoptError:
	print 'Unknown option, call using: python Classify.py -d <DataDir> -o <OutputFolder> -i <InputFile> -k <Kind> -j <JellyfishBinary> -q <QueryPerSequenceBinary> -Q <Quality> -x <SaveX Flag> -r <ReRun Flag> -n <Normalize Flag> -t <Threads>'
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print 'python Classify.py -d <DataDir> -o <OutputFolder> -i <InputFile> -k <Kind> -j <JellyfishBinary> -q <QueryPerSequenceBinary> -Q <Quality> -x <SaveX Flag> -r <ReRun Flag> -n <Normalize Flag> -t <Threads>'
		sys.exit(2)
	elif opt in ("-d", "--DataDir"):
		data_dir = arg
	elif opt in ("-o", "--OutputFolder"):
		output_folder = arg
	elif opt in ("-i","--InputFile"):
		input_file_name = arg
	elif opt in ("-k", "--Kind"):
		kind = arg
	elif opt in ("-j", "--JellyfishBinary"):
		jellyfish_binary = arg
	elif opt in ("-q", "--QueryPerSequenceBinary"):
		query_per_sequence_binary = arg
	elif opt in ("-Q", "--Quality"):
		quality = arg
#	elif opt in ("-y", "--SaveY"):
#		save_y = True 
	elif opt in ("-x", "--SaveX"):
		save_x = True 
	elif opt in ("-r", "--ReRun"):
		re_run = True 
	elif opt in ("-n", "--Normalize"):
		normalize = True
	elif opt in ("-t", "--Threads"):
		num_threads = int(arg) 

file_base_name = os.path.basename(input_file_name)
thresholds=[.90,.80,.70,.60,.50,.40,.30,.20,.10]
kmer_sizes = [30,50]

#Check args and file existence, etc.
if not os.path.isfile(input_file_name):
	print("Error: Input file name " + input_file_name + " does not exist.")
	sys.exit(2)
if not os.path.isdir(output_folder):
	print("Error: Output folder " + output_folder + " does not exist.")
	sys.exit(2)
if not os.path.isdir(data_dir):
	print("Error: Data directory " + data_dir + " does not exist.")
	sys.exit(2)
for kmer_size in kmer_sizes:
	if not os.path.isfile(os.path.join(data_dir,"CommonKmerMatrix-"+str(kmer_size)+"mers.h5")):
		print("Error: Missing Common Kmer Matrix file " + os.path.join(data_dir,"CommonKmerMatrix-"+str(kmer_size)+"mers.h5"))
		print("Please re-run Train.py or download the pre-trained data")
		sys.exit(2)
if not os.path.isfile(os.path.join(data_dir, "FileNames.txt")):
	print("Error, missing FileNames.txt file in " + data_dir)
	sys.exit(2)
if kind!="default" and kind!="sensitive" and kind!="specific":
	error("invalid kind (-k) option. Options are: default, sensitive, specific")
	sys.exit(2)
if not os.path.isfile(os.path.join(data_dir,"Taxonomy.txt")):
	print("Error: Missing taxonomy file: %s" % os.path.join(data_dir,"Taxonomy.txt"))
	sys.exit(2)

fid = open(os.path.join(data_dir,"FileNames.txt"),'r')
training_file_names = list()
for file in fid:
	training_file_names.append(os.path.basename(file.strip()))
	if not os.path.isfile(os.path.abspath(file.strip())):
		print("Error: File " + file + " given in " + os.path.join(data_dir,"FileNames.txt") + " but does not exist.")
		sys.exit(2)

#Form the sample jellyfish file
#Function
def count_kmers(file, kmer_size):
	extension = os.path.splitext(file)[1]
	if extension==".bz2" or extension==".bz":
		cmd = "bzip2 -c " + file + " | " + jellyfish_binary + " count /dev/fd/0 -m "+str(kmer_size)+" -t "+str(num_threads)+" -s 100M -C -Q "+str(quality)+" -o " + os.path.join(output_folder,file_base_name+"-"+str(kmer_size)+"mers.jf")
	elif extension==".gz" or extension==".z" or extension==".Z":
		cmd = "gunzip -c " + file + " | " + jellyfish_binary + " count /dev/fd/0 -m "+str(kmer_size)+" -t "+str(num_threads)+" -s 100M -C -Q "+str(quality)+" -o " + os.path.join(output_folder,file_base_name+"-"+str(kmer_size)+"mers.jf")
	else:
		cmd = jellyfish_binary + " count " + file + " -m "+str(kmer_size)+" -t "+str(num_threads)+" -s 100M -C -Q "+str(quality)+" -o " + os.path.join(output_folder,file_base_name+"-"+str(kmer_size)+"mers.jf")
	test = subprocess.check_output(cmd, shell = True)

#Function to form Y files
def form_y(file_name, kmer_size):
	cmd = query_per_sequence_binary +" "+ os.path.join(output_folder,file_base_name+"-"+str(kmer_size)+"mers.jf") + " " + os.path.join(data_dir,"Bcalms",file_name+"-30mers.bcalm.fa")
	res = subprocess.check_output(cmd, shell = True)
	return int(res)

def form_y_star(arg):
	return form_y(*arg)

#Form the Y files
Y_norms = list()
if re_run and all([os.path.join(output_folder,file_base_name+"-y"+str(kmer_size)+".txt") for kmer_size in kmer_sizes]):
	for kmer_size in kmer_sizes:
		fid = open(os.path.join(output_folder,file_base_name+"-y"+str(kmer_size)+".txt"),'r')
		Y = fid.readlines()
		fid.close()
		Y = np.array(map(lambda y: float(y), Y), dtype=np.float64)
		Y_norms.append(Y)
else:
	#Count the kmers
	for kmer_size in kmer_sizes:
		if os.path.isfile(os.path.join(output_folder,file_base_name+"-"+str(kmer_size)+"mers.jf")):
			print("Kmers already counted in file " + os.path.join(output_folder,file_base_name+"-"+str(kmer_size)+"mers.jf"))
			print("Skipping kmer counting step")
		else:
			count_kmers(input_file_name, kmer_size)
	#Form the y-files
	for kmer_size in kmer_sizes:
		pool = Pool(processes = num_threads)
		Y = np.array(pool.map(form_y_star, izip(training_file_names, repeat(kmer_size))), dtype=np.float64)
		pool.close()
		total_kmers = int(subprocess.check_output(jellyfish_binary + " stats " + os.path.join(output_folder,file_base_name+"-"+str(kmer_size)+"mers.jf"), shell = True).split()[5])
		Y_norm = Y/total_kmers
		Y_norms.append(Y_norm)
		fid = open(os.path.join(output_folder,file_base_name+"-y"+str(kmer_size)+".txt"),'w')
		for i in range(len(Y_norm)):
			fid.write(str(Y_norm[i])+"\n")
		fid.close()

#Load the common kmer matrices
CKM_matrices = list()
for kmer_size in kmer_sizes:
	fid = h5py.File(os.path.join(data_dir,"CommonKmerMatrix-"+str(kmer_size)+"mers.h5"),'r')
	CKM_matrices.append(np.array(fid["common_kmers"][:,:], dtype = np.float64))

#Do the classification
x = ClassifyPackage.Classify(training_file_names, CKM_matrices, Y_norms)

#Normalize the result
if normalize or (x.sum()>1):
	x=x/x.sum()

#Write x_file
fid = open(os.path.join(output_folder,file_base_name+"-x.txt"),'w')
for i in range(len(x)):
	fid.write(str(x[i]) + "\n")

fid.close()

#################################################################################################
#print("Doing format conversion")
#Convert to CAMI compatible output format
if kind=="sensitive":
	cutoff = .0001
else:
	cutoff = .001

#if not os.path.isfile(os.path.join(output_folder,file_base_name+"-x.txt")):
#	print("Error: Missing x-file: %s" % os.path.join(output_folder,file_base_name+"-x.txt"))
#	sys.exit(2)

#Read in the input file
#fid = open(os.path.join(output_folder,file_base_name+"-x.txt"),"r")
#input = []
#for line in fid:
#	input.append(float(line.strip()))
#fid.close()

input = [x[i] for i in range(len(x))]

#Next, read in the taxonomy file
fid = open(os.path.join(data_dir,"Taxonomy.txt"),"r")
taxonomy = []
for line in fid:
	taxonomy.append(line.strip().split()[2])

fid.close()

num_organisms = len(taxonomy)
support = np.where(np.array(input).flatten()> cutoff)[0]

#Do the non-hypothetical organisms first, populating the species and strains
#Then do the the hypothetical organisms, doing the LCA, populating some of the higher taxonomy levels
#Then starting at the bottom, increment (or create, as the case may be) the higher taxonomic levels
#Then write the output.
kmer_size = 30
fid = h5py.File(os.path.join(data_dir,"CommonKmerMatrix-"+str(kmer_size)+"mers.h5"),'r')
A = np.array(fid["common_kmers"][:,:], dtype = np.float64)
fid.close()
common_kmer_matrix_normalized = A/np.diag(A)
output_taxonomy = dict()
for support_index in support:
	if support_index < len(taxonomy): #If it's not a hypothetical organism
		if taxonomy[support_index] in output_taxonomy: #If it's in there, add to it.
			output_taxonomy[taxonomy[support_index]] = output_taxonomy[taxonomy[support_index]] + input[support_index]
		else:
			output_taxonomy[taxonomy[support_index]] = input[support_index]
	else: #it's a hypothetical organism, so do the LCA here
		corresponding_real_organism_index = support_index % num_organisms
		#fix zero guy
#		if corresponding_real_organism_index == 0:
#			corresponding_real_organism_index = num_organisms
		column = common_kmer_matrix_normalized[:,corresponding_real_organism_index];
		hyp_bin = int(math.floor((support_index)/num_organisms)); #This is the bin it belongs to 
		hyp_thresh = thresholds[hyp_bin-1]; #This is the corresponding threshold
		#LCA here
		temp = column;
		distances = abs(temp-hyp_thresh); #distances between column and thresh
		pos = np.argmin(distances) #This is the position of the nearest real organism to the threshold
		# If this distance is less than a bin's width, then do the LCA bit, otherwise go to a fixed taxonomy bit
		if kind=="default" or distances[pos]<=0.2: #Later, make this adaptable to given input bins
			#Now find the LCA between this organism and the hypothetical one
			hyp_taxonomy_split = taxonomy[corresponding_real_organism_index].split("|") #Hyp taxonomy based on corresponding organism
			candidate_LCA_taxnonmy_split = taxonomy[pos].split("|") #closest organism taxnonmy
			LCA = 0 #In case it only agrees to the kingdom level, and we need to declare LCA as a global variable
			for LCA_index in range(min([len(candidate_LCA_taxnonmy_split), len(hyp_taxonomy_split)])-1,-1,-1):
				if hyp_taxonomy_split[LCA_index] == candidate_LCA_taxnonmy_split[LCA_index]:
					LCA = LCA_index
					break
		else: #Fixed LCA based on bin
			hyp_taxonomy_split = taxonomy[corresponding_real_organism_index].split("|")
			candidate_LCA_taxnonmy_split = hyp_taxonomy_split
			LCA = 0;
			if kind=="sensitive":
				if hyp_thresh == .9:
					LCA = len(hyp_taxonomy_split)-1
				elif hyp_thresh == .8:
					LCA = len(hyp_taxonomy_split)-1
				elif hyp_thresh == .7:
					LCA = max([1, len(hyp_taxonomy_split)-1])-1
				elif hyp_thresh == .6:
					LCA = max([1, len(hyp_taxonomy_split)-1])-1
				elif hyp_thresh == .5:
					LCA = max([1, len(hyp_taxonomy_split)-2])-1
				elif hyp_thresh == .4:
					LCA = max([1, len(hyp_taxonomy_split)-2])-1
				elif hyp_thresh == .3:
					LCA = max([1, len(hyp_taxonomy_split)-2])-1
				elif hyp_thresh == .2:
					LCA = max([1, len(hyp_taxonomy_split)-3])-1
				elif hyp_thresh == .1:
					LCA = max([1, len(hyp_taxonomy_split)-4])-1
			elif kind=="specific":
				if hyp_thresh == .9:
					LCA = len(hyp_taxonomy_split)-1
				elif hyp_thresh == .8:
					LCA = len(hyp_taxonomy_split)-1
				elif hyp_thresh == .7:
					LCA = max([1, length(hyp_taxonomy_split)-1])-1
				elif hyp_thresh == .6:
					LCA = max([1, length(hyp_taxonomy_split)-2])-1
				elif hyp_thresh == .5:
					LCA = max([1, length(hyp_taxonomy_split)-2])-1
				elif hyp_thresh == .4:
					LCA = max([1, length(hyp_taxonomy_split)-2])-1
				elif hyp_thresh == .3:
					LCA = max([1, length(hyp_taxonomy_split)-2])-1
				elif hyp_thresh == .2:
					LCA = max([1, length(hyp_taxonomy_split)-3])-1
				elif hyp_thresh == .1:
					LCA = max([1, length(hyp_taxonomy_split)-4])-1
			else:
				error("invalid kind (-k) option. Options are: default, sensitive, specific")
		#Now update the taxonomy dictionary
		LCA_taxonomy = "|".join(candidate_LCA_taxnonmy_split[0:LCA+1])
		if LCA_taxonomy in output_taxonomy: #If it's in there, add to it.
			output_taxonomy[LCA_taxonomy] = output_taxonomy[LCA_taxonomy] + input[support_index]
		else:
			output_taxonomy[LCA_taxonomy] = input[support_index]

#Now sum up from coarser taxonomic ranks to higher ones (summing finer to coarser for each taxa level)
temp_dict = copy.deepcopy(output_taxonomy)
for taxonomic_level in range(1,8):
	for key in temp_dict: 
		split_key = key.split("|")
		if len(split_key)-1 == taxonomic_level:
			#loop through the higher taxonomic levels
			for higher_level in range((taxonomic_level),0,-1):
				higher_taxonomy = "|".join(split_key[0:higher_level])
				#update the dictionary value
				if higher_taxonomy in output_taxonomy: #If it's in there, add to it.
					output_taxonomy[higher_taxonomy] = output_taxonomy[higher_taxonomy] + output_taxonomy[key] #add the value at the base taxonomy
				else: #If not, make it equal to the base taxonomy value
					output_taxonomy[higher_taxonomy] = output_taxonomy[key]


output_file_handle = open(os.path.join(output_folder,file_base_name+".profile"),"w")
#Write the header
output_file_handle.write("# Taxonomic profile for file: %s\n" % os.path.abspath(input_file_name))
output_file_handle.write("@Version:0.9.1\n")
output_file_handle.write("@SampleID:%s\n" % file_base_name)
output_file_handle.write("@Ranks: superkingdom|phylum|class|order|family|genus|species|strain\n")
output_file_handle.write("\n")
output_file_handle.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")

taxa_names = [key for key in output_taxonomy]
taxa_names.sort(key=lambda x: len(x.split("|")))

for taxa_name in taxa_names:
	taxID = taxa_name.split("|")[-1].split("_")[2]
	rankAbvr = taxa_name.split("|")[-1].split("_")[0]
	if rankAbvr == "k":
		rank = "superkingdom"
	elif rankAbvr == "p":
		rank = "phylum"
	elif rankAbvr == "c":
		rank = "class"
	elif rankAbvr == "o":
		rank = "order"
	elif rankAbvr == "f":
		rank = "family"
	elif rankAbvr == "g":
		rank = "genus"
	elif rankAbvr == "s":
		rank = "species"
	elif rankAbvr == "t":
		rank = "strain"
	else:
		rank = "unknown"
	taxPath = map(lambda x: x.split("_")[2],taxa_name.split("|")) #Tax ID's
	taxPathSN = map(lambda x: "_".join(x.split("_")[3:]), taxa_name.split("|")) #Taxa names
	#If a Tax ID is repeated at a lower taxonomic rank, this means that that rank is missing, so let's just delete it.
	for i in range(len(taxPath)-1,1,-1):
		if i>=2:
			if taxPath[i] == taxPath[i-1]:
				taxPath[i] = ""
				taxPathSN[i] = ""
	#If it ends with a blank, we don't want to be including it
	if taxPath[-1] != "":
		output_file_handle.write("%s" % taxID)
		output_file_handle.write("\t")
		output_file_handle.write("%s" % rank)
		output_file_handle.write("\t")
		#Join back up the paths
		taxPath = "|".join(taxPath)
		taxPathSN = "|".join(taxPathSN)
		output_file_handle.write("%s" % taxPath)
		output_file_handle.write("\t")
		output_file_handle.write("%s" % taxPathSN)
		output_file_handle.write("\t")
		output_file_handle.write("%s" % str(100*output_taxonomy[taxa_name])) #Turn into a frequency (i.e. 90 not .9)
		output_file_handle.write("\n")

#Close the output file
output_file_handle.close()
#print("Done doing format conversion")
#print("Success! Profile contained in file: %s\n" % os.path.join(output_folder,file_base_name+".profile"))


#################################################################################################

#Clean up files
#for kmer_size in kmer_sizes:
#	if os.path.isfile(os.path.join(output_folder,file_base_name+"-"+str(kmer_size)+"mers.jf")):
#		os.remove(os.path.join(output_folder,file_base_name+"-"+str(kmer_size)+"mers.jf"))

if not save_x:
	if os.path.isfile(os.path.join(output_folder,file_base_name+"-x.txt")):
		os.remove(os.path.join(output_folder,file_base_name+"-x.txt"))






