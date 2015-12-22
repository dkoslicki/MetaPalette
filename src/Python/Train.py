#Should implement something that checks if all the counts are there, and only counts the ones that it needs to
import os, sys, shutil, subprocess, getopt
from itertools import *
from multiprocessing import Pool, freeze_support
import numpy as np
import h5py

#python Train.py -i /data/temp/Data/FileNamesFullPath10.txt -o /nfs1/Koslicki_Lab/koslickd/CommonKmers/Python/Output -b /nfs1/Koslicki_Lab/koslickd/Bcalm/bcalm/./bcalm -r /data/temp -j /home/pi/koslickd/jellyfish-2.2.3/bin/./jellyfish -c /nfs1/Koslicki_Lab/Backup/koslickd/CAMI/CommonKmers/src/CountInFile/./count_in_file -t 48 -k 30


bcalm_binary = 'bcalm'
jellyfish_binary = 'jellyfish'
count_in_file_binary ='count_in_file'
chunk_size = 200
try:
	opts, args = getopt.getopt(sys.argv[1:],"hi:o:b:r:j:c:t:k:s:",["Help=","InputListOfFiles=","OutputFolder=", "BcalmBinary=","RamdiskLocation=","JellyfishBinary=","CountInFileBinary=","Threads=","KmerCountingThreads=","ChunkSize="])
except getopt.GetoptError:
	print 'Unknown option, call using: python Train.py -i <InputListOfFiles> -o <OutputFolder> -b <BcalmBinary> -r <RamdiskLocation> -j <JellyfishBinary> -c <CountInFileBinary> -t <Threads> -k <KmerCountingThreads> -s <ChunkSize>'
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print 'python Train.py -i <InputListOfFiles> -o <OutputFolder> -b <BcalmBinary> -r <RamdiskLocation> -j <JellyfishBinary> -c <CountInFileBinary> -t <Threads> -k <KmerCountingThreads> -s <ChunkSize>'
		sys.exit(2)
	elif opt in ("-i", "--InputListOfFiles"):
		input_files = arg
	elif opt in ("-o", "--OutputFolder"):
		output_folder = arg
	elif opt in ("-b","--BcalmBinary"):
		bcalm_binary = arg
	elif opt in ("-r", "--RamdiskLocation"):
		ramdisk_location = arg
	elif opt in ("-j", "--JellyfishBinary"):
		jellyfish_binary = arg
	elif opt in ("-c", "--CountInFileBinary"):
		count_in_file_binary = arg
	elif opt in ("-t", "--Threads"):
		num_threads = int(arg)
	elif opt in ("-k", "--KmerCountingThreads"):
		kmer_counting_threads = int(arg)
	elif opt in ("-s", "--ChunkSize"):
		chunk_size = int(arg)

#These are the kmer sizes to train on
kmer_sizes = [30,50]

if chunk_size*num_threads>=1024:
	print("Warning: chunk_size("+str(chunk_size)+")*num_threads("+str(num_threads)+") is greater than or equal to 1024 (the typical value of ulimit -n, the maximum number of allowed open files). Please reduce the chunk_size (-s) or num_threads (-t) or risk the program throwing an error.")
if not os.path.isdir(output_folder):
	print("Error: Output folder " + output_folder + " does not exist.")
	sys.exit(2)
if not os.path.isfile(input_files):
	print("Error: List of input files " + input_files + " does not exist.")
	sys.exit(2)
if not os.path.isdir(ramdisk_location):
	print("Error: fast IO device location " + ramdisk_location + " does not exist.")
	sys.exit(2)

#Read in file names
fid = open(input_files,'r')
file_names = fid.readlines()
fid.close()
file_names = [name.strip() for name in file_names]

for file_name in file_names:
	if not os.path.isfile(os.path.abspath(file_name.strip())):
		print("Error: file " + file_name + " does not exist but given in input file " + input_files)
		sys.exit(2)

#Form k-mer counts in parallel
if not os.path.isdir(os.path.join(output_folder,"Counts")):
	os.makedirs(os.path.join(output_folder,"Counts"))

#Can put automatic decompression here
def count_kmers(file, kmer_size):
	extension = os.path.splitext(file)[1]
	if extension==".bz2" or extension==".bz":
		cmd = "bzip2 -c " + file + " | " + jellyfish_binary + " count /dev/fd/0 -m "+str(kmer_size)+" -t 1 -s 10M --out-counter-len 3 --disk -C -o " + os.path.join(output_folder,"Counts",os.path.basename(file)+"-"+str(kmer_size)+"mers.jf")
	elif extension==".gz" or extension==".z" or extension==".Z":
		cmd = "gunzip -c " + file + " | " + jellyfish_binary + " count /dev/fd/0 -m "+str(kmer_size)+" -t 1 -s 10M --out-counter-len 3 --disk -C -o " + os.path.join(output_folder,"Counts",os.path.basename(file)+"-"+str(kmer_size)+"mers.jf")
	else:
		cmd = jellyfish_binary + " count " + file + " -m "+str(kmer_size)+" -t 1 -s 10M --out-counter-len 3 --disk -C -o " + os.path.join(output_folder,"Counts",os.path.basename(file)+"-"+str(kmer_size)+"mers.jf")
	if not os.path.isfile(os.path.join(output_folder,"Counts",os.path.basename(file)+"-"+str(kmer_size)+"mers.jf")):
		test = subprocess.check_call(cmd, shell = True)

def count_kmers_star(arg):
	return count_kmers(*arg)

pool = Pool(processes = kmer_counting_threads)

#Count k-mers
for kmer_size in kmer_sizes:
	res = pool.map(count_kmers_star, izip(file_names, repeat(kmer_size)));

#Now to form the common kmer matrix
def count_in_file(file_list, kmer_size):
	#cmd = count_in_file_binary + " " + os.path.join(output_folder,"Counts",os.path.basename(file_tuple[0]))+"-"+str(kmer_size)+"mers.jf" + " " + os.path.join(output_folder,"Counts",os.path.basename(file_tuple[1]))+"-"+str(kmer_size)+"mers.jf"
	cmd = count_in_file_binary + " " + " ".join(file_list)
	out = subprocess.check_output(cmd, shell = True)
	res = np.fromstring(out, sep=" ").reshape((len(file_list), len(file_list)))
	return res

def count_in_file_star(arg):
	return count_in_file(*arg)

#Make the CKM's but chunk the indices into sublocks so we don't get a bunch of thrashing/memory issues
num_files = len(file_names)
for kmer_size in kmer_sizes:
	count_file_names = list()
	to_count_file_names = list()
	to_count_file_names_lengths = list()
	ijs = list()
	for file_name in file_names:
		count_file_names.append(os.path.join(output_folder,"Counts",os.path.basename(file_name)+"-"+str(kmer_size)+"mers.jf"))
	ckm = np.zeros((num_files,num_files),dtype=np.int64)
	for i in range(0,num_files+chunk_size,chunk_size):
		for j in range(0,num_files+chunk_size,chunk_size):
			if j>i:
				icount_file_names = list()
				jcount_file_names = list()
				for ii in range(i,i+chunk_size):
					if ii<num_files:
						icount_file_names.append(count_file_names[ii])
				for jj in range(j,j+chunk_size):
					if jj<num_files:
						jcount_file_names.append(count_file_names[jj])
				if len(icount_file_names)>0 and len(jcount_file_names)>0:
					to_count_file_names.append(icount_file_names + jcount_file_names)
					to_count_file_names_lengths.append((len(icount_file_names),len(jcount_file_names)))
					ijs.append((i,j))
	pool.close()
	pool = Pool(processes = num_threads)
	res = pool.map(count_in_file_star, izip(to_count_file_names, repeat(kmer_size)));
	#Turn the result into the Common Kmer Matrix. Put square submatrices of mat into proper place in ckm (draw a picture to see what's going on)
	for i in range(len(res)):
		mat = res[i]
		if i==0:
			ckm[0:mat.shape[0], 0:mat.shape[0]] = mat
		else:
			#(1,1)
			ckm[ijs[i][0]:(ijs[i][0]+to_count_file_names_lengths[i][0]), ijs[i][0]:(ijs[i][0]+to_count_file_names_lengths[i][0])] = mat[0:to_count_file_names_lengths[i][0], 0:to_count_file_names_lengths[i][0]]
			#(1,2)
			ckm[ijs[i][0]:(ijs[i][0]+to_count_file_names_lengths[i][0]), ijs[i][1]:(ijs[i][1]+to_count_file_names_lengths[i][1])] = mat[0:to_count_file_names_lengths[i][0], to_count_file_names_lengths[i][0]:(to_count_file_names_lengths[i][1]+to_count_file_names_lengths[i][0])]
			#(2,1)
			ckm[ijs[i][1]:(ijs[i][1]+to_count_file_names_lengths[i][1]), ijs[i][0]:(ijs[i][0]+to_count_file_names_lengths[i][0])] = mat[to_count_file_names_lengths[i][0]:(to_count_file_names_lengths[i][1]+to_count_file_names_lengths[i][0]), 0:to_count_file_names_lengths[i][0]]
			#(2,2)
			ckm[ijs[i][1]:(ijs[i][1]+to_count_file_names_lengths[i][1]), ijs[i][1]:(ijs[i][1]+to_count_file_names_lengths[i][1])] = mat[to_count_file_names_lengths[i][0]:(to_count_file_names_lengths[i][1]+to_count_file_names_lengths[i][0]), to_count_file_names_lengths[i][0]:(to_count_file_names_lengths[i][1]+to_count_file_names_lengths[i][0])]
	#Save the ckm matrices
	fid = h5py.File(os.path.join(output_folder,"CommonKmerMatrix-"+str(kmer_size)+"mers.h5"),'w')
	dset = fid.create_dataset("common_kmers", data=ckm)
	fid.close()

#Form the bcalms
#Single file bcalm function
counts_folder = os.path.join(output_folder,"Counts")
def form_bcalms(input_file, output_folder_bcalm, bcalm_binary, ramdisk_location, jellyfish_binary, counts_directory):
	if not os.path.isfile(os.path.join(output_folder_bcalm,input_file+"-30mers.bcalm.fa")):
		input_file = os.path.basename(input_file)
		#Make temporary directory
		if os.path.exists(os.path.join(ramdisk_location,input_file)):
			shutil.rmtree(os.path.join(ramdisk_location,input_file))
		os.makedirs(os.path.join(ramdisk_location,input_file))
		#Make .dot file
		cmd = jellyfish_binary +" dump "+os.path.join(counts_folder,input_file+"-30mers.jf")+" -c -t | cut -f 1 | tr '[:upper:]' '[:lower:]' | sed 's/$/;/g' > " + os.path.join(ramdisk_location,input_file,input_file+"-30mers.dot")
		temp = subprocess.check_call(cmd, shell = True)
		#Run Bcalm
		working_dir = os.getcwd()
		os.chdir(os.path.join(ramdisk_location,input_file))
		cmd = bcalm_binary +" "+ input_file +"-30mers.dot " + input_file + "-30mers.bcalm 5"
		FNULL = open(os.devnull, 'w')
		temp = subprocess.check_call(cmd, shell = True, stdout=FNULL)
		FNULL.close()
		#Move the file and delete the temp directory
		os.chdir(working_dir)
		cmd = "cat " + os.path.join(ramdisk_location, input_file, input_file+"-30mers.bcalm") + " | sed 's/;//g' | tr '[:lower:]' '[:upper:]'| sed 's/^/>seq\\n/g' > " + os.path.join(output_folder_bcalm,input_file+"-30mers.bcalm.fa")
		temp = subprocess.check_call(cmd, shell = True)
		shutil.rmtree(os.path.join(ramdisk_location, input_file))

def form_bcalms_star(arg):
	return form_bcalms(*arg)

#compute bcalms in parallel
if not os.path.isdir(os.path.join(output_folder,"Bcalms")):
	os.makedirs(os.path.join(output_folder,"Bcalms"))
pool = Pool(processes = kmer_counting_threads)
res = pool.map(form_bcalms_star, izip(file_names, repeat(os.path.join(output_folder,"Bcalms")), repeat(bcalm_binary), repeat(ramdisk_location), repeat(jellyfish_binary), repeat(os.path.join(output_folder,"Counts"))))
pool.close()

#Remove counts
#shutil.rmtree(counts_folder)

#Make FileNames.txt file
fid = open(os.path.join(output_folder,"FileNames.txt"),'w')
for file in file_names:
	fid.write(os.path.abspath(file)+"\n")
fid.close()


