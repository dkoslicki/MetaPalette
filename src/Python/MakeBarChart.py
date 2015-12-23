#This script will make a bar chart from a CAMI profile
import getopt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys
import numpy as np
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})


try:
	opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["Help=","InputCAMIProfile=","OutputFolder="])
except getopt.GetoptError:
	print 'Unknown option, call using: python MakeBarChart.py -i <InputCAMIProfile> -o <OutputFolder>'
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print 'python MakeBarChart.py -i <InputCAMIProfile> -o <OutputFolder>'
		sys.exit(2)
	elif opt in ("-i", "--InputCAMIProfile"):
		input_file = arg
	elif opt in ("-o", "--OutputFolder"):
		output_folder = arg

#input_file="all.fq-QC-default.profile"
profile_dict = dict()
fid = open(input_file,"r")
for line in fid:
	if line[0]!="@" and line[0]!="\n" and line[0]!="" and line[0]!="#":
		rank = line.strip().split('\t')[1]
		tax_path = line.strip().split('\t')[-2]
		tax_path_split = tax_path.split('|')
		for name in tax_path_split[::-1]:
			if name=='':
				pass
			else:
				organism = name
				break
		abundance = float(line.strip().split('\t')[-1])
		if rank not in profile_dict:
			profile_dict[rank]=dict()
			profile_dict[rank][organism] = abundance
		else:
			profile_dict[rank][organism] = abundance

fid.close()

ranks = profile_dict.keys()
for rank in ranks:
	names = profile_dict[rank].keys()
	vals = [profile_dict[rank][name] for name in names]
	#If you want it to be sorted in numerical (descending order), use:
	#names = [y for (x,y) in sorted(zip(vals,names),reverse=True)]
	#vals.sort(reverse=True)
	#If you want it to be sorted in alphabetical order of names, use:
	vals = [y for (x,y) in sorted(zip(names,vals))]
	names.sort()
	if len(names)>14:
		matplotlib.rc('xtick', labelsize=8)
	else:
		matplotlib.rc('xtick', labelsize=14)
	#If it's a strain, abbreviate the names to fit them on the whole figure
	if rank=="strain":
		for index in range(len(names)):
			name = names[index]
			name_split = name.split('_')
			name_split[0]=name_split[0][0]
			name = "_".join(name_split)
			names[index] = name
	xaxis = range(len(names))
	_=plt.figure();
	_=plt.xlim([-1,len(names)]);
	_=plt.bar(xaxis, vals, align='center');
	_=plt.xticks(xaxis, names, rotation='vertical');
	_=plt.xlabel('Taxa');
	_=plt.ylabel('Rel. abundance');
	_=plt.title(rank);
	_=plt.savefig(os.path.join(output_folder,os.path.basename(input_file)+"_"+rank+".png");
	_=plt.close();


