#This is the core plotting package to be called from Plot.py
import numpy as np
import h5py
from Bio import Phylo
from Bio.Phylo import NewickIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from ete2 import Tree, faces, TreeStyle, COLOR_SCHEMES, TextFace, BarChartFace, CircleFace, AttrFace, NodeStyle, RectFace
import math
import os

def MakePlot(x, org_names, ckm30, ckm50, outgroup, outfile):
	#x_file = "../DataForNJTreeWithTrue/test1-reads.fa-x.txt"
	#grinder_x_file = "../DataForNJTreeWithTrue/test1_grinder-x.txt"
	#file_name_file = "../DataForNJTreeWithTrue/AllTrainedData/FileNames.txt"
	#x_file_names_file = "../DataForNJTreeWithTrue/FullFileNames.txt"
	#ckm30_file = "../DataForNJTreeWithTrue/AllTrainedData/CommonKmerMatrix-30mers.h5"
	#ckm50_file = "../DataForNJTreeWithTrue/AllTrainedData/CommonKmerMatrix-50mers.h5"
	#ckm30_train_file = "../DataForNJTreeWithTrue/CommonKmerMatrix-30mers.h5"
	#ckm50_train_file = "../DataForNJTreeWithTrue/CommonKmerMatrix-50mers.h5"
	
	##Read in the x vector
	#fid = open(x_file,'r')
	#x = map(lambda y: float(y),fid.readlines())
	#fid.close()

	#Normalize the x vector
	x = map(lambda y: y/sum(x),x)

	##Read in the true x vector
	#fid = open(grinder_x_file,'r')
	#x_true = map(lambda y: float(y),fid.readlines())
	#fid.close()

	##Read in the basis for the ckm matrices
	#file_names = list()
	#fid = open(file_name_file,'r')
	#for line in fid:
	#	#file_names.append(os.path.basename(line.strip()))
	#	file_names.append(os.path.basename(line.strip()).split('.')[0])
	#fid.close()

	#Read in the taxonomy
	#taxonomy = list()
	#fid = open('CommonKmerTrainingOnAllTestData/Testtaxonomy.txt','r')
	#for line in fid:
	#    taxonomy.append('_'.join(line.split()[0].split("_")[1:]))
	
	#fid.close()

	#Taxonomy isn't needed for these figures, so just make the taxonomy the file names
	#taxonomy = file_names

	#Read in the basis for the reconstructed vector
	#x_file_names = list()
	#fid = open(x_file_names_file,'r')
	#for line in fid:
	#	#x_file_names.append(os.path.basename(line.strip()))
	#	x_file_names.append(os.path.basename(line.strip()).split('.')[0])
	#fid.close()


	#Read in the common kmer matrix
	#f=h5py.File(ckm30_file,'r')
	#ckm30=np.array(f['common_kmers'],dtype=np.float64)
	#f.close()
	#f=h5py.File(ckm50_file,'r')
	#ckm50=np.array(f['common_kmers'],dtype=np.float64)
	#f.close()
	ckm30_norm = np.multiply(ckm30,1/np.diag(ckm30))
	ckm50_norm = np.multiply(ckm50,1/np.diag(ckm50))
	num_rows = ckm30_norm.shape[0]
	num_cols = ckm30_norm.shape[1]
	names = org_names
	matrix=list()
	for i in range(num_rows):
		matrix.append([.5*(1-.5*ckm30_norm[i,j]-.5*ckm30_norm[j,i])+.5*(1-.5*ckm50_norm[i,j]-.5*ckm50_norm[j,i]) for j in range(i+1)])

	#Read in the training common kmer matrices
	#f=h5py.File(ckm30_train_file,'r')
	#ckm30_train=np.array(f['common_kmers'],dtype=np.float64)
	#f.close()
	#f=h5py.File(ckm50_train_file,'r')
	#ckm50_train=np.array(f['common_kmers'],dtype=np.float64)
	#f.close()
	#ckm30_train_norm = np.multiply(ckm30_train,1/np.diag(ckm30_train))
	#ckm50_train_norm = np.multiply(ckm50_train,1/np.diag(ckm50_train))
	ckm_ave_train = .5*ckm30_norm+.5*ckm50_norm
	ckm_ave_train_dist = dict()
	for i in range(len(org_names)):
		ckm_ave_train_dist[org_names[i]] = [.5*ckm_ave_train[i,j]+.5*ckm_ave_train[j,i] for j in range(len(org_names))]


	#Construct the tree. Note I could use RapidNJ here, but a few tests have shown that the trees that RapidNJ creates are rubbish.
	dm = _DistanceMatrix(names, matrix)
	constructor = DistanceTreeConstructor()
	tree = constructor.nj(dm)
	t=Tree(tree.format('newick'),format=1)
	#tree.format('newick')
	#Phylo.draw_ascii(tree)


	#Now I will put internal nodes in a certain phylogenetic distance between the root and a given node.
	#Function to insert a node at a given distance
	def insert_node(t, name_to_insert, insert_above, dist_along):
		insert_at_node = t.search_nodes(name=insert_above)[0]
		parent = (t&insert_above).up
		orig_branch_length = t.get_distance(insert_at_node,parent)
		if orig_branch_length < dist_along:
			raise ValueError("error: dist_along larger than orig_branch_length")
		removed_node = insert_at_node.detach()
		removed_node.dist = orig_branch_length - dist_along
		added_node = parent.add_child(name=name_to_insert, dist=dist_along)
		added_node.add_child(removed_node)

	#Function to insert a node some % along a branch, taking into account the ckm distances and nodes already created in the NJ tree (and what distance their descendants are from everyone else)
	def insert_hyp_node(t, leaf_name, percent, ckm_ave_train_dist, org_names):
		dists = map(lambda y: abs(y-percent), ckm_ave_train_dist[leaf_name])
		nearby_indicies = list()
	#	for i in range(len(dists)):
	#		if dists[i]<=.05:
	#			nearby_indicies.append(i)
		nearby_names = list()
		#If there are no nearby indicies, add the closest organism to the given percent
		if nearby_indicies==[]:
			nearby_names.append(org_names[dists.index(min(dists))])
		else:
			for i in range(len(nearby_indicies)):
				nearby_names.append(org_names[i])
		mean_dist = np.mean(map(lambda y: ckm_ave_train_dist[leaf_name][org_names.index(y)],nearby_names))
		nearby_names.append(leaf_name)
		LCA = t.get_common_ancestor(nearby_names)
		LCA_to_leaf_dist = t.get_distance(LCA,leaf_name)
		#divide the dist to the right/left of the LCA node by the number of percentage points in there
		if LCA.name==t.name:
			percent_dist = percent*LCA_to_leaf_dist
			if mean_dist <= percent:
				child_node = (t&leaf_name)
			else:
				child_node = (t&nearby_names[0])#This means "go up from root" in the direction of the nearest guy
			ancestor_node = (t&child_node.name).up
		elif mean_dist <= percent:
			#percent_dist = LCA_to_leaf_dist + abs(percent-mean_dist)*(t.get_distance(leaf_name)-LCA_to_leaf_dist)/(1-mean_dist)
			percent_dist = t.get_distance(LCA) + abs(percent-mean_dist)*(LCA_to_leaf_dist)/(1-mean_dist)
			child_node = (t&leaf_name)
			ancestor_node = (t&child_node.name).up
		else:
			#percent_dist = LCA_to_leaf_dist - abs(percent-mean_dist)*(LCA_to_leaf_dist)/(mean_dist)
			percent_dist = t.get_distance(LCA) - abs(percent-mean_dist)*(t.get_distance(LCA))/(mean_dist)
			child_node = (t&leaf_name)
			ancestor_node = (t&child_node.name).up
		while t.get_distance(t.name, ancestor_node) > percent_dist:
			child_node = ancestor_node
			ancestor_node = (t&child_node.name).up
		insert_node(t, leaf_name+"_"+str(percent), child_node.name, percent_dist-t.get_distance(t.name, ancestor_node))


	#set Bongori as outgroup
	#ancestor = t.get_common_ancestor("bongori_IP780492","bongori_RKI1373","bongori_CDC270376","bongori_RKI1786","bongori_IP477084","bongori_CEIM46049","bongori_CEIM24450","bongori_RKI1398","bongori_12419_ATCC_43975")
	#t.set_outgroup(ancestor)
	#t.set_outgroup(t&'Halobacterium_sp_DL1')
	t.set_outgroup(t&outgroup) #I will need to check that this outgroup is actually one of the names...

	#Insert hypothetical nodes
	hyp_node_names = dict()
	cutoffs = [.9,.8,.7,.6,.5,.4,.3,.2,.1]
	#cutoffs = [.9,.75,.6,.45,.3,.15]
	#cutoffs = map(lambda y: y**1.5,cutoffs)
	cutoffs = [-.5141*(val**3)+1.0932*(val**2)+0.3824*val for val in cutoffs]
	for i in range(len(org_names)):
		xi = x[i:len(x):len(org_names)]
		for j in range(1,len(cutoffs)+1):
			if xi[j]>0:
				insert_hyp_node(t, org_names[i], cutoffs[j-1],ckm_ave_train_dist, org_names)
				hyp_node_names[org_names[i]+"_"+str(cutoffs[j-1])] = [org_names[i], cutoffs[j-1], j-1] #in case there are "_" in the file names


	size_factor=250
	font_size=55

	#Now put the bubbles on the nodes
	def layout(node):
		node_style = NodeStyle()
		node_style["hz_line_width"] = 10
		node_style["vt_line_width"] = 10
		node.set_style(node_style)
		#print(node)
		if node.is_leaf():
			if node.name in org_names:
				#make reconstructed bubble
				size = x[org_names.index(node.name)]
				F = CircleFace(radius=size_factor*math.sqrt(size), color="RoyalBlue", style="sphere")
				F.border.width = None
				F.opacity = 0.6
				faces.add_face_to_node(F,node, 0, position="branch-right")
				#Denote that this was a training organism
				nameFace = AttrFace("name", fsize=font_size, fgcolor='black')
				faces.add_face_to_node(nameFace, node, 0, position="branch-right")
		elif node.name in hyp_node_names: #Otherwise it's a hypothetical node, just use recon x
			node_base_name = hyp_node_names[node.name][0]
			percent = hyp_node_names[node.name][1]
			if node_base_name in org_names:
				idx = hyp_node_names[node.name][2]
				size = x[org_names.index(node_base_name)+(idx+1)*len(x_file_names)]
				F = CircleFace(radius=size_factor*math.sqrt(size), color="RoyalBlue", style="sphere")
				F.border.width = None
				F.opacity = 0.6
				faces.add_face_to_node(F,node, 0, position="branch-top")
				#print node
				#print size
				#nameFace = AttrFace("name", fsize=font_size, fgcolor='black')
				#faces.add_face_to_node(nameFace, node, 0, position="branch-right")
			else:
				size=0
		else:
			size=0
		#print(size)
	
	ts = TreeStyle()
	ts.layout_fn = layout
	ts.mode = "r"
	#ts.mode = "c"
	ts.scale = 2*1000
	ts.show_leaf_name = False
	ts.min_leaf_separation = 50
	F = CircleFace(radius=.87*size_factor, color="RoyalBlue", style="sphere")
	F.border.width = None
	F.opacity = 0.6
	ts.legend.add_face(F,0)
	ts.legend.add_face(TextFace("  Inferred relative abundance",fsize=1.5*font_size,fgcolor="Blue"),1)
	ts.legend_position=4
	#t.show(tree_style=ts)
	t.render(outfile, w=550, units="mm", tree_style=ts)

