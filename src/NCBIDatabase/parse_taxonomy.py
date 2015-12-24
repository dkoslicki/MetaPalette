import os, sys, argparse

def read_params(args):
	parser = argparse.ArgumentParser(description='')
	arg = parser.add_argument
	arg( '--taxonomy', metavar='taxonomy', required = True, type=str, help="The taxonomy file")
	arg( '--out_dir', required = True, default = None, type = str, help="The output folder")
	arg( '--repophlan_dir', required = True, default = None, type = str, help="The repophland folder. eg: /RepoPhlAn/out")
	return vars(parser.parse_args())

if __name__ == '__main__':
	par = read_params(sys.argv)
	if not os.path.exists( par['out_dir'] ):
		os.mkdir( par['out_dir'] )
	outdir = par['out_dir']+ "/"
	repophlan_dir = par['repophlan_dir']+"/"
	Accession2TaxID = dict()
	Accession2FileName = dict()
	Accession2TaxIDFile = os.path.join(repophlan_dir,"repophlan_microbes.txt")
	taxonomyFile = os.path.join(repophlan_dir,"taxonomy_reduced_taxID.txt")
	fna_dir = os.path.join(repophlan_dir,"microbes/fna/")
	fid = open(Accession2TaxIDFile,"r")
	i=0
	for line in fid:
		if i==0:
			pass
		else:
			Accession = line.strip().split('\t')[0]
			TaxID = line.strip().split('\t')[27]
			#filename = line.strip().split('\t')[12]
			#filename = os.path.join(fna_dir,accession+".fna.bz2")
			Accession2TaxID[Accession] = TaxID
		i=1
	fid.close()
	
	taxids2taxonomy = dict()
	fid = open(taxonomyFile,"r")
	for line in fid:
		taxID = line.strip().split('\t')[1]
		taxonomy = line.strip()
		taxids2taxonomy[taxID] = taxonomy
	
	accessions = Accession2TaxID.keys()
	
	bacteria_taxonomies = list()
	bacteria_filenames = list()
	
	archaea_taxonomies = list()
	archaea_filenames = list()
	
	eukaryota_taxonomies = list()
	eukaryota_filenames = list()
	
	for accession in accessions:
		if Accession2TaxID[accession] in taxids2taxonomy:
			if os.path.isfile(os.path.join(fna_dir,accession+".fna.bz2")):
				taxonomy = taxids2taxonomy[Accession2TaxID[accession]]
				kingdom = taxonomy.split('\t')[2].split('|')[0].split('_')[-1]
				if kingdom=="Bacteria":
					bacteria_taxonomies.append(taxonomy)
					bacteria_filenames.append(os.path.abspath(os.path.join(fna_dir,accession+".fna.bz2")))
				elif kingdom=="Archaea":
					archaea_taxonomies.append(taxonomy)
					archaea_filenames.append(os.path.abspath(os.path.join(fna_dir,accession+".fna.bz2")))
				elif kingdom=="Eukaryota":
					eukaryota_taxonomies.append(taxonomy)
					eukaryota_filenames.append(os.path.abspath(os.path.join(fna_dir,accession+".fna.bz2")))
	
	fid = open(os.path.join(outdir,"BacteriaTaxonomy.txt"),"w")
	for taxonomy in bacteria_taxonomies:
		fid.write("%s\n" % taxonomy)
	fid.close()
	fid = open(os.path.join(outdir,"BacteriaFileNames.txt"),"w")
	for filename in bacteria_filenames:
		fid.write("%s\n" % filename)
	fid.close()
	
	fid = open(os.path.join(outdir,"ArchaeaTaxonomy.txt"),"w")
	for taxonomy in archaea_taxonomies:
		fid.write("%s\n" % taxonomy)
	fid.close()
	fid = open(os.path.join(outdir,"ArchaeaFileNames.txt"),"w")
	for filename in archaea_filenames:
		fid.write("%s\n" % filename)
	fid.close()
	
	fid = open(os.path.join(outdir,"EukaryotaTaxonomy.txt"),"w")
	for taxonomy in eukaryota_taxonomies:
		fid.write("%s\n" % taxonomy)
	fid.close()
	fid = open(os.path.join(outdir,"EukaryotaFileNames.txt"),"w")
	for filename in eukaryota_filenames:
		fid.write("%s\n" % filename)
	fid.close()
	
	#Then do the same thing for the viruses....
	virus_taxonomy_filename = os.path.join(repophlan_dir,"repophlan_viruses.txt")
	virus_dir = os.path.join(repophlan_dir,"viruses")
	fid = open(virus_taxonomy_filename,"r")
	virus_temp_accessions = list()
	virus_accession2TaxID = dict()
	for line in fid:
		accession = line.strip().split()[0]
		taxID=line.strip().split()[1]
		virus_temp_accessions.append(accession)
		virus_accession2TaxID[accession] = taxID
	
	virus_taxonomies = list()
	virus_filenames = list()
	
	for accession in virus_temp_accessions:
		if os.path.isfile(os.path.join(virus_dir,accession+".fna")):
			if accession in virus_accession2TaxID:
				if virus_accession2TaxID[accession] in taxids2taxonomy:
					taxonomy = taxids2taxonomy[virus_accession2TaxID[accession]]
					kingdom = taxonomy.split('\t')[2].split('|')[0].split('_')[-1]
					if kingdom == "Viruses":
						virus_taxonomies.append(taxonomy)
						virus_filenames.append(os.path.abspath(os.path.join(virus_dir,accession+".fna")))
	
	fid = open(os.path.join(outdir,"VirusTaxonomy.txt"),"w")
	for taxonomy in virus_taxonomies:
		fid.write("%s\n" % taxonomy)
	fid.close()
	
	fid = open(os.path.join(outdir,"VirusFileNames.txt"),"w")
	for filename in virus_filenames:
		fid.write("%s\n" % filename)
	fid.close()
