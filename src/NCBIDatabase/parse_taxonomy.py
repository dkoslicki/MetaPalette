import os, sys
Accession2TaxID = dict()
Accession2FileName = dict()
Accession2TaxIDFile = "/nfs1/Koslicki_Lab/koslickd/RepoPhlAn-12-20-14/out/repophlan_microbes_21122014.txt"
taxonomyFile = "/nfs1/Koslicki_Lab/koslickd/RepoPhlAn-12-20-14/out/taxonomy_reduced_taxID_22122014.txt"
fna_dir = "/nfs1/Koslicki_Lab/koslickd/RepoPhlAn-12-20-14/out/microbes_21122014/fna/"
fid = open(Accession2TaxIDFile,"r")
i=0
for line in fid:
	if i==0:
		pass
	else:
		Accession = line.strip().split('\t')[0]
		TaxID = line.strip().split('\t')[26]
		filename = line.strip().split('\t')[12]
		Accession2TaxID[Accession] = TaxID
	i=1

fid.close()

#Might want to use TaxID taxonomy
#Remove noname and norank
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
			

fid = open("BacteriaTaxonomy.txt","w")
for taxonomy in bacteria_taxonomies:
	fid.write("%s\n" % taxonomy)

fid.close()
fid = open("BacteriaFileNames.txt","w")
for filename in bacteria_filenames:
	fid.write("%s\n" % filename)

fid.close()

fid = open("ArchaeaTaxonomy.txt","w")
for taxonomy in archaea_taxonomies:
	fid.write("%s\n" % taxonomy)

fid.close()
fid = open("ArchaeaFileNames.txt","w")
for filename in archaea_filenames:
	fid.write("%s\n" % filename)

fid.close()

fid = open("EukaryotaTaxonomy.txt","w")
for taxonomy in eukaryota_taxonomies:
	fid.write("%s\n" % taxonomy)

fid.close()
fid = open("EukaryotaFileNames.txt","w")
for filename in eukaryota_filenames:
	fid.write("%s\n" % filename)

fid.close()

#Then do the same thing for the viruses....
# read in the file repophlan_viruses_21122014.txt (this will be the proper taxonomy after the modified run.sh)
# Then I just have to figure out which files correspond to which entries...
virus_taxonomy_filename = "/nfs1/Koslicki_Lab/Temp/RepoPhlAnTemp/Viruses/out/repophlan_viruses.txt"
virus_dir = "/nfs1/Koslicki_Lab/Temp/RepoPhlAnTemp/Viruses/out/viruses/"
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

fid = open("VirusTaxonomy.txt","w")
for taxonomy in virus_taxonomies:
	fid.write("%s\n" % taxonomy)

fid.close()

fid = open("VirusFileNames.txt","w")
for filename in virus_filenames:
	fid.write("%s\n" % filename)

fid.close()
