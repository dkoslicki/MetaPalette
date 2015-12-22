#~/bin/sh
#Change the "outDir" to the output directory you wish.
#Change the nProcs to the number of parallel processes to download from NCBI (too many will cause NCBI to lock you out)
#Comment out the line beginning with "python repophlan_get_microbes.py" if you do not wish the bacteria, archaea, and eukaryota to be downloaded
#Comment out the line beginning with "python repophlan_get_viruses.py" if you do not wish the viruses to be downloaded
#Changing anything else risks an error

outDir="out"
nProcs=15

python generate_taxonomy.py --output ${outDir}/taxonomy.txt --output_red ${outDir}/taxonomy_reduced.txt --pickle ${outDir}/taxonomy.pkl | tee ${outDir}/generate_taxonomy.txt
python repophlan_get_microbes.py --taxonomy ${outDir}/taxonomy_reduced.txt --out_dir ${outDir}/microbes --nproc ${nProcs} --out_summary ${outDir}/repophlan_microbes.txt | tee ${outDir}/repophlan_microbes.log
python generate_taxonomy_taxid.py --output ${outDir}/taxonomy_taxID.txt --output_red ${outDir}/taxonomy_reduced_taxID.txt --pickle ${outDir}/taxonomy_taxID.pkl | tee ${outDir}/generate_taxonomy_tax_ID.txt
python repophlan_get_viruses.py --taxonomy ${outDir}/taxonomy_reduced_taxID.txt --out_dir ${outDir}/viruses --out_summary ${outDir}/repophlan_viruses.txt | tee ${outDir}/repophlan_viruses.log
python parse_taxonomy.py --taxonomy ${outDir}/taxonomy_reduced_taxID.txt --out_dir ${outDir} --repophlan_dir ${outDir}
sed -i 's/norank__1_root|//g' ${outDir}/taxonomy_reduced_taxID.txt
sed -i 's/_noname//g' ${outDir}/taxonomy_reduced_taxID.txt
sed -i 's/norank__1_root|//g' ${outDir}/repophlan_viruses.txt
sed -i 's/_noname//g' ${outDir}/repophlan_viruses.txt