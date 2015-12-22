#~/bin/sh
outDir="out"
nProcs=15

python generate_taxonomy.py --output ${outDir}/taxonomy.txt --output_red ${outDir}/taxonomy_reduced.txt --pickle ${outDir}/taxonomy.pkl | tee ${outDir}/generate_taxonomy.txt
python repophlan_get_microbes.py --taxonomy ${outDir}/taxonomy_reduced.txt --out_dir ${outDir}/microbes --nproc ${nProcs} --out_summary ${outDir}/repophlan_microbes.txt | tee ${outDir}/repophlan_microbes.log
python generate_taxonomy_taxid.py --output ${outDir}/taxonomy_taxID.txt --output_red ${outDir}/taxonomy_reduced_taxID.txt --pickle ${outDir}/taxonomy_taxID.pkl | tee ${outDir}/generate_taxonomy_tax_ID.txt
python repophlan_get_viruses.py --taxonomy ${outDir}/taxonomy_reduced_taxID.txt --out_dir ${outDir}/viruses --out_summary ${outDir}/repophlan_viruses.txt | tee ${outDir}/repophlan_viruses.log