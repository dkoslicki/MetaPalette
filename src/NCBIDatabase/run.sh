#~/bin/sh
mkdir -p out
t=`date "+%d%m%Y"`

/local/cluster/bin/python generate_taxonomy.py --output out/taxonomy.txt --output_red out/taxonomy_reduced.txt --pickle out/taxonomy.pkl | tee out/generate_taxonomy.txt
/local/cluster/bin/python repophlan_get_microbes.py --taxonomy out/taxonomy_reduced.txt --out_dir out/microbes --nproc 15 --out_summary out/repophlan_microbes.txt | tee out/repophlan_microbes.log
/local/cluster/bin/python generate_taxonomy_taxid.py --output out/taxonomy_taxID.txt --output_red out/taxonomy_reduced_taxID.txt --pickle out/taxonomy_taxID.pkl | tee out/generate_taxonomy_tax_ID.txt

#I think this is how to get the correct taxonomy for the viruses (note the taxID taxonomy)
/local/cluster/bin/python repophlan_get_viruses.py --taxonomy out/taxonomy_reduced_taxID.txt --out_dir out/viruses --out_summary out/repophlan_viruses.txt | tee out/repophlan_viruses.log