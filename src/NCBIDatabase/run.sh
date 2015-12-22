#~/bin/sh
mkdir -p out
t=`date "+%d%m%Y"`

/local/cluster/bin/python generate_taxonomy.py --output out/taxonomy_${t}.txt --output_red out/taxonomy_reduced_${t}.txt --pickle out/taxonomy_${t}.pkl | tee out/generate_taxonomy_${t}.txt
/local/cluster/bin/python repophlan_get_microbes.py --taxonomy out/taxonomy_reduced_${t}.txt --out_dir out/microbes_${t} --nproc 15 --out_summary out/repophlan_microbes_${t}.txt | tee out/repophlan_microbes_${t}.log
#/local/cluster/bin/python repophlan_get_viruses.py --taxonomy out/taxonomy_reduced_${t}.txt --out_dir out/viruses_${t} --out_summary out/repophlan_viruses_${t}.txt | tee out/repophlan_viruses_${t}.log
/local/cluster/bin/python generate_taxonomy_taxid.py --output out/taxonomy_taxID_${t}.txt --output_red out/taxonomy_reduced_taxID_${t}.txt --pickle out/taxonomy_taxID_${t}.pkl | tee out/generate_taxonomy_tax_ID${t}.txt

#I think this is how to get the correct taxonomy for the viruses (note the taxID taxonomy)
/local/cluster/bin/python repophlan_get_viruses.py --taxonomy out/taxonomy_reduced_taxID_${t}.txt --out_dir out/viruses_${t} --out_summary out/repophlan_viruses_${t}.txt | tee out/repophlan_viruses_${t}.log