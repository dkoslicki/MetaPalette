#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import urllib2
import sys
import argparse
from StringIO import StringIO
import gzip
import os
import logging

def read_params(args):
    parser = argparse.ArgumentParser(description='')
    arg = parser.add_argument
    arg( '--refseq_virus_gbff', type=str,
         default = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz",
         help="The gbff NCBI file for viruses (default is ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz)")
    arg( '--taxonomy', metavar='taxonomy', required = True, type=str,
         help="The taxonomy file")
    arg( '--out_dir', required = True, default = None, type = str, 
         help="The output folder")
    arg( '--out_summary', required = True, default = None, type = str, 
         help="The output summary file")
    return vars(parser.parse_args())


if __name__ == '__main__':
    par = read_params(sys.argv)
    logging.basicConfig(level=logging.INFO, stream=sys.stdout, 
                        format = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logger = logging.getLogger(sys.argv[0])
    logger.info('Reading the taxonomy from '+par['taxonomy']+'... ')
    with open(par['taxonomy']) as inpf:
        taxids2taxonomy = dict([l.strip().split('\t')[1:] for l in inpf])
    logger.info('Done.')
    
    
    logger.info('Downloading and reading '+par['refseq_virus_gbff']+'... ')
    table = urllib2.urlopen( par['refseq_virus_gbff'] )
    f = gzip.GzipFile(fileobj= StringIO( table.read())  )
    logger.info('Done.')



    if not os.path.exists( par['out_dir'] ):
        os.mkdir( par['out_dir'] )
    outdir = par['out_dir']+ "/"

    summary = {}
    viruses = {}

    nrec = 0
    for seq_record in SeqIO.parse( f, "genbank"):
        accession = seq_record.annotations['accessions'][0]
        logger.info('Processing accession '+accession)

        ffn, faa = [], []
        taxid = None

        bioproject = None
        for dbxref in seq_record.dbxrefs:
            if 'BioProject' in dbxref:
                bioproject = dbxref.split("BioProject:")[-1]
        if bioproject is None:
            logger.warning('Accession '+accession+' without Bioproject, skipping it')
            continue

        if bioproject not in viruses:
            viruses[bioproject] = {'ffn':[],'faa':[],'fna':[],'taxon':""}
        
        viruses[bioproject]['fna'].append( seq_record )

        for feat in seq_record.features:
            if feat.type == 'source':
                if 'db_xref' in feat.qualifiers:
                    for tn in feat.qualifiers['db_xref']:
                        if 'taxon' in tn:
                            viruses[bioproject]['taxon'] = tn.split(":")[-1]
    
            if feat.type == 'gene':
                if 'db_xref' in feat.qualifiers:
                    gene_id = feat.qualifiers['db_xref'][0]
                    gene = feat.location.extract(seq_record)
                    gene.id = gene_id
                    viruses[bioproject]['ffn'].append( gene )
            if feat.type == 'CDS':
                if 'translation' in feat.qualifiers:
                    prot = SeqRecord( Seq(feat.qualifiers['translation'][0])  )
                    if 'protein_id' in feat.qualifiers:
                        prot.id = feat.qualifiers['protein_id'][0]
                    if 'product' in feat.qualifiers:
                        prot.description = feat.qualifiers['product'][0]
                    viruses[bioproject]['faa'].append( prot )
        nrec += 1
        logger.info('Processing '+accession+' [record # '+str(nrec)+'] ')
       
    for nrec,(k,v) in enumerate(viruses.items()):
        taxid = v['taxon']
        if taxid is None:
            logger.warning('No tax ID for '+accession.id+' skipping it')
            continue
        
        if taxid not in taxids2taxonomy:
            logger.info('TaxId not found in taxonomy:  '+taxid)
            continue
        
        SeqIO.write( v['fna'], outdir+k+".fna", "fasta" )
        SeqIO.write( v['ffn'], outdir+k+".ffn", "fasta" )
        SeqIO.write( v['faa'], outdir+k+".faa", "fasta" )
  
        summary[k] = {  'taxonomy':taxids2taxonomy[taxid],
        						'taxid': taxid,
                                'lfna': outdir+k+".fna", 
                                'lffn': outdir+k+".ffn", 
                                'lfaa': outdir+k+".faa" } 
        #print summary

        logger.info('Files exported for '+k+' [record # '+str(nrec)+'] '+",".join([outdir+k+".fna",outdir+k+".ffn",outdir+k+".faa"]))

    logger.info('Writing summary file to: '+par['out_summary'])
    with open( par['out_summary'], "w" ) as outf:
        for k,v in sorted(summary.items(),key=lambda y:y[0]):
            outf.write( "\t".join( [k,v['taxid'],v['taxonomy']] ) +"\n"  )
    logger.info('Done. Exiting.')
    

