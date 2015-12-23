#!/usr/bin/env python
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import urllib2
import argparse
import sys
from ftplib import FTP
from time import time, sleep
import pickle
import subprocess as sb
import multiprocessing as mp
from Bio import SeqIO
import math
import logging
import os
import tarfile
import StringIO
import itertools
import gzip
import bz2
import random

FTP_prot = "ftp://"
NCBI_ftp = "ftp.ncbi.nlm.nih.gov"
NCBI_assembly_folder = '/genomes/ASSEMBLY_BACTERIA/'
NCBI_assrep_folder = '/genomes/ASSEMBLY_REPORTS/Bacteria'
NCBI_prokaryotes_file = '/genomes/GENOME_REPORTS/prokaryotes.txt'
NCBI_ASREFSEQ_file = '/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt'
NCBI_ASGENBANK_file = '/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt'
GENBANK_draft_bacteria = '/genbank/genomes/Bacteria_DRAFT/'
NCBI_assemblies_all = '/genomes/all/'

logging.basicConfig(level=logging.INFO, stream=sys.stdout, 
                    format = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
logger = logging.getLogger(sys.argv[0])



toescl = ["p__Brachiopoda","p__Cnidaria","p__Annelida","p__Arthropoda","p__Chordata","p__Chytridiomycota","p__Echinodermata","p__Hemichordata","p__Mollusca","p__Nematoda","p__Platyhelminthes","p__Porifera","p__Streptophyta","p__Eukaryota_noname"]
        

def read_params(args):
    parser = argparse.ArgumentParser(description='')
    arg = parser.add_argument
    arg( '--nproc', default = 15, type = int, 
         help="The number of parallel download processes")
    arg( '--out_dir', required = True, type = str, 
         help="The output folder")
    arg( '--taxonomy', metavar='taxonomy', required = True, type=str,
         help="The taxonomy file")
    arg( '--out_summary', required = True, type = str, 
         help="The output summary file")
    return vars(parser.parse_args())

def add_protocol( fn ):
    if fn.startswith("ftp"):
        return "ftp://"+fn
    return "http://"+fn


def retry(tries, delay=3, backoff=2):
    if backoff <= 1:
        raise ValueError("backoff must be greater than 1")
    
    tries = math.floor(tries)
    if tries < 0:
        raise ValueError("tries must be 0 or greater")
    
    if delay <= 0:
        raise ValueError("delay must be greater than 0")
    
    def deco_retry(f):
        def f_retry(*args, **kwargs):
            mtries, mdelay = tries, delay # make mutable
            #rv = f(*args, **kwargs) # first attempt
            while mtries > 1:
                try:
                    return f(*args, **kwargs) 
                except Exception, e:
                    if "No such file or directory" in e.reason and "550" in e.reason:
                        msg = "No remote file found (some ffn and faa are know to be missing remotely). Aborting the download of %s. %s" % (str(args[0]), str(e) )
                        logger.warning(msg)
                        raise e 
                    else:
                        msg = "%s: %s, Retrying in %d seconds..." % (str(args[0]), str(e), mdelay)
                        logger.warning(msg)
                sleep( mdelay )
                mtries -= 1  # consume an attempt
                mdelay *= backoff  # make future wait longer
            return f(*args, **kwargs) # Ran out of tries :-(
        return f_retry # true decorator -> decorated function
    return deco_retry  # @retry(arg[, ...]) -> true decorator

@retry(tries=8, delay=20, backoff=2)
def get_remote_file_wr( url, seqtype = 'fasta' ):
    loc = urllib2.urlopen( url, timeout = 200 )

    logger.info('Parsing of '+url)

    if url.endswith( ".tgz" ):
        compressedFile = StringIO.StringIO(loc.read()) 
        tarf = tarfile.open(fileobj=compressedFile)
        seqs = []
        for m in tarf.getmembers():
            #seqs += list(StringIO.StringIO(tarf.extractfile( m )))
            seqs += list(SeqIO.parse(tarf.extractfile( m ),seqtype))
    elif url.endswith( ".gz" ):
        compressedFile = StringIO.StringIO(loc.read())
        gzipf = gzip.GzipFile(fileobj=compressedFile)
        seqs = list(SeqIO.parse(gzipf,seqtype))
    else:
        seqs = StringIO.StringIO(loc.read())
    return seqs


def is_plasmid( ppt ):
    if 'plasmid' in ppt[0]:
        return True
    return False

def is_contig( ppt ):
    if 'plasmid' in ppt[0]:
        return False
    if 'chromosome' in ppt[0]:
        return True
    if 'complete genome' in ppt[0]:
        return True
    return False

def dwl_summarize( assemblies ):
    for k,v in assemblies.items():
        allok, codingOK = True, True
        if not os.path.exists(v['fna_lname']):
            assemblies[k]['fna_lname'] = "" 
            allok, codingOK = False, False
        if not os.path.exists(v['faa_lname']):
            assemblies[k]['faa_lname'] = ""
            allok, codingOK = False, False
        if not os.path.exists(v['ffn_lname']):
            assemblies[k]['ffn_lname'] = ""
            allok, codingOK = False, False
        if not os.path.exists(v['frn_lname']):
            assemblies[k]['frn_lname'] = ""
            allok = False
        assemblies[k]['all_data'] = 'Y' if allok else 'N'
        assemblies[k]['all_coding_data'] = 'Y' if codingOK else 'N'
    return assemblies

def refseq_dwl( info ):
    try:
        fnagz = info['dwlf']+"/"+info['assembly_accession']+"_"+info['asm_name'].replace(" ","_")+"_genomic.fna.gz"
        ret = []
        try:
            ret = get_remote_file_wr( fnagz )
        except Exception, e:
            logger.error('Error in downloading of '+fnagz+' '+str(e))
        out_fna = info['fna_lname']
        fna_dir = os.path.dirname(out_fna)
        if not os.path.exists(fna_dir):
            os.mkdir( fna_dir )
       
        if ret:
            fo = bz2.BZ2File(out_fna,"w")
            SeqIO.write( ret, fo, "fasta")
            fo.close()
            fnad = SeqIO.to_dict(ret)
        
        faagz = info['dwlf']+"/"+info['assembly_accession']+"_"+info['asm_name'].replace(" ","_")+"_protein.faa.gz"
        try:
            ret = get_remote_file_wr( faagz )
        except Exception, e:
            logger.error('Error in downloading of '+faagz+' '+str(e))
       
        if not ret:
            logger.error('Error in downloading of '+faagz+': empty fna file!')
            return 

        out_faa = info['faa_lname']
        faa_dir = os.path.dirname(out_faa)
        if not os.path.exists(faa_dir):
            os.mkdir( faa_dir )
        fo = bz2.BZ2File(out_faa,"w")
        SeqIO.write( ret, fo, "fasta")
        fo.close()
        
        gbffgz = info['dwlf']+"/"+info['assembly_accession']+"_"+info['asm_name'].replace(" ","_")+"_genomic.gbff.gz"
        ret, ffn, frn = [], [], []
        try:
            ret = get_remote_file_wr( gbffgz, seqtype = 'genbank' )
        except Exception, e:
            logger.error('Error in downloading of '+gbffgz+' '+str(e))
        for s in ret:
            for f in s.features:
                if f.type == 'CDS':
                    gene = f.location.extract(fnad[s.id])
                    if 'protein_id' in f.qualifiers and 'locus_tag' in f.qualifiers:
                        gene.id = f.qualifiers['locus_tag'][0]+'__'+f.qualifiers['protein_id'][0]
                    elif 'pseudo' in f.qualifiers and 'locus_tag' in f.qualifiers:
                        gene.id = f.qualifiers['locus_tag'][0]+'__pseudo'
                    else:
                        gene.id = "randomID_"+str(random.randint(100000,999999))
                    gene.description = s.description
                    ffn.append(gene)
                if f.type in ['tRNA','rRNA','ncRNA']:
                    gene = f.location.extract(fnad[s.id])
                    if 'product' in f.qualifiers:
                        product = f.qualifiers['product'][0].replace(" ","_") if 'product' in f.qualifiers else "unknown"
                    else:
                        product = "unknown_product"
                    if 'locus_tag' in gene.id:
                        gene.id = f.qualifiers['locus_tag'][0]+'__'+product
                    else:
                        gene.id = "randomID_"+str(random.randint(100000,999999))+"__"+product
                    gene.description = s.description
                    frn.append(gene)

        out_ffn = info['ffn_lname']
        ffn_dir = os.path.dirname(out_ffn)
        if not os.path.exists(ffn_dir):
            os.mkdir( ffn_dir )
        if ffn:
            fo = bz2.BZ2File(out_ffn,"w")
            SeqIO.write( ffn, fo, "fasta")
            fo.close()
        out_frn = info['frn_lname']
        frn_dir = os.path.dirname(out_frn)
        if not os.path.exists(frn_dir):
            os.mkdir( frn_dir )
        if frn:
            fo = bz2.BZ2File(out_frn,"w")
            SeqIO.write( frn, fo, "fasta")
            fo.close()
        


    except Exception, e2: 
        logger.error("Something wrong in retrieving information for "+info['assembly_accession']+": "+str(e2) )



def get_assemblies( remote_file, outdir, sep = "\t" ):
    table = urllib2.urlopen( remote_file )
    table = [sline.strip().split(sep) for sline in table]
    table_header, table_content = [v.replace('#','').strip().lower() for v in table[0]],table[1:]
    assemblies = {}
    for line in table_content:
        line_d = dict(zip(table_header,line))
        if line_d['version_status'] != 'latest': continue
        ass_id = "G"+line_d['assembly_accession'].split("_")[1].split(".")[0]
        line_d['ass_id'] = ass_id
        line_d['fna_lname'] = outdir +"/fna/"+ass_id+".fna.bz2"
        line_d['ffn_lname'] = outdir +"/ffn/"+ass_id+".ffn.bz2"
        line_d['frn_lname'] = outdir +"/frn/"+ass_id+".frn.bz2"
        line_d['faa_lname'] = outdir +"/faa/"+ass_id+".faa.bz2"
        if line_d['taxid'] in taxids2taxonomy:
            line_d['taxonomy'] = taxids2taxonomy[line_d['taxid']]
        else:
            line_d['taxonomy'] = ""
            logger.warning(line_d['assembly_accession']+" ["+line_d['organism_name']+" "+line_d['infraspecific_name']+"] without known taxonomy!")
            continue
        if not line_d['taxonomy'] or line_d['taxonomy'].count("|") < 1 or line_d['taxonomy'].split('|')[1] in toescl:
            logger.warning(line_d['assembly_accession']+" ["+line_d['organism_name']+" "+line_d['infraspecific_name']+"] excluded from download because of uninteresting phyla!") 
            continue
        line_d['dwlf'] = add_protocol( NCBI_ftp + NCBI_assemblies_all + line_d['assembly_accession']  + "_" + line_d['asm_name'].replace(" ","_") ) 
        assemblies[ass_id] = line_d
    return assemblies


def get_table_by_assn( remote_file, key, sep = '\t' ):
    table = urllib2.urlopen( remote_file )
    table = [sline.strip().split(sep) for sline in table] 
    table_header, table_content = [v.replace('#','').strip().lower() for v in table[0]],table[1:] 
    ret_table = dict([(d[key],d) for d in [dict(zip(table_header,l)) for l in table_content]
                if key in d])
    return ret_table

if __name__ == '__main__':

    # ENSEMBLE ???? ftp://ftp.ensemblgenomes.org/pub/current/species.txt

    par = read_params(sys.argv)
    nproc = par['nproc'] 
 
    logger.info('Reading the taxonomy from '+par['taxonomy']+'... ')
    with open(par['taxonomy']) as inpf:
        taxids2taxonomy = dict([l.strip().split('\t')[1:] for l in inpf])
    logger.info('Done.')
    
    refseq_assemblies = get_assemblies( add_protocol(NCBI_ftp + NCBI_ASREFSEQ_file), par['out_dir'] )
    genbank_assemblies = get_assemblies( add_protocol(NCBI_ftp + NCBI_ASGENBANK_file), par['out_dir'] )
    
    if not os.path.exists( par['out_dir'] ):
        os.mkdir( par['out_dir'] )
        logger.warning(par['out_dir']+" does not exist. Creating it.")
    
    logger.info('Initializating the pool of '+str(nproc)+" downloaders for refseq genomes download")
    pool = mp.Pool( nproc )
    for ass, info in refseq_assemblies.items():
        info['outdir'] = par['out_dir']
        pool.map_async( refseq_dwl , [info] )
    pool.close()
    pool.join()
   
    logger.info('Producing the RefSeq output summary')
    dwl_summary = dwl_summarize(refseq_assemblies)

    pool = mp.Pool( nproc )
    for ass, info in genbank_assemblies.items():
        if info['ass_id'] in dwl_summary:
            logger.warning( info['ass_id'] + " already downloaded from RefSeq, skipping!")
            del genbank_assemblies[ass]
        else:
            info['outdir'] = par['out_dir']
            pool.map_async( refseq_dwl , [info] )
    pool.close()
    pool.join()
    
    logger.info('Producing the GenBank output summary')
    
    gb_dwl_summary = dwl_summarize(genbank_assemblies)
    
    logger.info('Writing the output summary')
    with open(par['out_summary'],"w") as outf:
        for i,(k,v) in enumerate(dwl_summary.items()+gb_dwl_summary.items()):
            if i == 0:
                outf.write( "\t".join(["#genome"] + list(sorted(v.keys()))) +"\n" )
            outf.write( "\t".join([k] + [v[kk] for kk in sorted(v.keys())]) +"\n" )

    
