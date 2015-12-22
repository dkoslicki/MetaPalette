#!/local/cluster/bin/python2.7
# !/usr/bin/env python

# ==============================================================================
# generate_taxonomy.py
#
# Authors: Nicola Segata (nicola.segata@unitn.it), Roman Stolyarov (r.m.stolyarov@gmail.com)
#
# Creates taxonomy file and/or reduced taxonomy file for all prokaryotes, eukaryotes, and viruses
# for which genome sequence data is present in Refseq database.
# ==============================================================================

__author__ = 'Nicola Segata (nicola.segata@unitn.it), Roman Stolyarov (r.m.stolyarov@gmail.com), additions by David Koslicki (david.koslicki@math.oregonstate.edu)'
__version__ = '1.1.2'
__date__ = '6 Aug 2013'

import sys
import collections
import re
import copy 
import argparse
import pickle
import tarfile
import urllib2
import StringIO
import logging

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree as BTree
from Bio.Phylo.BaseTree import Clade as BClade


def read_params(args):
    parser = argparse.ArgumentParser(description='')
    arg = parser.add_argument
    """
    arg( '--nodes', metavar='nodes', required = True, type=str,
         help="The nodes.dmp from the NCBI taxonomy")
    arg( '--names', metavar='names', required = True, type=str,
         help="The names.dmp from the NCBI taxonomy")
    """
    arg( '--ncbi_taxdump', type=str, required = False,
         default = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
         help="The remote NCBI taxdumpfile (default is ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)")
    arg( '--output', metavar='out', required = True, default = None, type = str, 
         help="The output taxonomy")
    arg( '--output_red', metavar='out_red', required = True, default = None, type = str, 
         help="The output taxonomy with reduced and fixed number of taxonomic levels")
    arg( '--pickle', metavar='pickle', required = True, default = None, type = str, 
         help="The pickle output file")
    return vars(parser.parse_args())




class Nodes:
    #
    # Format of nodes.dmp from RefSeq documentation
    #
    # ---------
    # 
    # This file represents taxonomy nodes. The description for each node includes 
    # the following fields:
    # 
    #   tax_id                  -- node id in GenBank taxonomy database
    #   parent tax_id               -- parent node id in GenBank taxonomy database
    #   rank                    -- rank of this node (superkingdom, kingdom, ...) 
    #   embl code               -- locus-name prefix; not unique
    #   division id             -- see division.dmp file
    #   inherited div flag  (1 or 0)        -- 1 if node inherits division from parent
    #   genetic code id             -- see gencode.dmp file
    #   inherited GC  flag  (1 or 0)        -- 1 if node inherits genetic code from parent
    #   mitochondrial genetic code id       -- see gencode.dmp file
    #   inherited MGC flag  (1 or 0)        -- 1 if node inherits mitochondrial gencode from parent
    #   GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
    #   hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
    #   comments                -- free-text comments and citations
    #

    def __init__( self ):
        pass

    def __init__( self, nodes_dmp, tax_ids_to_names = None ):
        #Go through every line of Nodes file to construct tree. tmp_nodes will be a dictionary pointing from the taxid to its clade
        tmp_nodes = {}
        #with open( nodes_dmp_file ) as inpf:
        for line in nodes_dmp:
            ( tax_id, parent_tax_id, rank, embl_code, division_id, inherited_div_flag,
            genetic_code_id, inherited_GC_flag, mitochondrial_genetic_code, inherited_MGC_flag,
            GenBank_hidden_flag, hidden_subtree_root_flag, comments ) = line[::2]
   
    #For every entry in Nodes (every location in the tree) create clade containing the scientific name and pointer to the parent node.
    #Specify the rank of the clade and the taxonomic ID of the root.
            name = (tax_ids_to_names[int(tax_id)] if tax_ids_to_names else None)

            clade = BClade( clades = [], name = name )
            clade.parent_tax_id = int(parent_tax_id)
            clade.rank = re.sub(r'\W+', '', rank).strip("_")
            clade.tax_id = int(tax_id)         
            #clade.accession = accessions[clade.tax_id] if clade.tax_id in accessions else []
            
    #Set clade status values to "True" for sequence data and "final" or "draft" if it appears in accessions (taxid -> name, status, accessions)
            #if clade.tax_id in accessions:
            #    clade.sequence_data = True
            #    clade.status = clade.accession['status']
        
            tmp_nodes[clade.tax_id] = clade 
                
                # can add any other info in node.dmp
    #Build the tree using all the clades (iterate through clades using tmp_nodes)
        self.tree = BTree()
        for node in tmp_nodes.values():
            # node = parent is the trick from NCBI to identify the root
            if node.tax_id == node.parent_tax_id:
                self.tree.root = node
                continue
            parent = tmp_nodes[node.parent_tax_id]
            parent.clades.append( node )

    #Recursively goes through all clades in the tree. Each clade root gets list of all accessions in the clade.
    def add_internal_accessions( self, clade = None ):
        if not clade:
            clade = self.tree.root

        clade.all_accessions = [] + ([clade.accession] if clade.accession else [])

        for child in clade.clades:
            clade.all_accessions += self.add_internal_accessions( child )
        return clade.all_accessions
       
    #Recursively go through tree, remove references to clades that have no accession information in any of their nodes.
    def remove_subtrees_without_accessions( self, clade = None ):
        if not clade:
            clade = self.tree.root
        clade.clades = [c for c in clade.clades if len(c.all_accessions)]
        for c in clade.clades:
            self.remove_subtrees_without_accessions( c )
   
    #Recursively go through the tree, and at each node remove references to child clades pertaining to plasmid DNA.
    def remove_plasmids( self, clade = None ):
        if not clade:
            clade = self.tree.root
        clade.clades = [c for c in clade.clades if 'plasmid' not in c.name]
        for c in clade.clades:
            self.remove_plasmids( c )
   
    def print_tree( self, out_file_name, reduced = False ):

        tree = self.reduced_tree if reduced else self.tree 

        #to_print = tree.find_clades({"sequence_data": True})
      
        ranks2code = {'superkingdom':'k','phylum':'p','class':'c','order':'o','family':'f','genus':'g','species':'s','taxon':'t'}

        def red_rank( rank ):
            if reduced and rank in ranks2code:
                return ranks2code[rank]
            return rank

        

        with open(out_file_name,"w") as outf:
            
            def trac_print_t( clade, names = None ):
                if names is None:
                    if clade.name == 'root':
                        names = ""
                    else:
                        names = red_rank(clade.rank)+'__'+clade.name 
                else:
                    names += ("|" if names else "")+red_rank(clade.rank)+'__'+clade.name 

                #if clade.is_terminal():
                if clade.tax_id is not None and clade.name != 'root':
                    outf.write("\t".join( [ clade.name,
                                            #t.accession['status'],
                                            #",".join(t.accession['gen_seqs']),
                                            str(clade.tax_id),
                                            #t.accession['code'],
                                            #",".join(t.accession['accession']),
                                            #str(t.accession['len']),
                                            names
                                            ]    )+"\n")

                if not clade.is_terminal():
                    for c in clade.clades:
                        trac_print_t( c, names )
            trac_print_t( tree.root )


            """
            for t in tree.get_terminals():
                tax = "|".join([red_rank(p.rank)+'__'+p.name for p in tree.get_path( t )])
                outf.write("\t".join( [ t.name,
                                        #t.accession['status'],
                                        #",".join(t.accession['gen_seqs']),
                                        str(t.tax_id),
                                        #t.accession['code'],
                                        #",".join(t.accession['accession']),
                                        #str(t.accession['len']),
                                        tax
                                        ]    )+"\n")
            """
    
    def get_tree_with_reduced_taxonomy( self, superkingdom = "Bacteria", logger = None ):
        reduced_tax_levels = ['superkingdom','phylum','class','order','family','genus','species']
        self.reduced_tree = copy.deepcopy(self.tree)
      

        def add_noranks( clade ):
            if not clade.rank in reduced_tax_levels:
                clade.rank = 'norank'
            for c in clade.clades:
                add_noranks( c )

        def remove_noranks( clade ):

            run = True

            while run:
                run = False
                new_clades = []
                for c in clade.clades:
                    #if len(c.clades) and c.rank not in reduced_tax_levels:
                    #if not hasattr(c,"sequence_data") and c.rank not in reduced_tax_levels:
                    if len(c.clades) and c.rank not in reduced_tax_levels:
                        run = True
                        c.rank = "norank"
                        new_clades += c.clades
                    else:
                        new_clades.append(c)
                    #if hasattr(c,"sequence_data") and c.rank not in reduced_tax_levels:
                    if c.rank not in reduced_tax_levels:
                        c.rank = "norank"
                clade.clades = new_clades
            for c in clade.clades:
                if len(c.clades):
                    remove_noranks( c )

        def add_taxa( clade ):
            #if clade.rank == "norank" and hasattr(clade,"sequence_data"):
            if clade.rank == "norank" and clade.is_terminal():
                clade.rank = "taxon"

            """
            if not len(clade.clades) and clade.accession:
                if clade.rank == 'species':
                    newclade = copy.deepcopy( clade )
                    clade.accession = [] 
                    clade.sequence_data = False
                    newclade.rank = "taxon"
                    clade.clades = [newclade]
            """

            for c in clade.clades:
                add_taxa( c )


        def add_internal_missing_levels( clade, lev = 1 ):
            if clade.rank == "taxon":
                return
            cur_lev = reduced_tax_levels.index(clade.rank) if lev > 0  else 0

            jumps, to_add = [], ""
            for i,c in enumerate(clade.clades):
                if c.rank == 'taxon':
                    continue
                next_lev = reduced_tax_levels.index(c.rank)
                if next_lev == cur_lev + 1 or c.rank == 'superkingdom':
                    add_internal_missing_levels( c )
                    continue

                for i,l in enumerate(reduced_tax_levels[:-1]):
                    if clade.rank == l and c.rank != reduced_tax_levels[i+1]:
                        jumps.append( c )
                        to_add =  reduced_tax_levels[i+1]
            if jumps:
                children_ok = [c for c in clade.clades if c not in jumps]
                newclade = copy.deepcopy( clade )
                newclade.clades = jumps
                clade.clades = [newclade]+children_ok
                newclade.rank = to_add 
                newclade.name = clade.name if "_noname" in clade.name else clade.name+"_noname"
                add_internal_missing_levels( newclade )

        def reduce_double_taxa( clade ):
            if clade.rank == 'species' and len(clade.clades):
                torem = []
                for c in clade.clades:
                    if c.rank == 'taxon':
                        if len(c.clades):
                            clade.clades += c.clades
                            torem.append( c )
                clade.clades = [c for c in clade.clades if c not in torem]
                return
            for c in clade.clades:
                reduce_double_taxa( c )

        #add_noranks( self.reduced_tree.root )
        
        logger.info("Removing noranks from the taxonomy")
        remove_noranks( self.reduced_tree.root )
        #add_noranks( self.reduced_tree.root )
        logger.info("Adding taxon names to the taxonomy")
        add_taxa( self.reduced_tree.root )
        logger.info("Adding nternal missing taxonomic levels")
        add_internal_missing_levels( self.reduced_tree.root, lev = -1 )
        logger.info("Removing duplicated taxa")
        reduce_double_taxa( self.reduced_tree.root )

    #def save( self, out_file_name ):
    #    self.tree = self.tree.as_phyloxml()
    #    Phylo.write( self.tree, out_file_name, "phyloxml")

class Names:
    #
    # Format of names.dmp from RefSeq documentation
    #
    # ---------
    #
    # Taxonomy names file has these fields:
    #
    #   tax_id                  -- the id of node associated with this name
    #   name_txt                -- name itself
    #   unique name             -- the unique variant of this name if name not unique
    #   name class              -- (synonym, common name, ...)
    #

    def __init__( self, names_dmp ):
        #Read from file names.dmp, get information in every field
        self.tax_ids_to_names = {}
        #with open( names_dmp_file ) as inpf:
        for line in names_dmp:
            tax_id, name_txt, unique, name_class = line[::2]
                
                # extracting scientific names only (at least for now) which are unique!
        #tax_ids_to_names relates taxid to the sole scientific name of the organism
            if name_class == "scientific name":
                name = re.sub(r'\W+', '_', name_txt).strip("_")
                #self.tax_ids_to_names[ int(tax_id) ] = name
                self.tax_ids_to_names[ int(tax_id) ] = tax_id + "_" + name #Let's hope this will print out the tax ID's instead of the names, it did!!

            
    def get_tax_ids_to_names( self ):
        return self.tax_ids_to_names


if __name__ == '__main__':
    par = read_params(sys.argv)

    logging.basicConfig(level=logging.INFO, stream=sys.stdout, 
                        format = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logger = logging.getLogger(sys.argv[0])

    logger.info('Downloading and reading the NCBI taxdump file from '+par['ncbi_taxdump'])
    try:
        loc = urllib2.urlopen( par['ncbi_taxdump'] )
        compressedFile = StringIO.StringIO(loc.read())
        tarf = tarfile.open(fileobj=compressedFile)
        for m in tarf.getmembers():
            if m.name == "names.dmp":
                names_buf = (l.strip().split('\t') for l in tarf.extractfile( m ))
            if m.name == "nodes.dmp":
                nodes_buf = (l.strip().split('\t') for l in tarf.extractfile( m ))
    except Exception, e:
        logger.error("Error in downloading, extracting, or reading "+par['ncbi_taxdump']+": "+str(e))
        sys.exit()

    logger.info('names.dmp and nodes.dmp succeffully downloaded, extracted, and read')

     
    logger.info('Processing the names.dmp file to extract the taxonomic names')
    names = Names( names_buf )
    logger.info('Finished reading names')
    logger.info('Loading the taxonomic tree from nodes.dmp')
    tax_tree = Nodes( nodes_buf, names.get_tax_ids_to_names() ) #, accessions.get_accessions() ) 
    logger.info('Finished reading the initial taxonomic tree')
    #tax_tree.print_tree( "aaa.txt"  )

    #tax_tree = Nodes( par['nodes'], names.get_tax_ids_to_names() ) #, accessions.get_accessions() ) 
    #tax_tree.add_internal_accessions()
    #tax_tree.remove_plasmids()
    #tax_tree.remove_subtrees_without_accessions()

    #tax_tree.tree.root = tax_tree.tree.root.clades[0]
    
    if par['output']:
        logger.info('Exporting the original unedited NCBI taxonomy to: '+par['output'])
        tax_tree.print_tree( par['output'] )
        logger.info(par['output'] + ' saved.')

    logger.info('Processing the taxonomy for consistency and a fixed number of taxonomic levels')
    tax_tree.get_tree_with_reduced_taxonomy(logger=logger)
    logger.info('Finished postprocessing the taxonomy')

    if par['output_red']:
        logger.info('Exporting the edited NCBI taxonomy to: '+par['output_red'])
        tax_tree.print_tree( par['output_red'], reduced = True)
        logger.info(par['output_red'] + ' saved.')
    if par['pickle']: 
        with open(par['pickle'], "w" ) as out:
            logger.info('Pickled taxonomy saved to: '+par['pickle'])
            pickle.dump(tax_tree, out, -1)
            logger.info(par['pickle'] + ' saved.')
        


