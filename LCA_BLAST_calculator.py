#! /usr/bin/env python3

################################################
################ IMPORT MODULES ################
################################################
import argparse
import collections
import os
from tqdm import tqdm
import subprocess as sp
from functions.module_extract_blast_info import format_info, extract_blast, extract_nodes, extract_names, tax_lineage, calculate_lca, read_otutable, combine_data, export_lca

################################################
######### DOWNLOAD NCBI TAXONOMY INFO ##########
################################################
def ncbi_taxdump(args):
    SOURCE = args.source
    
    print('\ndownloading taxonomy information')
    url_taxdump = 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
    wget_help = sp.check_output('wget -h', shell=True)
    helstr=wget_help.decode('utf-8')
    if 'show-progress' in helstr:
        results = sp.run(['wget', url_taxdump, '-q', '--show-progress'])
    else:
        results = sp.run(['wget', url_taxdump, '-q'])
    results = sp.run(['tar', '-zxvf', 'taxdump.tar.gz'], stdout = sp.DEVNULL, stderr = sp.DEVNULL)
    print('removing intermediary files')
    if SOURCE == 'all':
        files_to_remove = ['citations.dmp', 'delnodes.dmp', 'division.dmp', 'gencode.dmp', 'merged.dmp', 'gc.prt', 'readme.txt', 'taxdump.tar.gz']
        for file in files_to_remove:
            os.remove(file)
            print()
    elif SOURCE == 'nodes':
        files_to_remove = ['citations.dmp', 'delnodes.dmp', 'division.dmp', 'gencode.dmp', 'merged.dmp', 'gc.prt', 'readme.txt', 'taxdump.tar.gz', 'names.dmp']
        for file in files_to_remove:
            os.remove(file)
            print()
    elif SOURCE == 'names':
        files_to_remove = ['citations.dmp', 'delnodes.dmp', 'division.dmp', 'gencode.dmp', 'merged.dmp', 'gc.prt', 'readme.txt', 'taxdump.tar.gz', 'nodes.dmp']
        for file in files_to_remove:
            os.remove(file)
            print()
    else:
        print('Please specify which taxonomy files need to be downloaded from NCBI by using the "SOURCE" parameter\n')


################################################
################# CREATE LCA ###################
################################################
def lca(args):
    INPUT = args.input
    OUTPUT = args.output
    BLAST = args.blast
    FREQ = args.freq
    NODES = args.nodes
    NAMES = args.names
    RANKS = args.ranks
    PIDENT = args.pident
    QCOV = args.qcov

    print(f'\nextracting information from {INPUT}')
    pident, qcovs, staxid, qaccver = format_info(BLAST)
    if pident == 'NA' or qcovs == 'NA' or staxid == 'NA' or qaccver == 'NA':
        print('BLAST output does not contain the necessary information. Fields needed to be included are: "qaccver", "pident", "qcovs", and "staxid"\n')
        exit()
    LCA_dict, taxid_set = extract_blast(INPUT, pident, qcovs, staxid, qaccver)
    print(f'extracting information from {NODES}')
    taxid_info = extract_nodes(NODES)
    print(f'extracting information from {NAMES}')
    names_info = extract_names(NAMES)
    print(f'generating taxonomic lineages for all tax IDs found in {INPUT}')
    lineage = tax_lineage(taxid_set, RANKS, taxid_info, names_info)
    print(f'calculating the LCA based on taxonomic lineages for tax IDs that pass the query cover percentage threshold of {QCOV} and percent identity threshold of {PIDENT}')
    LCA_final = calculate_lca(LCA_dict, PIDENT, QCOV, RANKS, lineage)
    if FREQ == None:
        print(f'exporting LCA calculations to {OUTPUT}')
        fin = export_lca(LCA_final, OUTPUT)
    else:
        print(f'reading in the frequency table {FREQ}')
        otutable = read_otutable(FREQ)
        print(f'combine LCA calculation with frequency table {FREQ}')
        final = combine_data(otutable, LCA_final, OUTPUT)


################################################
################### ARGPARSE ###################
################################################
def main():
    parser = argparse.ArgumentParser(description = 'BLAST LCA interpreter')
    subparser = parser.add_subparsers()

    ncbi_taxdump_parser = subparser.add_parser('ncbi_taxdump', description = 'download taxonomy information from the NCBI servers')
    ncbi_taxdump_parser.set_defaults(func = ncbi_taxdump)
    ncbi_taxdump_parser.add_argument('-s', '--source', help = 'specify which files need to be downloaded, including "all", "nodes", or "names"', dest = 'source', type = str, required = True)

    lca_parser = subparser.add_parser('lca', description = 'generate LCA for each OTU')
    lca_parser.set_defaults(func = lca)
    lca_parser.add_argument('-i', '--input', help = 'input file name', dest = 'input', type = str, required = True)
    lca_parser.add_argument('-o', '--output', help = 'output file name', dest = 'output', type = str, required = True)
    lca_parser.add_argument('-b', '--blast', help = 'format of BLAST results', dest = 'blast', type = str, required = True)
    lca_parser.add_argument('-f', '--freq', help = 'frequency table file name', dest = 'freq', type = str, required = False)
    lca_parser.add_argument('-t', '--nodes', help = 'filename of the "nodes.dmp" NCBI file', dest = 'nodes', type = str, default = 'nodes.dmp')
    lca_parser.add_argument('-n', '--names', help = 'filename of the "names.dmp" NCBI file', dest = 'names', type = str, default = 'names.dmp')
    lca_parser.add_argument('-r', '--ranks', help = 'taxonomic ranks to be included in the lineage', dest = 'ranks', type = str, default = 'superkingdom+phylum+class+order+family+genus+species')
    lca_parser.add_argument('-p', '--pident', help = 'percent identity threshold', dest = 'pident', type = str, required = True)
    lca_parser.add_argument('-q', '--qcov', help = 'query coverage threshold', dest = 'qcov', type = str, required = True)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()

