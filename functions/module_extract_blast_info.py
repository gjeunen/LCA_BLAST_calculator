#! /usr/bin/env python3
import collections
from tqdm import tqdm
import os

def format_info(format):
    formatting = format.split('+')
    count = 0
    pident = 'NA'
    qcovs = 'NA'
    staxid = 'NA'
    qaccver = 'NA'
    for item in formatting:
        if item == 'pident':
            pident = count
        elif item == 'qcovs':
            qcovs = count
        elif item == 'staxid':
            staxid = count
        elif item == 'qaccver':
            qaccver = count
        count = count + 1

    return pident, qcovs, staxid, qaccver
        
def extract_blast(input, pident, qcovs, staxid, qaccver):
    LCA_dict = collections.defaultdict(list)
    taxid_list = []
    with tqdm(total = os.path.getsize(input)) as pbar:
        with open(input, 'r') as infile:
            for line in infile:
                pbar.update(len(line))
                pident_val = line.split('\t')[pident]
                qcovs_val = line.split('\t')[qcovs]
                staxid_val = line.split('\t')[staxid]
                qaccver_val = line.split('\t')[qaccver]
                LCA_dict[qaccver_val].append([staxid_val, pident_val, qcovs_val])
                taxid_list.append(staxid_val)
    taxid_set = set(taxid_list)

    return LCA_dict, taxid_set

def extract_nodes(nodes):
    taxids = {}
    with tqdm(total = os.path.getsize(nodes)) as pbar:
        with open(nodes, 'r') as infile:
            for line in infile:
                pbar.update(len(line))
                tax = line.split('\t|\t')[0]
                taxup = line.split('\t|\t')[1]
                rank = line.split('\t|\t')[2]
                taxids[tax] = [rank, taxup]

    return taxids


def extract_names(names):
    name_dict = {}
    with tqdm(total = os.path.getsize(names)) as pbar:
        with open(names) as infile:
            for line in infile:
                pbar.update(len(line))
                sc = line.split('\t')[6]
                if sc == 'scientific name':
                    taxid_name = line.split('\t|\t')[0]
                    name = line.split('\t|\t')[1].replace(' ', '_')
                    name_dict[taxid_name] = name

    return name_dict

def tax_lineage(taxid_set, ranks, taxid_info, names_info):
    rank_dict = {}
    rank_list = ranks.split('+')
    for item in rank_list:
        rank_dict[item] = 'yes'
    true_lineage = collections.defaultdict(list)
    true_lineage_test = {}
    for tax in tqdm(taxid_set):
        lineage = {}
        tax_line = {}
        correct_order_lineage = {}
        ktax = tax
        for i in range(10000000):
            if tax in taxid_info:
                lineage[taxid_info[tax][0]] = tax
                if taxid_info[tax][0] in ranks:
                    tax_line[taxid_info[tax][0]] = tax
                if tax == taxid_info[tax][1]:
                    break
                else:
                    tax = taxid_info[tax][1]
        for key in rank_dict:
            if key in tax_line:
                correct_order_lineage[key] = tax_line[key]
            else:
                correct_order_lineage[key] = 'nan'
        for k, v in correct_order_lineage.items():
            if v in names_info:
                true_lineage[ktax].append([k, v, names_info[v]])
            else:
                true_lineage[ktax].append([k, v, 'nan'])
        
    return true_lineage

def calculate_lca(LCA_dict, PIDENT, QCOV, RANKS, lineage):
    LCA_final = {}
    for item in tqdm(LCA_dict):
        temp = []
        count = 0
        for v in LCA_dict[item]:
            if float(v[1]) >= float(PIDENT) and float(v[2]) >= float(QCOV):
                spec_lin = ''
                ranks_test = RANKS.split('+')
                for index in range(len(ranks_test)):
                    spec_lin = spec_lin + lineage[v[0]][index][2] + ';'
                temp.append(spec_lin)
                count = count + 1
        blast_hit_number = count
        pract = os.path.commonprefix(temp)
        if pract.endswith(';'):
            result = pract.rstrip(';')
        else:
            result = pract.rsplit(';', 1)[0]
        if result == '':
            obtained_rank = 'NA'
        else:
            obtained_rank = ranks_test[len(result.split(';')) -1]
        LCA_final[item] = result + '___' + obtained_rank + '___' + str(blast_hit_number)
    
    return LCA_final
                
def read_otutable(FREQ):
    otutable = {}
    with tqdm(total = os.path.getsize(FREQ)) as pbar:
        with open(FREQ, 'r') as infile:
            for line in infile:
                pbar.update(len(line))
                line = line.rstrip('\n')
                otu_name = line.split('\t', 1)[0]
                samples = line.split('\t', 1)[1]
                otutable[otu_name] = samples
    
    return otutable

def combine_data(otutable, LCA_final, OUTPUT):
    LCA_total = {}
    for item in tqdm(otutable):
        if item.startswith('#'):
            LCA_total[item] = 'taxonomic_lineage___rank___blast_hit_number___' + otutable[item]
        elif item not in LCA_final:
            LCA_total[item] = 'na___na___na___' + otutable[item]            
        else:
            LCA_total[item] = LCA_final[item] + '___' + otutable[item]

    print(f'write data to {OUTPUT}')       
    with open(OUTPUT, 'w') as outfile:
        for item in tqdm(LCA_total):
            outfile.write(item + '\t' + LCA_total[item].split('___')[1] + '\t' + LCA_total[item].split('___')[2]+ '\t' + LCA_total[item].split('___')[0] + '\t' + LCA_total[item].split('___')[3] + '\n')
    