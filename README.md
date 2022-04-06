# LCA_BLAST_calculator

## Introduction

The LCA_BLAST_calculator program parses a BLAST output file and assigns a lowest common ancestor (LCA) to each submitted sequence using user-defined parameter thresholds for percent identity and query coverage. The LCA_BLAST_calculator program outputs the LCA taxonomic lineages in a tab-delimited document containing four columns, including (i) query sequence ID, (ii) taxonomic rank achieved, (iii) number of BLAST hits that passed the filter threshold and were used to calculate the LCA, and (iv) the taxonomic lineage of the LCA. If a frequency table (OTU table) is provided, the four columns can be inserted into the data table.

## Installing LCA_BLAST_calculator

LCA_BLAST_calculator is a command-line only toolkit running on typical Unix/Linux environments and is exclusively written in Python3. However, LCA_BLAST_calculator makes use of the subprocess module in python to run the `wget` command in bash syntax to circumvent python-specific idiosyncrasies and increase execution speed to download the NCBI taxonomy information files. Additionally, LCA_BLAST_calculator makes use of the folowing python modules, which might need to be installed separately:

1. argparse 1.1 or compatible
2. tqdm 4.59.0 or compatible

Thus far, LCA_BLAST_calculator is not incorporated in a pip or conda package. The easiest way to download LCA_BLAST_calculator would be to clone the GitHub repository using the code below. The LCA_BLAST_calculator code consist of a main python file and additional python functions that are stored in the `functions` subfolder. For proper execution of the code, it is important not to change the folder structure of LCA_BLAST_calculator.

```sh
git clone https://github.com/gjeunen/LCA_BLAST_calculator
cd LCA_BLAST_calculator
```

To check if the installation was successful, type in the following command to pull up the help information for LCA_BLAST_calculator.

```sh
./LCA_BLAST_calculator.py -h
```

Help information is also available for each of the modules included in LCA_BLAST_calculator and can be accessed by:

```sh
./LCA_BLAST_calculator.py MODULE -h
```

## Running LCA_BLAST_calculator

LCA_BLAST_calculator includes two modules:

1. `ncbi_taxdump`: download taxonomy information files from the NCBI website to be used for the LCA calculation
2. `lca`: calculates the LCA for each sequence in the BLAST output

### _1. ncbi_taxdump_
