# LCA_BLAST_calculator

The LCA_BLAST_calculator program parses a BLAST output file and assigns a lowest common ancestor (LCA) to each submitted sequence using user-defined parameter thresholds for percent identity and query coverage. The LCA_BLAST_calculator program outputs the LCA taxonomic lineages in a tab-delimited document containing four columns, including (i) query sequence ID, (ii) taxonomic rank achieved, (iii) number of BLAST hits that passed the filter threshold and were used to calculate the LCA, and (iv) the taxonomic lineage of the LCA. If a frequency table (OTU table) is provided, the four columns can be inserted into the data table.

## Download LCA_BLAST_calculator

The easiest way to download LCA_BLAST_calculator would be to clone the GitHub repository using the code below.

```sh
git clone https://github.com/gjeunen/LCA_BLAST_calculator
cd LCA_BLAST_calculator
```
