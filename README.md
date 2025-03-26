---
title: "BLAST Search with MMseqs2"
author: "Your Name"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Introduction

This workshop will introduce **MMseqs2**, a fast alternative to BLAST for sequence similarity searches.

### Objectives:

- Run a sequence search
- Interpret the results

## Installation

MMseqs2 is already installed on Dardel. 

```
module load bioinfo-tools
module load MMseqs2/15-6f452
module load MMseqs2_data/latest
```

## Creating a Sequence Database

Before searching, either you can choose any of the pre-build databases from Dardel : 

```
module help MMseqs2_data/latest
```
The list of databases available with 'mmseqs databases' with their local UPPMAX locations, if available:

| Local UPPMAX location            | Name                  | Type       | Taxonomy | URL |
|----------------------------------|----------------------|-----------|---------|-------------------------------------------------------------|
| $MMSEQS2_DATA/UniRef100          | UniRef100            | Aminoacid  | true    | [UniRef100](https://www.uniprot.org/help/uniref)            |
| $MMSEQS2_DATA/UniRef90           | UniRef90             | Aminoacid  | true    | [UniRef90](https://www.uniprot.org/help/uniref)             |
| $MMSEQS2_DATA/UniRef50           | UniRef50             | Aminoacid  | true    | [UniRef50](https://www.uniprot.org/help/uniref)             |
| $MMSEQS2_DATA/UniProtKB          | UniProtKB            | Aminoacid  | true    | [UniProtKB](https://www.uniprot.org/help/uniprotkb)         |
| $MMSEQS2_DATA/UniProtKB_TrEMBL   | UniProtKB/TrEMBL     | Aminoacid  | true    | [UniProtKB TrEMBL](https://www.uniprot.org/help/uniprotkb)  |
| $MMSEQS2_DATA/UniProtKB_Swiss-Prot | UniProtKB/Swiss-Prot | Aminoacid  | true    | [Swiss-Prot](https://uniprot.org)                           |
| $MMSEQS2_DATA/NR                 | NR                   | Aminoacid  | true    | [NR](https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA)          |
| $MMSEQS2_DATA/NT                 | NT                   | Nucleotide | false   | [NT](https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA)          |
| $MMSEQS2_DATA/GTDB               | GTDB                 | Aminoacid  | true    | [GTDB](https://gtdb.ecogenomic.org)                         |
| $MMSEQS2_DATA/PDB                | PDB                  | Aminoacid  | false   | [PDB](https://www.rcsb.org)                                  |
| $MMSEQS2_DATA/PDB70              | PDB70                | Profile    | false   | [PDB70](https://github.com/soedinglab/hh-suite)            |
| $MMSEQS2_DATA/Pfam-A.full        | Pfam-A.full          | Profile    | false   | [Pfam-A.full](https://pfam.xfam.org)                        |
| $MMSEQS2_DATA/Pfam-A.seed        | Pfam-A.seed          | Profile    | false   | [Pfam-A.seed](https://pfam.xfam.org)                        |
| $MMSEQS2_DATA/Pfam-B             | Pfam-B               | Profile    | false   | [Pfam-B](https://xfam.wordpress.com/2020/06/30/a-new-pfam-b-is-released) |
| $MMSEQS2_DATA/CDD                | CDD                  | Profile    | false   | [CDD](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml) |
| $MMSEQS2_DATA/eggNOG             | eggNOG               | Profile    | false   | [eggNOG](http://eggnog5.embl.de)                            |
| $MMSEQS2_DATA/VOGDB              | VOGDB                | Profile    | false   | [VOGDB](https://vogdb.org)                                  |
| $MMSEQS2_DATA/dbCAN2             | dbCAN2               | Profile    | false   | [dbCAN2](http://bcb.unl.edu/dbCAN2)                         |
| $MMSEQS2_DATA/SILVA              | SILVA                | Nucleotide | true    | [SILVA](https://www.arb-silva.de)                           |
| $MMSEQS2_DATA/Resfinder          | Resfinder            | Nucleotide | false   | [Resfinder](https://cge.cbs.dtu.dk/services/ResFinder)      |
| $MMSEQS2_DATA/Kalamari           | Kalamari             | Nucleotide | true    | [Kalamari](https://github.com/lskatz/Kalamari)             |


or create you own database :

```
mmseqs createdb your_file.fasta your_file_db
```

## Searching

# Objectif : find out what is the closest origin of this nucleotide sequence: 

```
>LC371381.1
ATGTTGGGTTCAGACAAATTCTCATGTTTCTCTGATCAGCACAGAGCTAGATCTCCAAGTCCTACTGACA
GGAAGGATAAGAAAAACCACACGAACAAGCTTCGAGAGCTGGCTTTGCTGATCCCTGTGACCATGAAGAC
CAGAGACAAGAAGTACACCAAGAAGGAGATCCTGTTGCGTGTCCTGCACTACATCCAGTACCTCCAGAGA
AACATTGACATGACCAAGGCCTTGCTCAAGCTCCACAGCAGCAATGGCAAAGGTAGATTTGTGGGGCCAG
GTTTGAACCCATCTGCTGGCCAGACACAGCAGCAACACTCCACTCCCTCCAGCTCTCAGAAGCCAAGCCT
TTGGAGCACCTCTTCAAAACCTCGAAAGAAGAAGTTTACCCGAGTGTCAGAGCATCCATCCTGGCCCTAT
AATCCTCGACGCTCTCTAGCTCTGGACCAGGCTGAAAATCCTAACACCATACATCCAGGCCTAAAGGAAG
AAAATGAGGAATGTGCCACCTATCCAGGAGTCCTCAGCCCTAGCACCTACCCCACAACTGAACCATCTGT
GTCTGAAGGCGATGGACAAGGGGCCCAGTTGGTGTTCTTGGACATGGCTCAGAACATCTTTGCCTATGAC
ATCCTAAGTGATCATGCTGTAGAAGTCCAGGGTGGAGAGCCCAATGCTGACATCAAGGTTCAGAGGTCCT
TTTTCCTCACCAGGGCACAGCCCTGTGTCAGTTCTTGCAGGCAGAAGCTATTCTTGTGTACTTCCAGTGA
GGCAGACAAAGAAGCCCCAGACTCTGACCCCTGGCTTCCTGTTTGGACCTCCGAGGACAGCCCCAATGGG
AGCCCGCTGGCTTTGGGGTCTTCCCAGATCAATACTTGGCATGTGGCAGACTACCTGAATGAGATCTTAG
GAGTCAGCTCTTCCCTCTTCAGCTCCCCAAGCAAAATCCTGCCGGATCATGTCCTAGAAGATGGCACCTA
CTTTCTGACTGAAGGTCTCTTGGAGTCTTCACCTGCTACCTGTGAAGTGGAGAGCCCACAAGAGAAGGAA
GTATCCTCTGAAGGCCCCACAGGCCCACCTAACTTCCAGTCCTCTGTCTCACTGGACCACTGCTACCTGT
CGCTGAGTGAGAATGTCAAGGTGCTGTCTAACTGTGGCTCCAGCTCAGAGTCCACAGACACAGAGTCTCT
GTGGGGACAGGAGGAGGTGAGAGGGGTGGCCACCTACCCAAGGAGGACAGGGGTCTTCTCCCTTGCTCAG
CCTTCTGTCCTGCAGCAGGCCAACCCTGAGGGATTACAGACCTCAAGTGATGAGGACAGAGACTACACAT
GGACACCTACTGGCCAGTCTTCTGGCCTGCCAGTAGCCAGCAAGAAGATCAAGAAGGTCCAGGCAAGCCA
GGGCCCCGTGAAGCCCAAAGACAGCAGAAAAGCCTGCCCTGGCCAGGTGAAGAAAAAGTGTGTTAATGGC
TTCATCATGTTCTGCAGGATGAACCGGAAGCAGTACATCCGAGCCTGTCCTGGAACTGCATCCACAGCTG
CCACCAAGGATCTGGCTCAACTGTGGCGAGGGATGACTCTGGAGGAAAAGAAACCATACTGCACTAAGGC
ACGCAGGTTCAGCCGCCAGAACAACCGCATTGTGAAGCAGGAGAACTCCAGCAGTGAGGACGACGATGGG
GAGACCCCCAAGCCCTTCTACCAGCTATTAGCTGAGAAGGCCCAAGTGTCTTTGGGCCTCACCTCACTGC
CCACCCCTAACTGCCAGTGA

```

Before searching, we need to convert the FASTA file containing query sequences and target sequences into a sequence DB. 
We can use the query database `examples/QUERY.fasta` and target database `UniProtKB_Swiss-Prot` to test the search workflow:

```
mmseqs createdb examples/QUERY.fasta queryDB
```

These calls should generate five files each, e.g. queryDB, queryDB_h and its corresponding index file queryDB.index, queryDB_h.index and queryDB.lookup from the FASTA QUERY.fasta input sequences.

The alignment consists of two steps the prefilter and alignment. To run the search, type:

```
mmseqs search queryDB $MMSEQS2_DATA/UniProtKB_Swiss-Prot resultDB tmp --threads 7
```

Search as standard does compute the score only. If you need the alignment information add the option “-a”. The speed and sensitivity of the search can be adjusted with -s parameter and should be adapted based on your use case (see setting sensitivity -s parameter). A very fast search would use a sensitivity of -s 1.0, while a very sensitive search would use a sensitivity of up to -s 7.0.

The output can be customized with the --format-output option

Then, convert the result database into a BLAST tab formatted file (option -m 8 in legacy blast, -outfmt 6 in blast+):

```
mmseqs convertalis --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tlen,qlen,tcov,qcov,taxid,taxname,taxlineage" queryDB $MMSEQS2_DATA/UniProtKB_Swiss-Prot resultDB resultDB.m8 --format-mode 0
```


## Results 

| Accession | Protein ID | Identity | Alignment Length | Mismatches | Gap Openings | Query Start | Query End | Subject Start | Subject End | E-value | Bit Score | Subject Length | Query Coverage | Subject Coverage | Taxonomy ID | Organism | Lineage |
|-----------|------------|----------|------------------|------------|--------------|-------------|-----------|---------------|-------------|---------|-----------|----------------|----------------|------------------|------------|----------|--------|
| LC371381.1 | A0A5K7RLP0 | 96.000   | 1767             | 24         | 0            | 1           | 1767      | 1             | 589        | 0.000E+00 | 1190     | 589            | 1.000          | 0.998          | 10090   | Mus musculus | -_cellular organisms;d_Eukaryota;-_Opisthokonta;k_Metazoa;-_Eumetazoa;-_Bilateria;-_Deuterostomia;p_Chordata;-_Craniata;-_Vertebrata;-_Gnathostomata;-_Teleostomi;-_Euteleostomi;-_Sarcopterygii;-_Dipnotetrapodomorpha;-_Tetrapoda;-_Amniota;c_Mammalia;-_Theria;-_Eutheria;-_Boreoeutheria;-_Euarchontoglires;-_Glires;o_Rodentia;-_Myomorpha;-_Muroidea;f_Muridae;-_Murinae;g_Mus;-_Mus;s_Mus musculus |
| LC371381.1 | C9JSJ3     | 46.300   | 1686             | 304        | 0            | 37          | 1722      | 64            | 630        | 9.136E-146 | 476      | 638            | 0.889          | 0.953          | 9606    | Homo sapiens | -_cellular organisms;d_Eukaryota;-_Opisthokonta;k_Metazoa;-_Eumetazoa;-_Bilateria;-_Deuterostomia;p_Chordata;-_Craniata;-_Vertebrata;-_Gnathostomata;-_Teleostomi;-_Euteleostomi;-_Sarcopterygii;-_Dipnotetrapodomorpha;-_Tetrapoda;-_Amniota;c_Mammalia;-_Theria;-_Eutheria;-_Boreoeutheria;-_Euarchontoglires;o_Primates;-_Haplorrhini;-_Simiiformes;-_Catarrhini;-_Hominoidea;f_Hominidae;-_Homininae;g_Homo;s_Homo sapiens |



we need to format our sequence file into an MMseqs2 database:

```bash
mmseqs createdb example.fasta exampleDB
```

## Running a Sequence Search

To search a query sequence against the database:

```bash
mmseqs search query.fasta exampleDB resultDB tmp
mmseqs convertalis query.fasta exampleDB resultDB result.tsv
```

This will generate a results file (`result.tsv`) containing hits.

## Interpreting Results in R

We can analyze the results in R:

```{r}
library(tidyverse)

# Load results
df <- read_tsv("result.tsv", col_names = c("query", "target", "evalue", "bit_score"))

# Inspect the top hits
head(df)

# Plot bit score distribution
ggplot(df, aes(x = bit_score)) +
  geom_histogram(binwidth = 5, fill = "blue", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of Bit Scores", x = "Bit Score", y = "Count")
```

## Conclusion

MMseqs2 is a powerful alternative to BLAST that allows rapid sequence similarity searches. This workshop demonstrated how to set up MMseqs2, create a database, perform searches, and analyze results using R.
