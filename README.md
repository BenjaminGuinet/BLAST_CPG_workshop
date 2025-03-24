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
- Install MMseqs2
- Create a sequence database
- Run a sequence search
- Interpret the results

## Installation

First, we need to install MMseqs2. You can do this by downloading the binary:

```bash
wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
mkdir -p ~/bin && tar xvf mmseqs-linux-avx2.tar.gz -C ~/bin
export PATH=~/bin/mmseqs/bin:$PATH
```

To make this change permanent, add it to your `~/.bashrc`:

```bash
echo 'export PATH=~/bin/mmseqs/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

## Creating a Sequence Database

Before searching, we need to format our sequence file into an MMseqs2 database:

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
