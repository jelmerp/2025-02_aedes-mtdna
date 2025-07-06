#!/bin/bash

micromamba activate entrez-direct
efetch -db protein -format fasta_cds_na -id UYE93284.1