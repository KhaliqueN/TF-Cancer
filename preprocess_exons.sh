#!/usr/bin/sh

# select only exons
# awk -F '\t' '$3 == "exon"' ../../../public_data/Homo_sapiens.GRCh38.113.gtf > ../data/processed_exons.tmp
# awk -F '\t' '$3 == "CDS"' ../../../public_data/Homo_sapiens.GRCh38.113.gtf > ../data/processed_cds.tmp

awk -F '\t' '$3 == "exon"' ../../../public_data/Homo_sapiens.GRCh37.87.gtf > ../data/processed_exons.tmp
awk -F '\t' '$3 == "CDS"' ../../../public_data/Homo_sapiens.GRCh37.87.gtf > ../data/processed_cds.tmp

