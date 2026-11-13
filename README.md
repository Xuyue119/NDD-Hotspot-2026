# NDD Missense Hotspot Analysis

This repository provides the R scripts used to generate the main and supplementary figures for the NDD missense-hotspot study, together with the derived plotting datasets required by those scripts. It also documents the command-line workflows used for sequence-based and structure-based clustering with Denovonear, mutation-rate estimation, and 3D-CLUMP.

The following figures were generated using the provided R scripts: Figure 2, Figure S2, Figure 3F, Figure 4, Figure S4, Figure 5, Figure S5, Figure 6, and Figure S6. The data/ directory contains the datasets required for figure generation.

The following figures were prepared with no code involved: Figure 1, Figure 3A-E, Figure 3G-H, Figure S1, Figure S3.


The following commands document the tool invocations used in this study.

Denovonear: 1D sequence-based clustering, 3D structure-based clustering, and mutation-rate estimation

GTF="${REF_DIR}/gencode.v38.annotation.gtf"
FASTA="${REF_DIR}/GRCh38.p13.genome.fa"
STRUCTURES="/path/to/UP000005640_9606_HUMAN_v6.tar"

1D: sequence-based clustering
Denovonear cluster \
  --in grch38.dnms.txt \
  --gencode gencode.v38.annotation.gtf \
  --fasta GRCh38.p13.genome.fa \
  --out output.txt \
  --genome-build grch38

3D: structure-based clustering
Denovonear cluster-structure \
  --in grch38.dnms.txt \
  --structures UP000005640_9606_HUMAN_v6.tar \
  --gencode gencode.v38.annotation.gtf \
  --fasta GRCh38.p13.genome.fa \
  --out output.txt \
  --genome-build grch38

Gene-level mutation-rate estimation
python gene_rates.py \
  --genes gene_ids.txt \
  --gencode gencode.v38.annotation.gtf \
  --fasta GRCh38.p13.genome.fa \
  --out gene_rates.tsv

Site-level missense mutation-rate estimation
python site_rates.py \
  --genes gene_ids.txt \
  --gencode gencode.v38.annotation.gtf \
  --fasta GRCh38.p13.genome.fa \
  --gene-rates gene_rates.tsv \
  --genome-build GRCh38 \
  --out missense_sites.tsv.gz


3D-CLUMP (3D-CLUstering by Mutation Position): 3D structure-based clustering

python 3D-CLUMP.modified.py \
  -f inputfile \
  -p protein_lengths \
  -s structurepath \
  -m 5 \
  -t \
  -z 100000