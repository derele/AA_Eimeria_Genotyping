# AA_Eimeria_Genotyping

Parasite data for this has been published in "Generalist Eimeria
species in rodents: Multilocus analyses indicate inadequate resolution
of established markers" (2020) Víctor Hugo Jarquín-Díaz, Alice Balard,
Anna Mácová, Jenny Jost, Tabea Roth von Szepesbéla, Karin Berktold,
Steffen Tank, Jana Kvičerová, Emanuel Heitlinger
https://doi.org/10.1002/ece3.5992

The repository for this manuscript is found at
https://github.com/VictorHJD/AA_Eimeria_Genotyping (the parent repo of
this fork)

# Host genotyping - the pipeline for this repository

One benefit of the multi-amplicon system is that we can simultaneously
target host and parasite genetic markers. Here we detail the host
genotyping performed in parallel to this parasite genotyping.

In essence we find no hybrid distinction in the ASVs. No hybrid
singature... This can be deduced from the (slightly messy) clustering
of ASVs in the single R script. The project is discontinued at this
point.

I integrated reading of the sequence files from SRA here. This is
making even the bioinformatics part of the analysis reproducible on
any computer. This approach might be worthwile to ingegrate into other
projects. 
