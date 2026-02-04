This was Levitan et al 2020 


cerevisiae:

cerevisiae_viable_annotations.txt, cerevisiae_inviable_annotations.txt - Sc viable/inviable data was taken from yeastgenome.org, date of retrieval inside the files - 26/09/2016.

pombe:

Pombe data was taken from the pombase.org website.
FYPOviability.tsv - retrieved from ftp://ftp.ebi.ac.uk/pub/databases/pombase/pombe/Phenotype_annotations/ on 5/10/16.
schizosaccharomyces_pombe.chr.gff3 - retrieved from ftp://ftp.ebi.ac.uk/pub/databases/pombase/pombe/Chromosome_Dumps/gff3/. Dated 31/05/2016.

albicans:

Calb orthologs retrieved from http://www.candidagenome.org/download/homology/ on 20/03/2017.

Albicans' gene2go was derived from http://www.candidagenome.org/download/go/ (retrieved from http://www.candidagenome.org/download/go/ on 21/10/2016). All non-albicans annotations were removed, goatools.associations.read_gaf was used on the remaining dataset, and a reverse dictionary was generated and stored as the above file.

go-basic.obo - retrieved from http://geneontology.org/ontology/go-basic.obo on 21/10/2016.


Updated by Joshua Ayelazuno 25/01/2026

InterProScan6 domain annotations for Candida glabrata BG2 (NCBI assembly)

======================================================================

WHY THIS FILE EXISTS
--------------------

The original Levitan et al. (2020) transposon-pipeline relies on several precomputed
annotation dependencies (e.g., feature tables and InterProScan outputs), but the
methods used to generate those dependencies were not fully documented. In
particular, the pipeline references an `iprscan.out` file for Candida albicans via
hard-coded paths, without describing how the InterProScan analysis was run or how
those files can be regenerated.

Because CGD-style resources are legacy and no longer actively maintained, we
standardize all reference dependencies on NCBI assemblies and annotations. This
ensures reproducibility, portability across species/strains, and long-term
maintainability of the pipeline.

This README documents how InterProScan domain annotations for Candida glabrata BG2
are generated using InterProScan 6 via Nextflow and Apptainer. These annotations
serve as a reproducible replacement for the legacy CGD-derived iprscan dependencies
used in the original pipeline.


WHAT IS GENERATED
-----------------

Input: all reference genome materials were downloaded from NCBI (25/01/2026)
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/217/725/GCA_014217725.1_ASM1421772v1/

- NCBI protein FASTA:
  GCA_014217725.1_ASM1421772v1_protein.faa

Outputs (written to this directory):
- InterProScan6 outputs (TSV and/or GFF, as produced by the workflow)
- GO term annotations (--goterms)
- Pathway annotations (--pathways)

All outputs are generated directly into:

~/workspace/HermesTnseq/transposon-pipeline/dependencies/glabrata/reference_genome/


COMPUTE ENVIRONMENT REQUIREMENTS
--------------------------------

InterProScan6 is run using a containerized workflow:

- Nextflow >= 25.04.6
- Apptainer (Singularity-compatible container runtime)

No manual installation of InterProScan binaries or databases is required. All
software and reference databases are downloaded automatically by the workflow.

Verified on Argon:
- apptainer version 1.3.2-1.el7
- nextflow >= 25.04.6


DIRECTORY LAYOUT AND CACHES
---------------------------

Large reference data and caches are stored outside of git-tracked code, under
~/workspace/, to avoid accidental deletion and unnecessary duplication.

Directories used:

- InterProScan databases:
  ~/workspace/interproscan6_data        (~15–20 GB, reused across runs)

- Apptainer container cache:
  ~/workspace/apptainer_cache

- Nextflow work directory (scratch space):
  ~/workspace/nextflow_work

- Nextflow internal cache/history:
  ~/workspace/.nextflow


ONE-TIME SETUP (PER SESSION)
----------------------------

Create required directories:

mkdir -p ~/workspace/interproscan6_data
mkdir -p ~/workspace/apptainer_cache
mkdir -p ~/workspace/nextflow_work
mkdir -p ~/workspace/.nextflow

Export environment variables so Nextflow and Apptainer use these locations:

export NXF_APPTAINER_CACHEDIR=~/workspace/apptainer_cache
export NXF_WORK=~/workspace/nextflow_work
export NXF_HOME=~/workspace/.nextflow

(Optional) These variables can be made persistent using:

conda activate tnseq
conda env config vars set NXF_APPTAINER_CACHEDIR=$HOME/workspace/apptainer_cache
conda env config vars set NXF_WORK=$HOME/workspace/nextflow_work
conda env config vars set NXF_HOME=$HOME/workspace/.nextflow
conda deactivate
conda activate tnseq


REPRODUCIBILITY TEST RUN (RECOMMENDED)
-------------------------------------

This test confirms that Nextflow, Apptainer, and InterProScan6 are correctly
configured and that databases and containers can be downloaded:

cd ~/workspace

nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0 \
  -profile apptainer,test \
  --datadir ~/workspace/interproscan6_data


FULL BG2 INTERPROSCAN RUN
------------------------

Run from ~/workspace (not from inside the datadir). On HPC systems, submit this
command to a compute node rather than running on a login node.

cd ~/workspace

nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0 \
  -profile apptainer \
  --input ~/workspace/HermesTnseq/transposon-pipeline/dependencies/glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_protein.faa \
  --datadir ~/workspace/interproscan6_data \
  --outdir ~/workspace/HermesTnseq/transposon-pipeline/dependencies/glabrata/reference_genome \
  --goterms \
  --pathways


NOTES ON IDENTIFIERS AND DOWNSTREAM USE
--------------------------------------

- InterProScan annotations are produced at the protein level; the primary key is
  the protein FASTA header / accession.
- For downstream Tn-seq analyses, proteins will be mapped to genes using NCBI
  annotations (GFF/feature table), typically via:
  - product_accession or non-redundant_refseq
  - locus_tag (used as the stable gene_id)
- CGD-specific identifiers are intentionally avoided.


PROVENANCE AND VERSIONING
-------------------------

- Workflow: ebi-pf-team/interproscan6
- Workflow release pinned to: 6.0.0
- Containers and databases are pulled automatically by Nextflow
- Reproducibility is ensured by pinning the workflow version and preserving the
  InterProScan data directory

To reproduce these results elsewhere, use the same workflow release (-r 6.0.0),
the same input protein FASTA, and either reuse or regenerate the InterProScan data
directory.
