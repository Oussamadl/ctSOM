# ctSOM
analysis pipeline for somatic diagnosis on ctDNA
# Introduction
ctSOM is a variant analysis pipeline designed for circulating tumor DNA (ctDNA) obtained through non-invasive liquid biopsy, serving diagnostic, prognostic, and therapeutic purposes.

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible.

The strength of this pipeline is attributed to the strategic integration of Unique Molecular Identifiers (UMIs). These molecular barcodes effectively mitigate the challenges posed by PCR duplicates and sequencing errors, contributing to the pipeline's enhanced accuracy and sensitivity in identifying genetic variants

# Pipeline summary
the pipeline can currently perform the following
. Extraction of the UMI from insert reads and generate consensus reads from them (fgbio )
. Adapter and quality trimming (fastp).
. Map Reads to Reference (bwa mem).
. Process BAM file (gatk SortSam, gatk MergeBamAlignment)
. convertion of file format ( gatk SamToFastq, gatk FastqToSam, vcf2tsv)
. quality control (samtools depth, picard collectHsMetrcis, multiQc)
. variant calling (vardict)
. Annotation (VEP)
