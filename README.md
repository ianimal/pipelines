# Pipelines

## Single Animal Variant Pipeline (SAVP)

* savp-refprep: Creates a tar file that includes the necessary index files
* savp-bwa-aln: BWA ALN to convert fastq files to bam
* savp-bwa-mem: BWA MEM to convert fastq files to bam
* savp-merge-sams: Uses Picard mergeSamFiles to join sam or bam files
* savp: Variant caller pipeline

### VERSION HISTORY

#### Version 0.3.0: Release 04/09/2015
---
* Adds clean-up options to prevent copying input and intermediate files to archive directory
* Adds BWA MEM app
* Adds BWA merge SAMs app
* Moves execution to Lonestar 4 supercomputer

#### Version 0.2.0
---
* Private testing only
* Executes on Stampede supercomputer



