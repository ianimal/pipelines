# Pipelines



## Single Animal Variant Pipeline (SAVP)

* savp-refprep: Creates a tar file that includes the necessary index files
* savp-bwa-aln: BWA ALN to convert fastq files to bam
* savp-bwa-mem: BWA MEM to convert fastq files to bam
* savp: Variant caller pipeline

### VERSION HISTORY

#### Version 0.3.0
---
* Adds clean-up options to prevent copying input and intermediate files to output directory
* Adds BWA MEM app
* Moves execution to Lonestar 4 supercomputer



#### Version 0.2.0
---
* Private testing only
* Executes on Stampede supercomputer



