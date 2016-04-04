# Pipelines

## Notes
The following applications are available in the iPlant Discovery Environment under the High Performance Computing section of the Apps catalog. They are accesible for scripting access via iPlant's Agave REST services platform, documented here http://preview.agaveapi.co/documentation/ .

Note also that the name iPlant is transitioning to be CyVerse.

## Single Animal Variant Pipeline (SAVP)

[Jump to SAVP Tutorial](/tutorial/tutorial_savp.md)

### VERSION HISTORY


#### Version 0.4.0: Coming soon
---
* Execution moved to Lonestar 5 supercomputer

| Agave Application ID | Label |
| -------------------- | ----- |
| savp-refprep-0.4.0u2 | Single Animal Variant Pipeline Reference Preparation Step |
| savp-bwa-aln-0.4.0u1 | Single Animal Variant Pipeline BWA ALN step |
| savp-bwa-mem-0.4.0u1 | Single Animal Variant Pipeline BWA MEM step |


#### Version 0.3.1: Released 01/07/2016
---
* Bug fixes in ianimal-savp
* ianimal-savp calls modules in preferred fashion

| Agave Application ID | Label |
| -------------------- | ----- |
| ianimal-savp-0.3.1u1 | Single Animal Variant Pipeline |

#### Version 0.3.0: Release 04/09/2015
---
* Adds clean-up options to prevent copying input and intermediate files to archive directory
* Adds BWA MEM app
* Adds BWA merge SAMs app
* Moves execution to Lonestar 4 supercomputer

| Agave Application ID | Label |
| -------------------- | ----- |
| ianimal-savp-refprep-0.3.0u1 | Single Animal Variant Pipeline Reference Preparation Step |
| ianimal-savp-merge-sams-1.129u1 | Single Animal Variant Pipeline merge SAM or BAM files |
| ianimal-savp-0.3.0u1 | Single Animal Variant Pipeline |
| ianimal-savp-bwa-aln-0.3.0u1 | Single Animal Variant Pipeline BWA ALN step |
| ianimal-savp-bwa-mem-0.3.0u1 | Single Animal Variant Pipeline BWA MEM step |

#### Version 0.2.0
---
* Private testing only
* Executes on Stampede supercomputer

### Source code

* savp-refprep-VERSION: Creates a tar file that includes the necessary index files for a reference genome
* savp-bwa-aln-VERSION: BWA ALN to align fastq files to reference genome
* savp-bwa-mem-VERSION: BWA MEM to align fastq files to reference genome
* savp-merge-sams-VERSION: Uses Picard mergeSamFiles to join sam or bam files
* savp-VERSION: Variant caller pipeline
* example-scripts

Each source application contains an example job file for submission via the jobs-submit command line tool
