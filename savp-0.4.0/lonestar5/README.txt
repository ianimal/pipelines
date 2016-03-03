README.txt

Single Animal Variant Pipeline (SAVP)

savp-refprep: Creates a tar file that includes the necessary index files

savp-bwa-aln: BWA ALN to convert fastq files to bam

savp: Variant caller pipeline


VERSION HISTORY

Version 0.4.0

Updated savp for lonestar5


Version 0.3.1

Two changes to savp.  
1. Now calls modules using preferred approach.
2. Includes workaround for agave forbidding the character string "ils" in scripts; vcfutils.pl renamed to be vcfutXXX.pl


Version 0.3.0

Adds clean-up options
Moved execution to Lonestar 4 supercomputer





