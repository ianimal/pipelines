The cow fasta files are there and should be ready for you to test.  Please note though that some files may be missing and that what is currently on iPlant is not “all” of the data because I’m still generating and processing data.  I also need to reconcile everything to make sure that I have copied everything and give you guys a detailed description of what everything represents and what has been done to it to get it to the point you currently see it in.

Having said that, you should be able to use it to do some testing but I wouldn’t try processing everything until I can get the above done.  Under my schnabelr home directory I have the following directories:

* fasta_cow
* fasta_dog
* fastq_cow
* fastq_dog
* results_fasta_cow
* results_fasta_dog
* results_fastq_cow
* results_fastq_dog

The first 4 contain the input data and I gave you read permission on the parent dirtectories.  I assume the permissions are inherited to the child objects.  You have write permissions on the results directories which is where you can put whatever output you have.  I assume that granting write permission also implies read permissions, is this correct?

File nomenclature
=================
All are in *FASTQ* format and named BBB.LLLLL.AF.RR.[1/2].fastq

BBB is the breed abbreviation
LLLLL is my unique lab_id
A represents the library the reads came from.  In most cases we generate two libraries per animal and they are designated A & B.  If additional libraries are used they are D, E etc.
F is the read format where P=paired, S=single, M=mate pair (large insert).

In most cases these two characters will be in the format AF but due to legacy reasons they may also be FA.  Additionally, as a general rule, any files for the same animal with the same library designation [A] came from the same library and any reads with a different library designation came from different libraries.  As usual, due to getting data from multiple sources who did not provide me with library info this may not always be true.  I have other metadata where I have a more accurate representation as to whether reads from different files came from the same source library.

RR represents the replicate run for a given library.  So if the “A” library was run three times there will be three files with “01”, “02”, “03”.
[1/2] represents the read direction where “1” are forward reads and “2” are reverse reads.

*_Unique means that the file has had exact duplicate reads removed.

*FASTA*
The format for the fasta files is very similar.
LLLLL.[pe.cor.fa].AA.[1/2].fasta

LLLLL is lab_id as before.
[pe.cor.fa] this part means that the file went through error correction using MaSuRCA.  For paired end data (small insert) this will be represented as pe.cor.fa.  For mate pair (large insert) this will be represented as sj.cor.clean.fa OR sj.cor.clean2.fa.  Long story as to why this is but it has to do with changes in the MaSuRCA program and when the sample was run.

AA is what I call the MaSuRCA abbreviation.  If you are familiar with the assembler, this is the two character abbreviation that must be supplied in the config file for each set of input files.  This is based on a combination of the AF fields from the fastq description above.  The first character represents the library the read originated from and the second character is the replicate number.  My MaSuRCA abbreviation is where I do a better job of identifying reads that originated from the same source library.  In general, any reads with the same first character came from the same library.  The problem is that data has been generated over the course of several years and I’ve received data from multiple collaborators so that’s why these are all generally applicable and not absolutes.

[1/2] read direction, same as fastq.

*_Unique same as above.  After the reads come out of error correction I again remove exact duplicates.  Keep in mind that this is quite different than the way other programs identify duplicates based on mapping info.  My definition of a duplicate read is one where the forward and reverse from pair1 are the *exact* same as the forward and reverse of read2.

*One last very important point concerning the FASTQ format.* Because our data has been generated literally over the course of the last 6-7 years, the format of the read name “@” line and the quality values are not uniform.  The read name portion runs the gamut from the different versions of the Illumina pipelines to having nothing to do with the standard Illumina format because someone renamed the read name portion before they sent me the data.  There are also some of the dog files that are actually SoLiD data that I ran through a conversion tool to convert from the SoLiD format to fastq.  The quality values (QV) are also an issue.  Again, since this data spans all the Illumina pipeline versions and I needed all the QV for an animal to be on the same scale for MaSuRCA I converted a lot of these to ASCII+33.  I did this by simply converting each QV to its ASCII integer and subtracting to make it base 33.  So if the Illumina version was ASCII+64 I would subtract 31 to make it ASCII+33.  If it was ASCI+66 I would subtract 33, if it was ASCII+35 I would subtract 2.  Therefore, if you are using the special designation of the ASCII+35 character # for anything you need to know that in general, all the QVs have been scaled to 33 ! as the min QV.  So if the original file was coded with a minimum QV of 35 # those got changed to 33 ! and therefore what is present in the file as 35 # actually was 37 in the original file.  For what you are doing it shouldn’t matter that much because in effect what I have done is make my QV more conservative.
