# Tutorials

## Single Animal Variant Pipeline (SAVP)

The SAVP can be run via the web-based Discovery Environment or from the command line. In either case, you will need to first [register a CyVerse/iPlant] (https://user.iplantcollaborative.org/register/) account. Once your account is active, you can launch the [Discovery Environment] (https://de.iplantcollaborative.org/de) from your web browser.  

If you have files to upload, select the **Data** icon to open a file browser window.  then under *Navigation* select your username and then select *New Folder* from the *File* menu option.   Name your folder, browse into the folder, and then you can use *Simple Upload from Desktop* in the *Upload* menu.  Selecting *Bulk Upload* instead will send you to information on how to upload a large number of files and/or large files.

You do not need to upload any files to work through this tutorial. This tutorial demonstrates its steps on *E. coli* datasets, all of which are readily accessible in the Discovery Environment.

*NB: SAVP currently is running on TACC Lonestar4 and will be transitioning to Lonestar5 in early 2016.  Downtime during transition is expected to be minimal.* 

### Access and execution through the CyVerse (formerly iPlant) Discovery Environment

The different executable steps of SAVP are accessible as Apps in the DE. Select the **Apps** icon to open an analysis window. Under *Category*, scroll down and select *High-Performance Computing*.  

#### Step 1 - Preparing the reference bundle

Different parts of the SAVP require a reference fasta and associated files. As a first step, these associated files are created from the reference fasta and bundle into a tar file. This tar file only needs to be created once for a specific reference fasta.

In the **Apps** window under *Workspace* select *Single Animal Variant Pipeline Reference Preparation Step 0.3.0*.  This will open a window for the app. 

Select **Inputs** and then select the **Browse** button on the right to open a navigation window.  From this window you can select a reference fasta you previously uploaded, or the sample file provided for this tutorial which is found by navigating into *Community Data* then *iplant_training* then *savp_example_data* and selecting *e-coli-K-12.fa.gz*.  Select the *OK* button at the bottom of the window.  

Now select the *Launch Analysis* button.  A new folder for the output of this analysis will be created after launching the analysis.  Unless you specified a different location, it will be found in */iplant/home/yourusername/analyses/* which you can browse to using your **Data** window.

Upon successful completion, the new reference tar bundle file will appear in the output folder.  Job-related output and error text files will also appear.

Summary of Step 1:

* App name: *Single Animal Variant Pipeline Reference Preparation Step 0.3.0*
* Input: reference fasta
* Output: tar bundle containing reference fasta and associated files
* Parameters: clean-up


#### Step 2 - Creating bam files using bwa mem or bwa aln

Both bwa mem and bwa aln apps are currently supported.  As in Step 1, use the **Apps** window.  Choose either *Single Animal Variant Pipeline BWA ALN step 0.3.0* or *Single Animal Variant Pipeline BWA MEM step 0.3.0* and select it so that its app window pops up.

Select **Inputs** and then select the **Browse** button on the right to open a navigation window.  From this window you can select a single fasta/fastq or paired set of fasta/fastq files you previously uploaded.  The sample files provided for this tutorial which are found by navigating into *Community Data* then *iplant_training* then *savp_example_data* and selecting first *SRR2601715_1.fastq* and then *SRR2601715_2.fastq* as its pair.  Also select *e-coli-K-12.tar* which will be found in the output folder created by Step 1.

Select **Parameters** and specify a *Barcode*.  The shared base name of the inputs is usually a good choice.

Now select the *Launch Analysis* button.  As with Step 1, a new folder for the output of this analysis will be created after launching the analysis.  Unless you specified a different location, it will be found in */iplant/home/yourusername/analyses/* which you can browse to using your **Data** window.

Upon successful completion, *barcode*.bam and *barcode*.bam.bai files will appear in the output folder.  Job-related output and error text files will also appear.

Summary of Step 2:

* App name: *Single Animal Variant Pipeline BWA ALN step 0.3.0* or *Single Animal Variant Pipeline BWA MEM step 0.3.0*
* Inputs: single fasta files or paired fasta files; reference tar bundle from Step 1
* Outputs: bam and bam.bai files
* Parameters: barcode; clean-up


#### Step 3 - Variant callers

Three different variant callers are currently supported through the DE: GATK UnifiedGenotyper; Platypus CallVariants; and Samtools mpileup/vcfutils.

Select **Inputs** and then select the **Browse** button on the right to open a navigation window.  From this window you will need to select in turn the outputs of previous steps: the bam file, the bam.bai file, and the reference tar. You can also use a known variant file, of which a sample for this tutorial which are found by navigating into *Community Data* then *iplant_training* then *savp_example_data* and selecting *known.variants.vcf*.

Select **Parameters** and, if desired, specify a specific region (e.g. Chromosome) to variant call. specify a *Barcode*.  You can also modify the default selection of variant callers.  However, this is not recommended.

Now select the *Launch Analysis* button.  As with Steps 1 and 2, a new folder for the output of this analysis will be created after launching the analysis.  Unless you specified a different location, it will be found in */iplant/home/yourusername/analyses/* which you can browse to using your **Data** window.

Upon successful completion, depending on which variant callers you selected to run, appropriate *.vcf* files will appear in the output folder.  Job-related output and error text files will also appear.

Summary of Step 3:

* App name: *Single Animal Variant Pipeline 0.3.1* 
* Inputs: bam and bam.bai files from Step 2; reference tar bundle from Step 1; known variants vcf (optional)
* Outputs: vcf files depending on variant caller selection
* Parameters: region; variant callers selected; clean-up


### Running SAVP workflows from the command line

* [Go to command line workflow instructions example #1](example1/tutorial_savp_commandline_workflow.md)

* [Go back to README](../README.md)
