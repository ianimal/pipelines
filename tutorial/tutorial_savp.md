# Tutorials  (coming soon)

[Jump back to README](README.md)

## Single Animal Variant Pipeline (SAVP)

The SAVP can be run via the web-based Discovery Environment or from the command line.  
In either case, you will need to first [register a CyVerse/iPlant] (https://user.iplantcollaborative.org/register/) account.  
Once your account is active, you can launch the [Discovery Environment] (https://de.iplantcollaborative.org/de) from your web browser.  
Select the **Data** icon to open a file browser window.  
If you have files to upload, then under *Navigation* select your username and then select *New Folder* from the *File* menu option.   
Name your folder, browse into the folder, and then you can use *Simple Upload from Desktop* in the *Upload* menu.  
Selecting *Bulk Upload* instead will send you to information on how to upload a large number of files and/or large files.
You do not need to upload any files to work through this tutorial.
This tutorial demonstrates its steps on e coli datasets, all of which are readily accessible in the Discovery Environment.

** NB: SAVP currently is running on TACC Lonestar4 and will be transitioning to Lonestar5 in early 2016.  Downtime during transition is expected to be minimal. **


### Access and execution through the CyVerse (formerly iPlant) Discovery Environment

Select the **Analysis** icon to open an analysis window.
Under *Navigation, scroll down and select High Performance Computing apps.

#### Step 1 - Preparing the reference bundle

Different parts of the SAVP require a reference fasta and associated files.  
As a first step, these associated files are created from the reference fasta and bundle into a tar file.
This tar file only needs to be created once for a specific reference fasta.
In the analysis window, look for ______ and select it (may require a double click) so that its app window pops up.
To specify the input file,  _______
For the purposes of this tutorial, you can select the e_coli.fa file.  
Select ______
Select ______ to submit the job.   You will receive a notification when it has finished.

App name: 
Input: reference fasta
Output: tar bundle containing reference fasta and associated files
Parameters: clean up



#### Step 2 - Creating bam files using bwa mem or bwa aln

Both bwa mem and bwa aln apps are currently supported.
As in the previous step, in the analysis window look for ________ and select it so that its app window pops up.

App name: 
Input: single fasta files or paired fasta files
Input: reference tar bundle from Step 1
Input: known variants vcf file (optional)
Outputs: bam and bam.bai files
Parameters:



#### Step 3 - Running the pipeline









### Running SAVP workflows from the command line

Running from the command line requires an account and likely a compute allocation on the system on which you are running.  
A key advantage of running from the command line is that outputs between steps can also be stored locally on the execution system.  
This eliminates the need to copy the files back to the execution system for subsequent steps, saving significant time and bandwidth.

First, however, you will need to set up your local environment.  

#### Step 0 - Prepare local environment

Many of the steps to do this overlap with the steps required to develop your own apps for the DE.

You will need to work through the following once:
* [Initial Assumptions](https://github.com/iPlantCollaborativeOpenSource/iplant-agave-sdk/blob/master/docs/iplant-assumptions.md)
* [Installing the SDK](https://github.com/iPlantCollaborativeOpenSource/iplant-agave-sdk/blob/master/docs/install-sdk.md)
* [Client and key creation](https://github.com/iPlantCollaborativeOpenSource/iplant-agave-sdk/blob/master/docs/client-create.md)
* [Obtain a token](https://github.com/iPlantCollaborativeOpenSource/iplant-agave-sdk/blob/master/docs/set-token.md)

If done correctly, then for subsequent actively, you will only need to refresh your tokens using the following command:

```sh
auth-tokens-refresh -S
```


#### Step 1 - Preparing the reference bundle (coming soon)


#### Step 2 - Creating bam files (coming soon)


#### Step 3 - Running the pipeline (coming soon)


