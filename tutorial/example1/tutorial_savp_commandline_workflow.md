
Running from the command line requires an account and likely a compute allocation on the system on which you are running.  A key advantage of running from the command line is that outputs between steps can also be stored locally on the execution system.  This eliminates the need to copy the files back to the execution system for subsequent steps, saving significant time and bandwidth.

First, however, you will need to set up your local environment.  

#### Preparing your local environment

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

#### Command line job submission

You will want to familiarize yourself with the general job submission process

* [Running a job with Agave](https://github.com/iPlantCollaborativeOpenSource/iplant-agave-sdk/blob/master/docs/iplant-first-app-job.md)

#### Workflow example 1 part 1 - preparing the reference bundle 

Use the following json file for your job submission.  

```sh
{
    "jobName": "savp-refprep-030-local",
    "softwareName": "ianimal-savp-refprep-0.3.0u1",
    "processorsPerNode": 12,
    "requestedTime": "08:00:00",
    "memoryPerNode": 24,
    "nodeCount": 1,
    "batchQueue": "normal",
    "archive": false,
    "archivePath": "YOUR-ARCHIVE-PATH",
    "inputs": {
		"referenceFasta": "agave://data.iplantcollaborative.org/shared/iplant_training/savp_example_data/e-coli-K-12.fa"
    },
    "parameters": {
    	"cleanupParameter": true
    },
    "outputs": {
    }
}
```

Fill in *YOUR-ARCHIVE-PATH* or leave the string empty.

Launch the job will the following command:

```sh
jobs-submit -F savp-refprep-030-local-job.json
```

Copy the job ID reported back to you.  Subsequent steps refer to it as *JOBID-REFPREP* 


#### Workflow example 1 part 2 - creating the bam file

Use the following json file for your job submission.  

```sh
{
    "jobName": "savp-bwa-mem-030-local",
    "softwareName": "ianimal-savp-bwa-mem-0.3.0u1",
    "processorsPerNode": 12,
    "requestedTime": "08:00:00",
    "memoryPerNode": 24,
    "nodeCount": 1,
    "batchQueue": "normal",
    "archive": false,
    "archivePath": "YOUR-ARCHIVE-PATH",
    "inputs": {
        "inputSequence1": "agave://data.iplantcollaborative.org/shared/iplant_training/savp_example_data/SRR2601691_1.fastq",
        "inputSequence2": "agave://data.iplantcollaborative.org/shared/iplant_training/savp_example_data/SRR2601691_2.fastq",
        "referenceBundle": "https://agave.iplantc.org/jobs/v2/JOBID-REFPREP/outputs/media/e-coli-K-12.tar"
    },
    "parameters": {
        "inputBarcode": "SRR2601691",
        "cleanupParameter": true
    },
    "outputs": {
    }
}
```

Fill in *YOUR-ARCHIVE-PATH* or leave the string empty.

Replace *JOBID-REFPREP* with the job ID you copied in the previous step.

Launch the job will the following command:

```sh
jobs-submit -F savp-bwa-mem-030-local-job.json
```

Copy the job ID reported back to you.  Subsequent steps refer to it as *JOBID-BWA* 


#### Workflow example 1 part 3 - variant calling

Use the following json file for your job submission.  

```sh
{
    "jobName": "savp-031-local",
    "softwareName": "ianimal-savp-0.3.1u1",
    "processorsPerNode": 12,
    "requestedTime": "08:00:00",
    "memoryPerNode": 24,
    "nodeCount": 1,
    "batchQueue": "normal",
    "archive": false,
    "archivePath": "YOUR-ARCHIVE-PATH",
    "inputs": {
        "inputBam": "https://agave.iplantc.org/jobs/v2/JOBID-BWA/outputs/media/SRR2601691.mem.bam",
        "inputBamBai": "https://agave.iplantc.org/jobs/v2/JOBID-BWA/outputs/media/SRR2601691.mem.bam.bai",
        "referenceBundle": "https://agave.iplantc.org/jobs/v2/JOBID-REFPREP/outputs/media/e-coli-K-12.tar",
	    "knownVariants": "agave://data.iplantcollaborative.org/shared/iplant_training/savp_example_data/known.variants.vcf"
    },
    "parameters": {
    	"regionName": "",
    	"runUnifiedGenotyper": true,
    	"runPlatypus": true,
    	"runMpileup": false,
    	"cleanupInputs": true,
    	"cleanupIntermediates": true
    },
    "outputs": {
    }
}
```

Fill in *YOUR-ARCHIVE-PATH* or leave the string empty.

Replace *JOBID-REFPREP* and *JOBID-BWA* with the appropriate job IDs you previously copied.

Launch the job will the following command:

```sh
jobs-submit -F savp-031-local-job.json
```

* [Go back to SAVP Discovery Environment tutorial](../tutorial_savp.md)

* [Go back to SAVP README](../../README.md)
