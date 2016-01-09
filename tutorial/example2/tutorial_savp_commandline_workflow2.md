# Data Locality with Agave

See json job files in this directory
[AN.156.AP.01.mem.job.json](AN.156.AP.01.mem.job.json)
[AN.156.BP.01.mem.job.json](AN.156.BP.01.mem.job.json)
[AN.156.BP.02.mem.job.json](AN.156.BP.02.mem.job.json)
[AN.156.BP.03.mem.job.json](AN.156.BP.03.mem.job.json)
[AN.156.merge.job.json](AN.156.merge.job.json)
[AN.156.savp.job.json](AN.156.savp.job.json)

## Copy the reference data to tacc-global.iplantcollaborative.org

## Submit stage 1

Each of 4 bwa-mem jobs, with reference data coming from tacc-global.iplantcollaborative.org. Archive is set to false.

Capture job IDs on submit...

* AN.156.AP.01.mem.job.json   7041323658071961115-e0bd34dffff8de6-0001-007
* AN.156.BP.01.mem.job.json   6888036275277721115-e0bd34dffff8de6-0001-007
* AN.156.BP.02.mem.job.json   6531553989709721115-e0bd34dffff8de6-0001-007
* AN.156.BP.03.mem.job.json   6395661224464281115-e0bd34dffff8de6-0001-007

## Create stage 2: Merge

Use outputs as input, rather than archived BAMs at iDS. Archive is set to false.

* "https://agave.iplantc.org/jobs/v2/7041323658071961115-e0bd34dffff8de6-0001-007/outputs/media/AN.156.AP.01.mem.bam",
* "https://agave.iplantc.org/jobs/v2/6888036275277721115-e0bd34dffff8de6-0001-007/outputs/media/AN.156.BP.01.mem.bam",
* "https://agave.iplantc.org/jobs/v2/6531553989709721115-e0bd34dffff8de6-0001-007/outputs/media/AN.156.BP.02.mem.bam",
* "https://agave.iplantc.org/jobs/v2/6395661224464281115-e0bd34dffff8de6-0001-007/outputs/media/AN.156.BP.03.mem.bam"

Submit, get job ID

## Same for stage 3

Use outputs as input, rather than archived merged BAM. Archive set to true if you wish!

## Now, manually move the appropriate files back to iDS via a set of files-import commands

* [Go to command line workflow instructions example #1](../example1/tutorial_savp_commandline_workflow.md)

* [Go back to SAVP Discovery Environment tutorial](../tutorial_savp.md)

* [Go back to SAVP README](../../README.md)
