#!/bin/env python

import os
import re
import json
import uuid
import hashlib
from irods.session import iRODSSession

# This is a really annoying way to grab all the files from an iPlant directory
# Cleaner way is to use the Agave files endpoint. I was just trying out the iRODs bindings for Python
sess = iRODSSession(host='data.iplantcollaborative.org', port=1247, user='vaughn', password='CHANGEIT', zone='iplant')
coll = sess.collections.get("/iplant/home/schnabelr/fastq_cow")

formats={'P' : 'PE', 'S': 'SE', 'M': 'MP'}

# Find the *1_Unique.fastq.gz. We assume they are anchors.
manifest = {}
for obj in coll.data_objects:
    fname = obj.name
    match = re.search("1_Unique.fastq.gz$", fname)
    if match:
        fields = fname.split(".")
        animal = {  'breed' : fields[0].upper() ,
                    'lab_id' : fields[1].upper(),
                    'library' : fields[2][0:1].upper(),
                    'format' : formats [ fields[2][1:2].upper() ],
                    'replicate' : 'R' + fields[3].upper() }

        experiment_list = [ animal['lab_id'], animal['library'], animal['format'], animal['replicate'] ]
        experiment_id = "-".join(experiment_list)
        manifest[ experiment_id ] = animal
        manifest[ experiment_id ]['file1'] = fname

# Find the *2_Unique.fastq.gz. We assume they are optional
for obj in coll.data_objects:
    fname = obj.name
    match = re.search("2_Unique.fastq.gz$", fname)
    if match:
        fields = fname.split(".")
        animal = {  'breed' : fields[0].upper() ,
                    'lab_id' : fields[1].upper(),
                    'library' : fields[2][0:1].upper(),
                    'format' : formats [ fields[2][1:2].upper() ],
                    'replicate' : 'R' + fields[3].upper() }
        experiment_list = [ animal['lab_id'], animal['library'], animal['format'], animal['replicate'] ]
        experiment_id = "-".join(experiment_list)
        if experiment_id in manifest:
            manifest[ experiment_id ]['file2'] = fname

# Write out job scripts for savp-bwa-aln
appid1='savp-bwa-aln'
fpath='./schnabelr/cow/' + appid1
try:
    os.makedirs(fpath)
except OSError:
    if not os.path.isdir(fpath):
        raise

with open('lonestar-savp-bwa-aln-0.2.1.jsonx') as json_data:
    bwa_template = json.load(json_data)
    json_data.close()

for expt in manifest.keys():
    bwa_job = bwa_template
    # We will use this as a global identifier for downstream events
    id = hashlib.md5(expt).hexdigest()
    bwa_job['name'] = 'bwa.' + id
    bwa_job['archiveSystem'] = 'data.iplantcollaborative.org'
    bwa_job['archivePath'] = 'vaughn/lonestar-savp/schnabelr/cow/' + manifest[expt]['lab_id']
    bwa_job['inputs']['referenceBundle'] = "agave://lonestar-vaughn-work/iAnimal/umd_3_1_Y_Mito.tar"
    bwa_job['inputs']['inputSequence1'] = "agave://data.iplantcollaborative.org/schnabelr/fastq_cow/" + manifest[expt]['file1']
    if manifest[expt]['file2']:
        bwa_job['inputs']['inputSequence2'] = "agave://data.iplantcollaborative.org/schnabelr/fastq_cow/" + manifest[expt]['file2']
    bwa_job['parameters']['inputBarcode'] = expt

    with open(fpath + '/' + id + '.json', 'w') as outfile:
        json.dump(bwa_job, outfile, sort_keys = True, indent = 4)

