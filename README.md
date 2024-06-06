# scCounting
 A pipeline for read counting of multiple chemistries.

# Directory Structure

```
scripts/
├──count_reads.sh
├── utils/
│   ├──count_dropseq_read.sh
│   ├──count_smartseq_read.sh
│   ├──count_tenx_read.sh
temp/
results/
data/
logs/

```

# Input Data

Input reads must be the following format:
```
[Dataset ID]-[Sample_ID]_[OPTIONAL_IDs]_S[0-9]_L[0-9]*_R[1,2]_[0-9]*.fastq.gz
```

- must be .fastq.gz
- must be in bcl2 format
- Dataset ID must be of pattern: `v[0-9]-[0-9]*`
- Dataset ID must be found in manifest.json see `example_manifest.json`

e.g. `v1-01-SAMN33906275_SRR23958444_S1_L001_R1_001.fastq.gz`

# Output Data

Link to the raw and filetered .mtx directories generated in `outs/sample`