## Installation


Add directories to $PATH in ~/.profile or ~/.bash_profile

```
PATH = /path/to/transloc_pipeline/bin:/path/to/transloc_pipeline/R:$PATH
```

## Pre-processing libraries

```
$ cd data
$ TranslocPreprocess.pl metadata.txt preprocess/ --read1 pooled_R1.fq.gz --read2 pooled_R2.fq.gz
```

## Running the pipeline

```
$ TranslocWrapper.pl metadata.txt preprocess/ results/ --threads 2
```