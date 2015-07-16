## Installation

### Dependencies

Install [Bioperl](http://www.bioperl.org/wiki/Installing_BioPerl_on_Unix)

```
$ cpan
```
Find the most recent version

```
cpan> d /bioperl/
Distribution    CJFIELDS/BioPerl-1.6.901.tar.gz
Distribution    CJFIELDS/BioPerl-1.6.923.tar.gz
Distribution    CJFIELDS/BioPerl-1.6.924.tar.gz
```

Now install:

```
cpan> install CJFIELDS/BioPerl-1.6.924.tar.gz
```
or
```
cpan> force install CJFIELDS/BioPerl-1.6.924.tar.gz
```



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


