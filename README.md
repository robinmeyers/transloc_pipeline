# Installation

## Download code

```
cd ~
git clone https://github.com/robinmeyers/transloc_pipeline
```

Add directories to $PATH in ~/.profile or ~/.bash_profile

```
echo 'export PATH = ~/transloc_pipeline/bin:~/transloc_pipeline/R:$PATH' >> ~/.bash_profile
```

## Reference Genomes

```
mkdir -p ~/genomes/bowtie2_indexes
cd ~/genomes/bowtie2_indexes
wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip
unzip hg19.zip
mkdir ~/genomes/hg19
bowtie2-inspect hg19 > ~/genomes/hg19/hg19.fa 
```

Add environment variables to in ~/.profile or ~/.bash_profile

```
echo 'export BOWTIE2_INDEXES = ~/genomes/bowtie2_indexes' >> ~/.bash_profile
echo 'export GENOME_DB = ~/genomes' >> ~/.bash_profile
```



## Dependencies

### Bowtie2

Install [Bowtie2](http://bowtie-bio.sourceforge.net//bowtie2/index.shtml)

### Perl >= 5.16

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

Other modules:

- 


# Running the Pipeline

## Pre-processing libraries

```
$ cd data
$ TranslocPreprocess.pl metadata.txt preprocess/ --read1 pooled_R1.fq.gz --read2 pooled_R2.fq.gz
```

## Running the pipeline

```
$ TranslocWrapper.pl metadata.txt preprocess/ results/ --threads 2
```


