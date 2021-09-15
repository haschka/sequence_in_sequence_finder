

Please cite: [![DOI](https://zenodo.org/badge/406458347.svg)](https://zenodo.org/badge/latestdoi/406458347)

https://doi.org/10.5281/zenodo.5510252

# A FASTA - Fourier Transform based Sequence in Sequence Finder

(c) Thomas Haschka 2021 - for a license see the provided license file.

This tool finds a fasta sequence in a fasta sequence using
the fourier transform.

# Installation

This tool requires the FFTW3 library with double precision support. In
order to build this tool also the -dev package (i.e. header files) are
required. On debian/ubuntu based distributions you may install this package
using:

```sudo apt-get install build-essential fftw3 fftw3-dev```

Once all dependencies are installed you may compile this tool on your system
like this:

```gcc -O3 -march=native binary_array.c dataset.c seq_in_seq.c -lfftw3_threads -lfftw3 -lm -o seq_in_seq```

Once you have compiled this tool you can call it with the following arguments:
```
  [fasta] sequence to search sequence in 
  [fasta] sequence to be searched for 
  [int] number of sequence to search in 
        in the first fasta file. Beginning with 1
  [int] 1 created a wisdom file during this run (if you do not have one) 
        0 use a wisdom file already available 
  [wisdom] path of a wisdom file, always to be specified 
  [int] cutoff - number of sequence changes to accept 
  [string] chromosome specifier 
  [int] number of threads fftw can use 
```

The first run of the program might be long as you are supposed to generate a
*wisdom file* for the sequence to be searched in. This allows you to find the
optimal fourier transformation for your system as well as the size of your
problem. Once this is done a sequence in a sequence can be found very quickly
using such a *wisdom file*. I.e. if you are searching frequently in hg38
I suggest you create wisdom files for each chromosome.
The *chromosome specifier* is the first string that will be written in your
bed coordinates output. I.e. if you search on hg38-chr1 and type chr1 you can
upload your output to the UCSC Genome Browser and visulize the locations of
the sequence that you are searching for.

# Memory usage

The current version requires a lot of memory. A machine with 32 GB of memory
should nevertheless be sufficant to search all of Hg38-chr1. On a 64 GB
searching Hg38-chr1 poses no problem. On machines with less ram you might cut
your fasta files into pieces and perform the search piecwise. A fix for this
might come up in a future version that performs this cut automatically on low
memory machines.
