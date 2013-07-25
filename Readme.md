Burrows-Wheeler Aligner (BWA) is an efficient program that aligns
relatively short nucleotide sequences against a long reference sequence
such as the human genome. This app runs BWA to map **letterspace** reads to
a reference genome and produce mappings.

This app accepts one or more gtable object of type Reads. The reads in all
of these tables will be mapped to the genome and output in a single gtable
of type Mappings. (See http://wiki.dnanexus.com/Types/Mappings for more
information on the Mappings type.) Mappings objects can be then used as
inputs to certain variation calling apps, mappings QC apps, etc.; they can
also be visualized in the genome browser.

## Genome indexing

BWA requires a special processing on the genome, called indexing; this
processing can take several hours for long genomes. If the genome given in
the input is not indexed for BWA, the app will automatically index it and
include an indexed version in the output, for future use. When possible,
please run this app with an indexed genome to avoid re-indexing. DNAnexus
provides several pre-indexed genomes in the 'Reference Genomes' public
project.

## Mapping algorithm

BWA implements two different algorithms, both based on the Burrows-Wheeler
Transform (BWT). The first algorithm, called 'aln', is designed for short
query sequences up to ~200bp with low error rate (<3%). It does gapped
global alignment with respect to the queries, supports paired-end reads,
and is one of the fastest short read alignment algorithms to date while
also visiting suboptimal hits. The second algorithm, called 'bwasw', is
designed for long reads with more errors. It performs heuristic
Smith-Waterman-like alignment to find high-scoring local hits (and thus
chimera). On low-error short queries, 'bwasw' is slower and less accurate
than the first algorithm, but on long queries, it is better. Using a value
of 'auto' will automatically choose a suitable algorithm based on the
length of the reads in the input (however, this is currently not
implemented, and using 'auto' will lead to the 'aln' algorithm always being
chosen).

## Parallelization

This app parallelizes itself by dividing the input into chunks of a certain
size (a certain number of reads), and mapping each chunk individually.
Lower chunk sizes lead to higher levels of parallelization, reducing the
wall-clock time that one has to wait for the app to finish. However, lower
chunks sizes may also increase the cost of running the app in the cloud, as
they lead to a higher number of chunks, each of which adds a constant
processing overhead. The default value is 25 million reads per chunk;
DNAnexus suggests caution when experimenting with this parameter.

## Low-level parameters

Users familiar with the BWA executable can directly manipulate the
parameters that are used for the `bwa aln`, `bwa samse`, `bwa sampe` and
`bwa bwasw` calls. These parameters are: `aln_n`, `aln_o`, `aln_e`,
`aln_i`, `aln_d`, `aln_l`, `aln_k`, `aln_m`, `aln_M`, `aln_O`, `aln_E`,
`aln_R`, `aln_q`, `aln_N`, `sampe_a`, `sampe_o`, `sampe_n`, `sampe_N`,
`sampe_c`, `sampe_s`, `samse_n`, `sw_a`, `sw_b`, `sw_q`, `sw_r`, `sw_w`,
`sw_m`, `sw_T`, `sw_c`, `sw_z`, `sw_s`, `sw_N`. Each one of this parameters
directly correspond to the respective command-line argument, e.g. `aln_o`
corresponds to the `-o` option of `bwa aln` (maximum number of gap opens).
Certain options, such as the `-t` option of `bwa aln`, are not exposed to
users because they are set by the app, based on the kind of cloud
environment that the app runs on.
