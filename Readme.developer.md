Letterspace BWA Aligner, Developer Readme
=========================================

Introduction
------------

This app runs the BWA aligner (letterspace only, version 0.6.2-r126).

The app defines a type called BwaLetterContigSetV3 which represents an
indexed genome in letterspace, appropriate for this version of the bwa. If
such an indexed genome is given in the input, it will be used as-is;
otherwise a (non-indexed) genome of type ContigSet must be given instead,
and it will be indexed.

The app receives an array of LetterReads objects, which may be
heterogeneous. These will be mapped and represented in a single output of
type Mappings.

Inputs
------

The app MUST receive in the "reads" parameter an array of LetterReads
objects. The app MUST receive the in the "reference" parameter a record of
type ContigSet or BwaLetterContigSetV3.

The objects in the "reads" parameter MUST be homogeneous in terms of the
presence of the quality column (either present in all or missing in all).

The app can be given the algorithm to run ("aln" or "bwasw"); it will
choose automatically if algorithm is set to "auto". If the app ends up
using the "aln" algorithm (either because that was set in the input, or
because that was chosen by the app in the case of "auto"), it will use any
options specified in the aln\_* and sampe\_* (for paired inputs)
parameters. If the app ends up using the "bwasw" algorithm, it will use any
options specified in the sw\_* parameters. The aln\_\*, sampe\_* and sw\_*
parameters reflect precisely the arguments that the "bwa aln", "bwa sampe"
and "bwa bwasw" commands can take.

This app parallelizes itself into using multiple mapping jobs. The
chunk\_size argument governs the maximum number of reads per mapping job
launched by the map. It is optional and should not need modification in
normal operation.

Outputs
-------

The app outputs an object of type LetterMappings containing the mappings,
and an object of type BwaLetterContigSetV3.

Implementation
--------------

This section provides a sketch of the implementation, which may serve as
inspiration for those writing their own alignment apps:

* Validate all inputs (Reads, ContigSet, BwaLetterContigSetV3, etc.). Ensure
  that the reads are homogeneous with respect to the presence of the quality
  field. Ensure that the reads are LetterReads, etc.

* If the input reference is a ContigSet, fetch the genome from the "reference"
  (as fasta) and index it with either "bwa index -a is" (up to 2Gbases genome)
  or "bwa index -a bwtsw" (&gt;2Gbases genome). Tar and xzip the
  results (including the fasta) and upload it back to the platform; output
  "indexed\_reference" as a record of type BwaLetterContigSetV3 with the following
  fields: "index\_archive" linking to the file above, and "original\_contigset"
  linking to the original ContigSet from which it was created.

* Determine the sets of columns that are required in the output -- for
  LetterReads, these are: sequence (always), name (if at least one Reads object
  has name), quality (if present), status, chr, lo, hi, negative\_strand,
  error\_probability, qc, cigar, template\_id. If at least one Reads object is
  paired: mate\_id, status2, chr2, lo2, hi2, negative\_strand2, proper\_pair.
  Create a gtable with a genomic range index.

* Align each of the objects in the "reads" input. Start by determining the
  number of rows of each Reads input, and calculate an "offset template id".
  The first reads object will be assigned an offset template id of 0, the
  second reads object will be assigned an offset template id equal to the
  number of rows in the first reads object, etc. For each Reads object R, do
  the following:

    * Decide on an appropriate chunk size to parallelize over. A good chunk
      might be 25M rows. For each chunk C, do the following:

        * Convert the particular piece (chunk C of reads object R) into one or
          two fasta or fastq files (say 1.fq and 2.fq), depending on whether
          the input is paired or not and whether it has quality scores or not.
          The conversion is as follows:

            * Original read names are ignored. Instead, for read names, use a
              number representing the row index in R that the entry corresponds
              to. Thus, read names are expected to be "1", "2", "3" etc.

            * Paired reads are named in the same way, but with "/1" and "/2"
              appended to the name for the first and second mate, respectively.

            * The sequence is taken from the "sequence" and "sequence2" fields

            * The quality is taken from the "quality" field. If there isn't
              such a field, output a fasta instead.

        * Fetch and extract the indexed\_reference (either from the input, or
          from the one generated at the earlier step)

        * If the "aln" algorithm is chosen, call "bwa aln" on 1.fq (and,
          separately, "bwa aln" on 2.fq, if paired), followed by "bwa samse"
          or "bwa sampe" (if paired). The result will be a sam file, sorted by
          read index.

        * If the "bwasw" algorithm is chosen, call "bwa bwasw" (on 1.fq, or on
          "1.fq 2.fq" if paired; either way it is run only once). The result
          will be a sam file, sorted by read index.

        * Start reading (streaming) in parallel from the sam file and from the
          chunk C of the Reads object R. These should correspond, as the sam
          file is expected to be sorted by read index and contain exactly the
          information as to how reads in C were mapped. For each row Q in C, do
          the following:

            * Consume all the consecutive rows S in the sam file that
              correspond to the row Q in C. It could be one or more rows
              depending on whether the row Q in C represents paired input or
              not. It should be easy to figure out which rows to consume from
              the sam file, because the read names contain the row indices of
              the rows in C.

            * Process the set S of rows, and the row Q, and push some
              respective info to the output mappings. More specifically, fill
              in the following values:

                * name, sequence, quality, flowgram, flow\_indices,
                  clip\_qual\_left, clip\_qual\_right, clip\_adapter\_left, and
                  clip\_adapter\_right are taken from Q; the rest are taken
                  from S
                * chr and lo are taken from RNAME and POS
                * hi is calculated by counting the number of M/D/N/=/X in CIGAR
                  and adding it to lo
                * status ::= (flags &amp; 0x4) ? UNMAPPED : (flags &amp; 0x100)
                  ? SECONDARY : PRIMARY
                * negative\_strand ::= (flags &amp; 0x10)
                * error\_probability is taken from MAPQ
                * qc is taken from (flags &amp; 0x200 and/or flags &amp; 0x400)
                * cigar is taken from CIGAR
                * mate\_id should be 0, 1, or -1 based on S/Q.
                * status2, chr2, lo2, hi2, negative\_strand and proper\_pair
                  should be reconciled from S.
                * template\_id should be "**offset template id**" plus the row
                  index of Q.

* Close and emit the output gtable.
