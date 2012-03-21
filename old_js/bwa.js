// The BWA app
//
// To create this app, call the create app route with the contents of this file as
// the code, and the following specs:
//
// input_spec:
//
// {
//   leftReads: table,
//   rightReads: table,
//   bwaGenomeArchive: file
//   chromosomes: object
//   runtimeArgs: object
// }
//
// output_spec:
//
// {
//   mappings: coords_table
// }
//
//
// The two reads tables must be made from a fastq file (this app does not deal
// with the exact schema of the leftReads and rightReads tables as it expects
// to call a helper for converting to fastq)
//
// The bwaGenomeArchive must be a .tar.gz file which contains the following
// files:
//
// genome.fa
// genome.fa.amb
// genome.fa.ann
// genome.fa.bwt
// genome.fa.fai
// genome.fa.pac
// genome.fa.sa
//
// The chromosomes object must have a property called chromosomes that contains an
// array of "chromosomes" (where a "chromosome" is an array of three values:
// a chromosome name ("chr1"), index (0), and size (249250621)).
// See ../../push_hg19.js for an example.
//
// runtimeArgs: 
// {
//   rowFetchChunk: <split input reads into chunks of approx this many rows>, 
//   rowLimit: <process at most this many rows. 0 = process all table rows> 
// }

function main() {
  var runtimeArgs = {};
  if (job.input.runtimeArgs !== undefined) {
    runtimeArgs = DNAnexus.getObjectProperties(job.input.runtimeArgs);
  }

  // default = 15_000_000
  var rowFetchChunk = (runtimeArgs.rowFetchChunk !== undefined) ?  +runtimeArgs.rowFetchChunk : 15000000.0;
  var rowLimit = (runtimeArgs.rowLimit !== undefined) ?  +runtimeArgs.rowLimit : 0;

  // Create a table for the mappings

  // See http://samtools.sourceforge.net/SAM1.pdf
  var mappings_schema = [
    "QNAME:string",
    "FLAG:int32",
    "chr:int32", // RNAME
    "start:double", // POS
    "MAPQ:int32",
    "CIGAR:string",
    "RNEXT:string",
    "PNEXT:int32",
    "TLEN:int32",
    "SEQ:string",
    "QUAL:string",
    "OPT:string",
    "stop:double"
  ];

  var tableId = DNAnexus.createTable(mappings_schema).id;

  var rowCount = DNAnexus.getTableInfo(job.input.leftReads).numRows;
  if (rowLimit > 0 ) {
    rowCount = Math.min(rowCount,rowLimit);
  }
  var chunkCount = Math.ceil(rowCount / rowFetchChunk);
  var chunkSize = Math.floor(rowCount / chunkCount);

  // Iterate over chunks
  var i, from, to, mapInput, reduceInput = {}, mapJobId, reduceJobId;
  for (i = 0; i < rowCount; i += chunkSize) {
    from = i;
    to = Math.min(from + chunkSize - 1, rowCount - 1);
    mapInput = {
      leftReads: job.input.leftReads,
      rightReads: job.input.rightReads,
      bwaGenomeArchive: job.input.bwaGenomeArchive,
      chromosomes: job.input.chromosomes,
      from: from,
      to: to,
      tableId: tableId
    };
    // Run a "map" job for each chunk
    mapJobId = DNAnexus.runJob({input: mapInput, func_name: "map"});
    reduceInput["mapJob" + i.toString() + "TableId"] = {job: mapJobId, field: 'id'};
  }

  // Run a "reduce" job
  reduceInput.tableId = tableId;
  reduceJobId = DNAnexus.runJob({input: reduceInput, func_name: "reduce"});
  job.output = {mappings: {job: reduceJobId, field: 'mappings'}};
}

function map() {
  DNAnexus.system("dx_fetchReadsTableAsFastq " + job.input.leftReads.toString() + " " + job.input.from.toString() + " " + job.input.to.toString() + " >  left.fq");

  DNAnexus.system("dx_fetchReadsTableAsFastq " + job.input.rightReads.toString() + " " + job.input.from.toString() + " " + job.input.to.toString() + " > right.fq");

  DNAnexus.system("dx_fetchFile " + job.input.bwaGenomeArchive.toString() + " > genome.tgz");

  DNAnexus.system("tar zxf genome.tgz");
  DNAnexus.system("bwa aln genome.fa left.fq >left.fq.sai");
  DNAnexus.system("bwa aln genome.fa right.fq >right.fq.sai");
  DNAnexus.system("bwa sampe genome.fa left.fq.sai right.fq.sai left.fq right.fq >out.sam");

  DNAnexus.system("dx_storeSamAsMappingsTable out.sam " + job.input.tableId.toString() + " " + job.input.chromosomes.toString());

  job.output.id = job.input.tableId;
}

function reduce() {
  DNAnexus.closeTable(job.input.tableId);
  var id = DNAnexus.createCoordsTable("/tables/" + job.input.tableId).id;
  job.output.mappings = id;
}
