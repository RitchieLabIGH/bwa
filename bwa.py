'''

TODO: update docs for python port

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
'''

import dxpy
import subprocess
from math import floor, ceil

def main():
    # Create a table for the mappings
    # See http://samtools.sourceforge.net/SAM1.pdf
    mappings_schema = [
        "QNAME:string",
        "FLAG:int32",
        "chr:int32", # RNAME
        "start:double", # POS
        "MAPQ:int32",
        "CIGAR:string",
        "RNEXT:string",
        "PNEXT:int32",
        "TLEN:int32",
        "SEQ:string",
        "QUAL:string",
        "OPT:string",
        "stop:double"
    ]

    tableId = dxpy.new_dxtable(mappings_schema).get_id()
    # rowCount = dxpy.open_dxtable(job['input']['leftReads']).describe()['numRows']
    rowCount = dxpy.open_dxtable(job['input']['leftReads']).describe()['size']
    rowFetchChunk = job['input']['rowFetchChunk']
    rowLimit = job['input']['rowLimit']
    if (rowLimit > 0):
        rowCount = min(rowCount, rowLimit)
    chunkCount = int(ceil(float(rowCount) / float(rowFetchChunk)))
    chunkSize = int(floor(float(rowCount) / float(chunkCount)))

    reduceInput = {}

    # Iterate over chunks
    #var i, from, to, mapInput, reduceInput = {}, mapJobId, reduceJobId;
    for i in range(0, rowCount, chunkSize):
        _from = i
        to = min(_from + chunkSize - 1, rowCount - 1)
        mapInput = {
            'leftReads': job['input']['leftReads'],
            'rightReads': job['input']['rightReads'],
            'bwaGenomeArchive': job['input']['bwaGenomeArchive'],
            'chromosomes': job['input']['chromosomes'],
            'from': _from,
            'to': to,
            'tableId': tableId
        }
        # Run a "map" job for each chunk
        mapJobId = dxpy.new_dxjob(fn_input=mapInput, fn_name="map").get_id()
        reduceInput["mapJob" + str(i) + "TableId"] = {'job': mapJobId, 'field': 'id'}

    # Run a "reduce" job
    reduceInput['tableId'] = tableId
    reduceJobId = dxpy.new_dxjob(fn_input=reduceInput, fn_name="reduce").get_id()
    job['output'] = {'mappings': {'job': reduceJobId, 'field': 'mappings'}}

def map():
    subprocess.check_call("dx_tableToFastq %s --start_row %d --end_row %d > left.fq" % (job['input']['leftReads'], job['input']['from'], job['input']['to']), shell=True)
    subprocess.check_call("dx_tableToFastq %s --start_row %d --end_row %d > right.fq" % (job['input']['rightReads'], job['input']['from'], job['input']['to']), shell=True)
    dxpy.download_dxfile(job['input']['bwaGenomeArchive'], 'genome.tgz')

    subprocess.check_call("tar zxf genome.tgz", shell=True)
    subprocess.check_call("bwa aln genome.fa left.fq >left.fq.sai", shell=True)
    subprocess.check_call("bwa aln genome.fa right.fq >right.fq.sai", shell=True)
    subprocess.check_call("bwa sampe genome.fa left.fq.sai right.fq.sai left.fq right.fq >out.sam", shell=True)
    
    subprocess.check_call("dx_storeSamAsMappingsTable out.sam " + job['input']['tableId'] + " " + job['input']['chromosomes'], shell=True)
    job['output']['id'] = job['input']['tableId']

def reduce():
    t = dxpy.open_dxtable(job['input']['tableId'])
    t.close()
    job['output']['mappings'] = t.get_id()
