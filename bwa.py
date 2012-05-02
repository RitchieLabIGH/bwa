import dxpy
import subprocess, logging
from math import floor, ceil

logging.basicConfig(level=logging.DEBUG)

def main():
    t = dxpy.new_dxgtable([{"name": "leMapping", "type": "string"}])
    t.add_types(["LetterMappings", "Mappings"])
    t.add_types(["BwaLetterContigSetV1"])
    t.close(block=True)
    f = dxpy.upload_local_file("/bin/ls", wait_on_close=True)
    job['output'] = {'mappings': [dxpy.dxlink(t)], 'indexed_reference': dxpy.dxlink(f)}

def old_main():
    # Create a table for the mappings
    # See http://samtools.sourceforge.net/SAM1.pdf
    mappings_schema = [
        {"name": "QNAME", "type": "string"},
        {"name": "FLAG", "type": "int32"},
        {"name": "chr", "type": "string"}, # RNAME
        {"name": "start", "type": "int32"}, # POS
        {"name": "MAPQ", "type": "int32"},
        {"name": "CIGAR", "type": "string"},
        {"name": "RNEXT", "type": "string"},
        {"name": "PNEXT", "type": "int32"},
        {"name": "TLEN", "type": "int32"},
        {"name": "SEQ", "type": "string"},
        {"name": "QUAL", "type": "string"},
        {"name": "OPT", "type": "string"},
        {"name": "stop", "type": "int32"},
    ]

    tableId = dxpy.new_dxgtable(mappings_schema).get_id()
    # rowCount = dxpy.open_dxgtable(job['input']['leftReads']).describe()['numRows']
    rowCount = dxpy.open_dxgtable(job['input']['leftReads']).describe()['size']
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
            'leftReads': job['input']['leftReads']['$dnanexus_link'],
            'rightReads': job['input']['rightReads']['$dnanexus_link'],
            'bwaGenomeArchive': job['input']['bwaGenomeArchive']['$dnanexus_link'],
            'chromosomes': job['input']['chromosomes']['$dnanexus_link'],
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

def old_map():
    subprocess.check_call("dx_tableToFastq %s --start_row %d --end_row %d > left.fq" % (job['input']['leftReads'], job['input']['from'], job['input']['to']), shell=True)
    subprocess.check_call("dx_tableToFastq %s --start_row %d --end_row %d > right.fq" % (job['input']['rightReads'], job['input']['from'], job['input']['to']), shell=True)
    dxpy.download_dxfile(job['input']['bwaGenomeArchive'], 'genome.tgz')

    subprocess.check_call("tar zxf genome.tgz", shell=True)
    subprocess.check_call("bwa aln genome.fa left.fq >left.fq.sai", shell=True)
    subprocess.check_call("bwa aln genome.fa right.fq >right.fq.sai", shell=True)
    subprocess.check_call("bwa sampe genome.fa left.fq.sai right.fq.sai left.fq right.fq >out.sam", shell=True)
    
    subprocess.check_call("dx_storeSamAsMappingsTable out.sam " + job['input']['tableId'] + " " + job['input']['chromosomes'], shell=True)
    job['output']['id'] = job['input']['tableId']

def old_reduce():
    t = dxpy.open_dxgtable(job['input']['tableId'])
    t.close()
    job['output']['mappings'] = t.get_id()
