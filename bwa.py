import dxpy
import subprocess, logging, os
from math import floor, ceil

logging.basicConfig(level=logging.DEBUG)

def make_indexed_reference():
    logging.info("Indexing reference genome")
    # TODO: get implementation of fasta fetcher from Joe
    # This is a debug stub only
    ref_details = dxpy.DXRecord(job['input']['reference']['$dnanexus_link']).get_details()
    seq_name = ref_details['contigs']['names'][0]
    seq_file_id = ref_details['flat_sequence_file']['$dnanexus_link']
    with open("reference.fasta", "wb") as fh:
        fh.write(">"+seq_name+"\n")
    dxpy.download_dxfile(seq_file_id, "reference.fasta", append=True)

    # TODO: test if the genomes near the boundary work OK
    if sum(ref_details['contigs']['sizes']) < 2*1024*1024*1024:
        subprocess.check_call("bwa index -a is reference.fasta", shell=True)
    else:
        subprocess.check_call("bwa index -a bwtsw reference.fasta", shell=True)

    subprocess.check_call("tar -cJf reference.tar.xz reference.fasta*", shell=True)
    indexed_ref_dxfile = dxpy.upload_local_file("reference.tar.xz", keep_open=True)
    indexed_ref_dxfile.add_types(["BwaLetterContigSetV1"])
    indexed_ref_dxfile.set_details({'originalContigSet': job['input']['reference']})
    indexed_ref_dxfile.close()
    return indexed_ref_dxfile

def main():
    reads_inputs = job['input']['reads']
    reads_ids = [r['$dnanexus_link'] for r in reads_inputs]
    reads_descriptions = {r: dxpy.DXGTable(r).describe() for r in reads_ids}
    reads_columns = {r: [col['name'] for col in desc['columns']] for r, desc in reads_descriptions.items()}
    
    all_reads_have_FlowReads_tag = all(['FlowReads' in desc['types'] for desc in reads_descriptions.values()])
    all_reads_have_LetterReads_tag = all(['LetterReads' in desc['types'] for desc in reads_descriptions.values()])
    reads_have_names = any(['name' in columns for columns in reads_columns.values()])
    reads_have_qualities = any(['quality' in columns for columns in reads_columns.values()])
    reads_are_paired = any(['sequence2' in columns for columns in reads_columns.values()])
    
    assert(all_reads_have_FlowReads_tag or all_reads_have_LetterReads_tag)
    
    if 'indexed_reference' in job['input']:
        job['output']['indexed_reference'] = job['input']['indexed_reference']
    else:
        job['output']['indexed_reference'] = dxpy.dxlink(make_indexed_reference())
    
    table_columns = [("sequence", "string")]
    if reads_have_names:
        table_columns.append(("names", "string"))
    if reads_have_qualities:
        table_columns.append(("quality", "string"))
    table_columns.extend([("status", "string"),
                          ("chr", "string"),
                          ("lo", "int32"),
                          ("hi", "int32"),
                          ("negative_strand", "boolean"),
                          ("error_probability", "uint8"),
                          ("qc", "string"),
                          ("cigar", "string"),
                          ("template_id", "int64")])
    if reads_are_paired:
        table_columns.extend([("mate_id", "string"),
                              ("status2", "string"),
                              ("chr2", "string"),
                              ("lo2", "int32"),
                              ("hi2", "int32"),
                              ("negative_strand2", "boolean"),
                              ("proper_pair", "boolean")])
    if all_reads_have_FlowReads_tag:
        table_columns.extend([("flowgram", "string"),
                              ("flow_indices", "string"),
                              ("clip_qual_left", "int32"),
                              ("clip_qual_right", "int32"),
                              ("clip_adapter_left", "int32"),
                              ("clip_adapter_right", "int32")])
    
    column_descriptors = [dxpy.DXGTable.make_column_desc(name, type) for name, type in table_columns]
    
    gri_index = dxpy.DXGTable.genomic_range_index("chr", "lo", "hi")
    t = dxpy.new_dxgtable(column_descriptors, indices=[gri_index])
    t.add_types(["LetterMappings", "Mappings"])

    row_offsets = {}; row_cursor = 0
    for reads_id, description in reads_descriptions.items():
        row_offsets[reads_id] = row_cursor
        row_cursor += description["size"]

    # FIXME: chunking

    chunk_size = 25000000
    
    

    map_job_inputs = job["input"].copy()
    map_job_inputs["start_row"] = 0
    map_job_inputs["num_rows"] = chunk_size
    map_job_inputs["table_id"] = t.get_id()
    map_job_inputs["indexed_reference"] = job['output']['indexed_reference']
    
    map_job = dxpy.new_dxjob(map_job_inputs, "map")

    postprocess_job_inputs = job["input"].copy()
    postprocess_job_inputs["table_id"] = t.get_id()
    postprocess_job_inputs["chunk1result"] = {'job': map_job.get_id(), 'field': 'ok'}

    postprocess_job = dxpy.new_dxjob(postprocess_job_inputs, "postprocess")

    # Scan through reads, computing start row offset
    # Save start row offsets in a record
    # Spawn map jobs, parameters: passthrough all params plus start_row, num_rows

    # In map job:
    # Convert with dx_tableToFastq or dx_tableToFasta depending on quality column presence
    #   - Rewrite to support trimming of FlowReads
    # Map with params; use the apparent number of cpus in bwa aln
    # 

    # (TODO: how do JBORs interact with array outputs?)
    job['output']['mappings'] = {'job': postprocess_job.get_id(), 'field': 'mappings'}

    print "MAIN OUTPUT:", job['output']

def write_reads_to_fastq(reads_id, filename, seq_col='sequence', qual_col='quality', start_row=0, end_row=None):
    row_id = start_row
    with open(filename, "w") as fh:
        for row in dxpy.open_dxgtable(reads_id).iterate_rows(columns=[seq_col, qual_col], start=start_row, end=end_row):
            fh.write("\n".join(['>'+str(row_id), row[0], "+", row[1], ""]))
            row_id += 1

def write_reads_to_fasta(reads_id, filename, seq_col='sequence', start_row=0, end_row=None):
    row_id = start_row
    with open(filename, "w") as fh:
        for row in dxpy.open_dxgtable(reads_id).iterate_rows(columns=[seq_col], start=start_row, end=end_row):
            fh.write("\n".join(['>'+str(row_id), row[0], ""]))
            row_id += 1

def run_shell(command):
    logging.debug("Running "+command)
    subprocess.check_call(command, shell=True)

def run_alignment(algorithm, filename1, filename2=None):
    if algorithm == "bwasw":
        if filename2 is not None:
            run_shell("bwa bwasw reference.fasta %s %s > %s.sai" % (filename1, filename2, filename1))
            run_shell("bwa sampe reference.fasta %s.sai %s.sai %s %s > %s.sam" % (filename1, filename2, filename1, filename2, filename1))
        else:
            run_shell("bwa bwasw reference.fasta %s > %s.sai" % (filename1, filename1))
            run_shell("bwa samse reference.fasta %s.sai %s > %s.sam" % (filename1, filename1))
    else:
        run_shell("bwa aln reference.fasta %s > %s.sai" % (filename1, filename1))
        if filename2 is not None:
            run_shell("bwa aln reference.fasta %s > %s.sai" % (filename2, filename2))
            run_shell("bwa sampe reference.fasta %s.sai %s.sai %s %s > %s.sam" % (filename1, filename2, filename1, filename2, filename1))
        else:
            run_shell("bwa samse reference.fasta %s.sai %s > %s.sam" % (filename1, filename1))

def map():
    print "Map:", job["input"]
    reads_inputs = job['input']['reads']
    reads_ids = [r['$dnanexus_link'] for r in reads_inputs]
    reads_descriptions = {r: dxpy.DXGTable(r).describe() for r in reads_ids}
    reads_columns = {r: [col['name'] for col in desc['columns']] for r, desc in reads_descriptions.items()}

    reads_are_paired = any(['sequence2' in columns for columns in reads_columns.values()])

    dxpy.download_dxfile(job["input"]["indexed_reference"], "reference.tar.xz")

    # TODO: Async everything below
    subprocess.check_call("tar -xJf reference.tar.xz", shell=True)

    if job["input"]["algorithm"] == "bwasw":
        bwa_algorithm = "bwasw"
    else:
        # algorithm = aln or auto. TODO: check what auto should do
        bwa_algorithm = "aln"
    
    subchunk_id = 0
    for reads_id in reads_ids:
        subchunk_id += 1
        if 'quality' in reads_columns[reads_id]:
            if reads_are_paired:
                filename1 = "input"+str(subchunk_id)+"_1.fastq"
                filename2 = "input"+str(subchunk_id)+"_2.fastq"
                write_reads_to_fastq(reads_id, filename1, seq_col='sequence', qual_col='quality')
                write_reads_to_fastq(reads_id, filename2, seq_col='sequence2', qual_col='quality2')
                run_alignment(bwa_algorithm, filename1, filename2)
            else:
                filename1 = "input"+str(subchunk_id)+".fastq"
                write_reads_to_fastq(reads_id, filename1)
                run_alignment(bwa_algorithm, filename1)
        else:
            if reads_are_paired:
                filename1 = "input"+str(subchunk_id)+"_1.fasta"
                filename2 = "input"+str(subchunk_id)+"_2.fasta"
                write_reads_to_fasta(reads_id, filename1, seq_col='sequence')
                write_reads_to_fasta(reads_id, filename2, seq_col='sequence2')
                run_alignment(bwa_algorithm, filename1, filename2)
            else:
                filename1 = "input"+str(subchunk_id)+".fasta"
                write_reads_to_fasta(reads_id, filename1)
                run_alignment(bwa_algorithm, filename1)
    
        cmd = "dx_storeSamAsMappingsTable_bwa"
        cmd += " --alignments '%s.sam'" % filename1
        cmd += " --table_id '%s'" % job["input"]["table_id"]
        cmd += " --reads_id '%s'" % reads_id
        cmd += " --start_row '%d'" % 0
        if False:
            cmd += " --end_row '%d'" % 0
        
        logging.debug("Would run "+cmd)

    job["output"]["ok"] = True

def postprocess():
    print "Postprocess:", job["input"]
    t = dxpy.DXGTable(job["input"]["table_id"])
    t.close(block=True)
    job['output']['mappings'] = [dxpy.dxlink(t)]

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
