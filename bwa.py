import dxpy
import subprocess, logging, os, time, re
from multiprocessing import Pool, cpu_count

logging.getLogger().setLevel(logging.DEBUG)

def run_shell(command):
    logging.debug("Running "+command)
    subprocess.check_call(command, shell=True)

def make_indexed_reference():
    logging.info("Indexing reference genome")

    run_shell("contigset2fasta %s reference.fasta" % job['input']['reference']['$dnanexus_link'])
    ref_details = dxpy.DXRecord(job['input']['reference']['$dnanexus_link']).get_details()

    # TODO: test if the genomes near the boundary work OK
    if sum(ref_details['contigs']['sizes']) < 2*1024*1024*1024:
        subprocess.check_call("bwa index -a is reference.fasta", shell=True)
    else:
        subprocess.check_call("bwa index -a bwtsw reference.fasta", shell=True)

    subprocess.check_call("XZ_OPT=-1 tar -cJf reference.tar.xz reference.fasta*", shell=True)
    indexed_ref_dxfile = dxpy.upload_local_file("reference.tar.xz", keep_open=True)
    indexed_ref_dxfile.add_types(["BwaLetterContigSetV1"])
    indexed_ref_dxfile.set_details({'original_contigset': job['input']['reference']})
    indexed_ref_dxfile.close(block=True)
    return indexed_ref_dxfile

def main():
    reads_inputs = job['input']['reads']
    reads_ids = [r['$dnanexus_link'] for r in reads_inputs]
    reads_descriptions = {r: dxpy.DXGTable(r).describe() for r in reads_ids}
    reads_columns = {r: [col['name'] for col in desc['columns']] for r, desc in reads_descriptions.items()}
    
    all_reads_have_FlowReads_tag = all(['FlowReads' in desc['types'] for desc in reads_descriptions.values()])
    all_reads_have_LetterReads_tag = all(['LetterReads' in desc['types'] for desc in reads_descriptions.values()])
    reads_have_names = any(['name' in columns for columns in reads_columns.values()])
    reads_are_paired = any(['sequence2' in columns for columns in reads_columns.values()])
    reads_have_qualities = any(['quality' in columns for columns in reads_columns.values()])
    if reads_have_qualities:
        assert(all(['quality' in columns for columns in reads_columns.values()]))
    
    if job["input"]["algorithm"] == "bwasw":
        assert(not reads_are_paired) # bwasw does not support paired inputs
    
    assert(all_reads_have_FlowReads_tag or all_reads_have_LetterReads_tag)
    
    if 'indexed_reference' in job['input']:
        job['output']['indexed_reference'] = job['input']['indexed_reference']
    else:
        job['output']['indexed_reference'] = dxpy.dxlink(make_indexed_reference())
    
    table_columns = [("sequence", "string")]
    if reads_have_names:
        table_columns.append(("name", "string"))
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
        table_columns.extend([("mate_id", "int32"), # TODO: int8
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
    
    if 'indexed_reference' in job['input']:
        original_contigset = dxpy.DXFile(job['input']['indexed_reference']).get_details()['original_contigset']
    else:
        original_contigset = job['input']['reference']
    t.set_details({'original_contigset': original_contigset})

    t.add_types(["LetterMappings", "Mappings", "gri"])

    row_offsets = []; row_cursor = 0
    for i in range(len(reads_ids)):
        row_offsets.append(row_cursor)
        row_cursor += reads_descriptions[reads_ids[i]]["size"]
    
    chunk_size = job["input"]["chunk_size"]

    map_job_inputs = job["input"].copy()
    map_job_inputs["row_offsets"] = row_offsets
    map_job_inputs["num_rows"] = chunk_size
    map_job_inputs["table_id"] = t.get_id()
    map_job_inputs["indexed_reference"] = job['output']['indexed_reference']

    postprocess_job_inputs = job["input"].copy()
    postprocess_job_inputs["table_id"] = t.get_id()
    
    # Partition gtable part id range between jobs
    num_jobs = max(1, 1 + row_cursor/chunk_size)
    
    for start_row in xrange(0, row_cursor, chunk_size):
        map_job_inputs["start_row"] = start_row
        map_job = dxpy.new_dxjob(map_job_inputs, "map")
        print "Launched map job with", map_job_inputs
        postprocess_job_inputs["chunk%dresult" % start_row] = {'job': map_job.get_id(), 'field': 'ok'}
        postprocess_job_inputs["chunk%ddebug" % start_row] = {'job': map_job.get_id(), 'field': 'debug'}
        

    postprocess_job = dxpy.new_dxjob(postprocess_job_inputs, "postprocess")

    # (TODO: how do JBORs interact with array outputs?)
    job['output']['mappings'] = {'job': postprocess_job.get_id(), 'field': 'mappings'}

    print "MAIN OUTPUT:", job['output']

def write_reads_to_fastq(reads_id, filename, seq_col='sequence', qual_col='quality', start_row=0, end_row=None):
    row_id = start_row
    with open(filename, "w") as fh:
        for row in dxpy.open_dxgtable(reads_id).iterate_rows(columns=[seq_col, qual_col], start=start_row, end=end_row):
            for line in '@%d' % row_id, row[0], "+", row[1]:
                print >>fh, line
            row_id += 1

def write_reads_to_fasta(reads_id, filename, seq_col='sequence', start_row=0, end_row=None):
    row_id = start_row
    with open(filename, "w") as fh:
        for row in dxpy.open_dxgtable(reads_id).iterate_rows(columns=[seq_col], start=start_row, end=end_row):
            for line in '>%d' % row_id, row[0]:
                print >>fh, line
            row_id += 1

def run_alignment(algorithm, reads_file1, reads_file2=None, aln_opts='', sampe_opts='', sw_opts=''):
    commands = []
    if algorithm == "bwasw":
        if reads_file2 is None:
            commands.append("bwa bwasw reference.fasta {r1} {sw_opts} > {r1}.sam")
        else: # Paired read data
            commands.append("bwa bwasw reference.fasta {r1} {r2} {sw_opts} > {r1}.sam")
    else: # algorithm is "aln"
        commands.append("bwa aln reference.fasta {r1} {aln_opts} > {r1}.sai")
        if reads_file2 is not None:
            commands.append("bwa aln reference.fasta {r2} {aln_opts} > {r2}.sai")
            commands.append("bwa sampe reference.fasta {r1}.sai {r2}.sai {r1} {r2} {sampe_opts} > {r1}.sam")
        else:
            commands.append("bwa samse reference.fasta {r1}.sai {r1} > {r1}.sam")

    for command in commands:        
        run_shell(command.format(r1=reads_file1, r2=reads_file2, aln_opts=aln_opts, sampe_opts=sampe_opts, sw_opts=sw_opts))

def parse_bwa_cmd_opts(input):
    aln_opts, sampe_opts, sw_opts = '', '', ''
    for opt in 'n', 'o', 'e', 'i', 'd', 'l', 'k', 'm', 'M', 'O', 'E', 'R', 'q':
        if 'aln_'+opt in input:
            aln_opts += " -"+opt+" "+str(input['aln_'+opt])
    for opt in 'a', 'o', 'n', 'N', 'c':
        if 'sampe_'+opt in input:
            sampe_opts += " -"+opt+" "+str(input['sampe_'+opt])
    for opt in 'a', 'b', 'q', 'r', 'w', 'm', 'T', 'c', 'z', 's', 'N':
        if 'sw_'+opt in input:
            sw_opts += " -"+opt+" "+str(input['sw_'+opt])
    return aln_opts, sampe_opts, sw_opts

def map():
    print "Map:", job["input"]
    times = [('start', time.time())]
    reads_inputs = job['input']['reads']
    reads_ids = [r['$dnanexus_link'] for r in reads_inputs]
    reads_descriptions = {r: dxpy.DXGTable(r).describe() for r in reads_ids}
    reads_columns = {r: [col['name'] for col in desc['columns']] for r, desc in reads_descriptions.items()}

    reads_are_paired = any(['sequence2' in columns for columns in reads_columns.values()])

    times.append(('preamble', time.time()))
    dxpy.download_dxfile(job["input"]["indexed_reference"], "reference.tar.xz")
    times.append(('download reference', time.time()))

    # TODO: Async everything below
    # subprocess.check_call("pixz -d reference.tar.xz && tar -xf reference.tar", shell=True)
    subprocess.check_call("tar -xJf reference.tar.xz", shell=True)

    if job["input"]["algorithm"] == "bwasw":
        bwa_algorithm = "bwasw"
    else:
        # algorithm = aln or auto. TODO: check what auto should do
        bwa_algorithm = "aln"
    
    aln_opts, sampe_opts, sw_opts = parse_bwa_cmd_opts(job['input'])
    
    # Set the number of threads BWA parameter to the apparent number of CPUs.
    aln_opts += " -t " + str(cpu_count())
    sw_opts += " -t " + str(cpu_count())
    
    row_offsets = job['input']['row_offsets']
    start_row = job['input']['start_row']
    num_rows = job['input']['num_rows']
    subjobs = []
    for i in range(len(reads_ids)):
        reads_length = reads_descriptions[reads_ids[i]]["size"]
        if start_row >= row_offsets[i] and start_row < row_offsets[i] + reads_length:
            rel_start = start_row - row_offsets[i]
            rel_end = min(reads_length, start_row - row_offsets[i] + num_rows) # Using half-open intervals: [start, end)
            subjobs.append({'reads_id': reads_ids[i], 'start_row': rel_start, 'end_row': rel_end})
    
    times.append(('parse parameters', time.time()))
    print 'SUBJOBS:', subjobs
    
    for subchunk_id in range(len(subjobs)):
        subjob = subjobs[subchunk_id]
        reads_id = subjob['reads_id']
        # TODO: FlowReads trimming support
        if 'quality' in reads_columns[reads_id]:
            if reads_are_paired:
                reads_file1 = "input"+str(subchunk_id)+"_1.fastq"
                reads_file2 = "input"+str(subchunk_id)+"_2.fastq"
                write_reads_to_fastq(reads_id, reads_file1, seq_col='sequence', qual_col='quality', start_row=subjob['start_row'], end_row=subjob['end_row'])
                write_reads_to_fastq(reads_id, reads_file2, seq_col='sequence2', qual_col='quality2', start_row=subjob['start_row'], end_row=subjob['end_row'])
                times.append(('fetch reads (subchunk %d)' % subchunk_id, time.time()))
                run_alignment(bwa_algorithm, reads_file1, reads_file2, aln_opts=aln_opts, sampe_opts=sampe_opts, sw_opts=sw_opts)
                times.append(('run alignment (subchunk %d)' % subchunk_id, time.time()))
            else:
                reads_file1 = "input"+str(subchunk_id)+".fastq"
                write_reads_to_fastq(reads_id, reads_file1, start_row=subjob['start_row'], end_row=subjob['end_row'])
                run_alignment(bwa_algorithm, reads_file1, aln_opts=aln_opts, sampe_opts=sampe_opts, sw_opts=sw_opts)
        else: # No qualities, use plain fasta
            if reads_are_paired:
                reads_file1 = "input"+str(subchunk_id)+"_1.fasta"
                reads_file2 = "input"+str(subchunk_id)+"_2.fasta"
                write_reads_to_fasta(reads_id, reads_file1, seq_col='sequence', start_row=subjob['start_row'], end_row=subjob['end_row'])
                write_reads_to_fasta(reads_id, reads_file2, seq_col='sequence2', start_row=subjob['start_row'], end_row=subjob['end_row'])
                run_alignment(bwa_algorithm, reads_file1, reads_file2, aln_opts=aln_opts, sampe_opts=sampe_opts, sw_opts=sw_opts)
            else:
                reads_file1 = "input"+str(subchunk_id)+".fasta"
                write_reads_to_fasta(reads_id, reads_file1, start_row=subjob['start_row'], end_row=subjob['end_row'])
                run_alignment(bwa_algorithm, reads_file1, aln_opts=aln_opts, sampe_opts=sampe_opts, sw_opts=sw_opts)
        
        cmd = "dx_storeSamAsMappingsTable_bwa"
        cmd += " --alignments '%s.sam'" % reads_file1
        cmd += " --table_id '%s'" % job["input"]["table_id"]
        cmd += " --reads_id '%s'" % reads_id
        cmd += " --start_row %d" % subjob['start_row']
        
        if job['input'].get('discard_unmapped_rows'):
            cmd += " --discard_unmapped_rows"
        run_shell(cmd)
        times.append(('run table upload (subchunk %d)' % subchunk_id, time.time()))

    job["output"]["ok"] = True
    
    timing_report = {}
    for i in range(len(times)-1):
        timing_report[times[i+1][0]] = times[i+1][1] - times[i][1]
    job["output"]["debug"] = {'times': timing_report}

def postprocess():
    print "Postprocess:", job["input"]
    
    time_report = {k: v for k, v in job["input"].iteritems() if re.match("chunk\d+debug", k)}
    
    t = dxpy.DXGTable(job["input"]["table_id"])
    d = t.get_details()
    d['time_report'] = time_report
    t.set_details(d)
    t.close(block=True)
    job['output']['mappings'] = dxpy.dxlink(t)
