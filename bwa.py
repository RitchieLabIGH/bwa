import dxpy
import subprocess, logging, os, time, re
from multiprocessing import Pool, cpu_count

logging.getLogger().setLevel(logging.DEBUG)

def run_shell(command):
    logging.debug("Running "+command)
    subprocess.check_call(command, shell=True)

def make_indexed_reference(job_inputs):
    logging.info("Indexing reference genome")

    run_shell("contigset2fasta %s reference.fasta" % job_inputs['reference']['$dnanexus_link'])
    ref_details = dxpy.DXRecord(job_inputs['reference']['$dnanexus_link']).get_details()
    ref_name = dxpy.DXRecord(job_inputs['reference']['$dnanexus_link']).describe()['name']

    # TODO: test if the genomes near the boundary work OK
    if sum(ref_details['contigs']['sizes']) < 2*1024*1024*1024:
        subprocess.check_call("bwa index -a is reference.fasta", shell=True)
    else:
        subprocess.check_call("bwa index -a bwtsw reference.fasta", shell=True)

    subprocess.check_call("XZ_OPT=-0 tar -cJf reference.tar.xz reference.fasta*", shell=True)
    indexed_ref_dxfile = dxpy.upload_local_file("reference.tar.xz", hidden=True, wait_on_close=True)
    
    indexed_ref_record = dxpy.new_dxrecord(name=ref_name + " (indexed for BWA)",
                                           types=["BwaLetterContigSetV2"],
                                           details={'index_archive': dxpy.dxlink(indexed_ref_dxfile.get_id()),
                                                    'original_contigset': job_inputs['reference']})
    indexed_ref_record.close()
    
    # TODO: dxpy project workspace convenience functions
# FIXME
#    if "projectWorkspace" in job:
#        indexed_ref_record.clone(job["projectWorkspace"])
    
    return indexed_ref_record

@dxpy.entry_point('main')
def main(**job_inputs):
    job_outputs = {}
    reads_inputs = job_inputs['reads']
    reads_ids = [r['$dnanexus_link'] for r in reads_inputs]
    reads_descriptions = {r: dxpy.DXGTable(r).describe() for r in reads_ids}
    reads_columns = {r: [col['name'] for col in desc['columns']] for r, desc in reads_descriptions.items()}
    
    print reads_inputs
    print reads_ids
    print reads_descriptions
    print reads_columns

    all_reads_have_FlowReads_tag = all(['FlowReads' in desc['types'] for desc in reads_descriptions.values()])
    all_reads_have_LetterReads_tag = all(['LetterReads' in desc['types'] for desc in reads_descriptions.values()])
    reads_have_names = any(['name' in columns for columns in reads_columns.values()])
    reads_are_paired = any(['sequence2' in columns for columns in reads_columns.values()])
    reads_have_qualities = any(['quality' in columns for columns in reads_columns.values()])
    if reads_have_qualities:
        assert(all(['quality' in columns for columns in reads_columns.values()]))
    
    if job_inputs["algorithm"] == "bwasw":
        assert(not reads_are_paired) # bwasw does not support paired inputs
    
    assert(all_reads_have_FlowReads_tag or all_reads_have_LetterReads_tag)
    
    reference_record_types = dxpy.describe(job_inputs['reference'])['types']
    if "BwaLetterContigSetV2" in reference_record_types:
        input_ref_is_indexed = True
    elif "ContigSet" in reference_record_types:
        input_ref_is_indexed = False
    else:
        raise dxpy.ProgramError("Unrecognized object passed as reference. It must be a ContigSet record or a BwaLetterContigSetV2 file")
    
    if input_ref_is_indexed:
        job_outputs['indexed_reference'] = job_inputs['reference']
    else:
        found_cached_idx = False
        for result in dxpy.find_data_objects(classname='record',
                                             typename='BwaLetterContigSetV2',
                                             link=job_inputs['reference']['$dnanexus_link']):
            job_outputs['indexed_reference'] = dxpy.dxlink(result['id'])
            found_cached_idx = True
            break
        if not found_cached_idx:
            job_outputs['indexed_reference'] = dxpy.dxlink(make_indexed_reference(job_inputs))
    
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
                          ("template_id", "int64"),
                          ("read_group", "int32")])

    # optional sam fields: RG BC XC XT NM CM XN SM AM XM X0 X1 XG MD XA

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

    table_columns.extend([("sam_field_BC", "string"),
                          ("sam_field_XC", "int32"),
                          ("sam_field_XT", "string"),
                          ("sam_field_NM", "int32"),
                          ("sam_field_CM", "int32"),
                          ("sam_field_XN", "int32"),
                          ("sam_field_SM", "int32"),
                          ("sam_field_AM", "int32"),
                          ("sam_field_XM", "int32"),
                          ("sam_field_X0", "int32"),
                          ("sam_field_X1", "int32"),
                          ("sam_field_XG", "int32"),
                          ("sam_field_MD", "string"),
                          ("sam_field_XA", "string"),
                          ("sam_optional_fields", "string")])

    
    column_descriptors = [dxpy.DXGTable.make_column_desc(name, type) for name, type in table_columns]
    
    gri_index = dxpy.DXGTable.genomic_range_index("chr", "lo", "hi")
    t = dxpy.new_dxgtable(column_descriptors, indices=[gri_index])
    
    if input_ref_is_indexed:
        original_contigset = dxpy.get_details(job_inputs['reference'])['original_contigset']
    else:
        original_contigset = job_inputs['reference']
    t.set_details({'original_contigset': original_contigset})

    t.add_types(["LetterMappings", "Mappings", "gri"])

    # name table
    if 'output name' in job_inputs:
        t.rename( job_inputs['output name'] )
    else:
        first_reads_name = dxpy.DXGTable( job_inputs['reads'][0] ).describe()['name']
        contig_set_name = dxpy.describe(job_inputs['reference'])['name']
        # if we're working on an indexed_reference we're not guaranteed to have access to original_contigset
        if input_ref_is_indexed:
            contig_set_name = contig_set_name.split('(index')[0]
        t.rename(first_reads_name + " mapped to " + contig_set_name)

    # declare how many paired or single reads are in each reads table
    read_group_lengths = []
    for i in range(len(reads_ids)):
        current_length = reads_descriptions[reads_ids[i]]["length"] 
        if 'sequence2' in dxpy.DXGTable(reads_ids[i]).get_col_names():
            num_pairs = current_length
            num_singles = 0
        else:
            num_pairs = 0
            num_singles = current_length

        read_group_lengths.append( {"num_singles":num_singles, "num_pairs":num_pairs} )
    
    details = t.get_details()
    details['read_groups'] = read_group_lengths
    t.set_details(details)

    row_offsets = []; row_cursor = 0
    for i in range(len(reads_ids)):
        row_offsets.append(row_cursor)
        row_cursor += reads_descriptions[reads_ids[i]]["length"]

    chunk_size = job_inputs["chunk_size"]

    map_job_inputs = job_inputs.copy()
    map_job_inputs["row_offsets"] = row_offsets
    map_job_inputs["num_rows"] = chunk_size
    map_job_inputs["table_id"] = t.get_id()
    map_job_inputs["indexed_reference"] = job_outputs['indexed_reference']

    postprocess_job_inputs = job_inputs.copy()
    postprocess_job_inputs["table_id"] = t.get_id()
    
    for start_row in xrange(0, row_cursor, chunk_size):
        map_job_inputs["start_row"] = start_row
        map_job = dxpy.new_dxjob(map_job_inputs, "map")
        print "Launched map job with", map_job_inputs
        postprocess_job_inputs["chunk%dresult" % start_row] = {'job': map_job.get_id(), 'field': 'ok'}
        postprocess_job_inputs["chunk%ddebug" % start_row] = {'job': map_job.get_id(), 'field': 'debug'}
        
    postprocess_job = dxpy.new_dxjob(postprocess_job_inputs, "postprocess")

    job_outputs['mappings'] = {'job': postprocess_job.get_id(), 'field': 'mappings'}

    print "MAIN OUTPUT:", job_outputs
    return job_outputs

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

@dxpy.entry_point('map')
def map(**job_inputs):
    print "Map:", job_inputs
    job_outputs = {}
    times = [('start', time.time())]
    reads_inputs = job_inputs['reads']
    reads_ids = [r['$dnanexus_link'] for r in reads_inputs]
    reads_descriptions = {r: dxpy.DXGTable(r).describe() for r in reads_ids}
    reads_columns = {r: [col['name'] for col in desc['columns']] for r, desc in reads_descriptions.items()}

    reads_are_paired = any(['sequence2' in columns for columns in reads_columns.values()])

    times.append(('preamble', time.time()))

    dxpy.download_dxfile(dxpy.get_details(job_inputs["indexed_reference"])['index_archive'], "reference.tar.xz")
    times.append(('download reference', time.time()))

    # TODO: Async everything below
    # subprocess.check_call("pixz -d reference.tar.xz && tar -xf reference.tar", shell=True)
    subprocess.check_call("tar -xJf reference.tar.xz", shell=True)

    if job_inputs["algorithm"] == "bwasw":
        bwa_algorithm = "bwasw"
    else:
        # algorithm = aln or auto. TODO: check what auto should do
        bwa_algorithm = "aln"
    
    aln_opts, sampe_opts, sw_opts = parse_bwa_cmd_opts(job_inputs)
    
    # Set the number of threads BWA parameter to the apparent number of CPUs.
    aln_opts += " -t " + str(cpu_count())
    sw_opts += " -t " + str(cpu_count())
    
    row_offsets = job_inputs['row_offsets']   # starting row for each reads table if you added them all up
    start_row = job_inputs['start_row']       # the position in this chunk relative to the row_offsets 'total'
    num_rows = job_inputs['num_rows']         # size of chunk to do this time
    subjobs = []
    for i in range(len(reads_ids)):
        reads_length = reads_descriptions[reads_ids[i]]["length"]
        read_group = i
        # see if the reads table is part of this chunk

        # if start is inside this reads table, add it
        # doing this in the form:   (A_start < B_end) and (A_end > B_start)
        # A is the reads tables
        # B is the current chunk
        # A_start = row_offsets[i]
        # A_end = row_offsets[i] + reads_length
        # B_start = start_row
        # B_end = start_row + num_rows
        if row_offsets[i] < (start_row+num_rows) and (row_offsets[i]+reads_length) > start_row:

            rel_start = max(start_row - row_offsets[i], 0)
            rel_end = min(reads_length, start_row - row_offsets[i] + num_rows) # Using half-open intervals: [start, end)
            subjobs.append({'reads_id': reads_ids[i], 'start_row': rel_start, 'end_row': rel_end, 'read_group':read_group})
    
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
        cmd += " --table_id '%s'" % job_inputs["table_id"]
        cmd += " --reads_id '%s'" % reads_id
        cmd += " --start_row %d" % subjob['start_row']
        cmd += " --read_group %d" % subjob['read_group']
        
        if job_inputs.get('discard_unmapped_rows'):
            cmd += " --discard_unmapped_rows"
        run_shell(cmd)
        times.append(('run table upload (subchunk %d)' % subchunk_id, time.time()))

    job_outputs["ok"] = True
    
    timing_report = {}
    for i in range(len(times)-1):
        timing_report[times[i+1][0]] = times[i+1][1] - times[i][1]
    job_outputs["debug"] = {'times': timing_report}
    return job_outputs

@dxpy.entry_point('postprocess')
def postprocess(**job_inputs):
    print "Postprocess:", job_inputs
    job_outputs = {}

    time_report = {k: v for k, v in job_inputs.iteritems() if re.match("chunk\d+debug", k)}
    
    t = dxpy.DXGTable(job_inputs["table_id"])
    d = t.get_details()
    d['time_report'] = time_report
    t.set_details(d)
    t.close(block=True)
    job_outputs['mappings'] = dxpy.dxlink(t)
    return job_outputs

dxpy.run()
