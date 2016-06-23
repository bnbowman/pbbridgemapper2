"""smrtview_output contains methods for writing an output file in the
Aformat that SMRTview expects. It's not super elegant.
"""

import gzip
from pbbridgemapper import get_full_version

OUTPUT_COLUMNS = [
    'read',                   # Read name, movieName/holeNumber/rStart_rEnd
    'subreadLength',          # Length of the subread
    'prologSubreadStart',     # Start of prolog alignment wrt the subread
    'prologSubreadEnd',       # End of prolog alignment wrt the subread
    'prologRefid',            # Ref string of prolog alignment, e.g. 'chr1'
    'prologTstrand',          # Strand of prolog alignment, 0 = Forward 1 = RC
    'prologStart',            # Start of prolog alignment wrt the target
    'prologEnd',              # End of prolog alignment wrt the target
    'prologScore',            # Prolog alignment score (from blasr)
    'prologPctsimilarity',    # Prolog alignment percent similarity
    'prologMapqv',            # Prolog map QV (from blasr)
    'primarySubreadStart',    # and so forth...
    'primarySubreadEnd',
    'primaryRefid',
    'primaryTstrand',
    'primaryStart',
    'primaryEnd',
    'primaryScore',
    'primaryPctsimilarity',
    'primaryMapqv',
    'epilogSubreadStart',
    'epilogSubreadEnd',
    'epilogRefid',
    'epilogTstrand',
    'epilogStart',
    'epilogEnd',
    'epilogScore',
    'epilogPctsimilarity',
    'epilogMapqv']


def output_line(subread_alns):
    """Create a line of the split_reads.bridgemapper file from a
    dict of subread alignments.
    """

    output = []
    output.append(subread_alns['primary']['read_name'])
    output.append(subread_alns['primary']['subread_length'])

    ordered_keys = ['subread_start', 'subread_end', 'ref_name',
                    'strand', 'target_start', 'target_end', 'score',
                    'pct_similarity', 'map_qv']

    for affix_type in ('prefix', 'primary', 'suffix'):
        if affix_type not in subread_alns:
            for i in xrange(len(ordered_keys)):
                output.append('_')
        else:
            for key in ordered_keys:
                output.append(subread_alns[affix_type][key])
    return '\t'.join([str(k) for k in output])

def count_subreads(subread_dict):
    """Count the subreads that will be written out so it can be put
    in the header.
    """
    counter = 0
    for key in subread_dict:
        for subread_bounds in subread_dict[key]:
            subread_alns = subread_dict[key][subread_bounds]
            if subread_alns is None:
                continue
            counter +=1
    return counter

def write_split_reads_file(subread_dict, split_reads_filename):
    """Write out the split_reads.bridgemapper file."""
    split_reads_file = gzip.open(split_reads_filename, 'w')
    
    num_total_subreads = count_subreads(subread_dict)
    split_reads_file.write("##pbbridgemapper version " + get_full_version())
    split_reads_file.write('\n')
    split_reads_file.write("##Total Subreads: " + str(num_total_subreads))
    split_reads_file.write('\n')
    split_reads_file.write('##' + '\t'.join(OUTPUT_COLUMNS))
    split_reads_file.write('\n')

    for key in subread_dict:
        for subread_bounds in subread_dict[key]:
            subread_alns = subread_dict[key][subread_bounds]
            if subread_alns is None:
                continue
            split_reads_file.write(output_line(subread_alns))
            split_reads_file.write('\n')


def remove_nonunique_alignments(subread_dict):
    """Pop and prefix or suffix alignment with MapQV < 254."""

    for key in subread_dict:
        for bounds in subread_dict[key]:
            affix_dict = subread_dict[key][bounds]
            for affix_type in affix_dict:
                if affix_type in ("prefix", "suffix"):
                    if affix_dict[affix_type]['map_qv'] < 254:
                        affix_dict.pop(affix_type)


def subread_dictionary(bash5_reader):
    """Return a dictionary of dictionaries. First level keys are
    (movie_name, hole_number). They point to a second dictionary with
    keys (rStart, rEnd) where rStart and rEnd are from the bax.h5, not
    the cmp.h5.

    """
    subread_dict = {}
    for zmw in bash5_reader:
        first_key = (bash5_reader.movieName, zmw.holeNumber)
        second_keys = [tuple(k) for k in zmw.insertRegions]

        if first_key not in subread_dict:
            subread_dict[first_key] = {}

        for second_key in second_keys:
            subread_dict[first_key][second_key] = {}

    return subread_dict


def alignment_to_output_dict(cmph5_alignment, subread_bounds):
    """Create a dictionary of features relevant to bridgemapper out
    of a pbcore.io.CmpH5IO.CmpH5Alignment object.
    """
    output_dict = {}
    read_name = (cmph5_alignment.movieInfo.Name + '/' +
                 str(cmph5_alignment.HoleNumber))
    read_name += '/' + str(subread_bounds[0]) + '_' + str(subread_bounds[1])
    output_dict['read_name'] = read_name
    output_dict['subread_length'] = subread_bounds[1] - subread_bounds[0]
    output_dict['subread_start'] = cmph5_alignment.rStart - subread_bounds[0]
    output_dict['subread_end'] = cmph5_alignment.rEnd - subread_bounds[0]
    output_dict['ref_name'] = cmph5_alignment.referenceName
    output_dict['strand'] = cmph5_alignment.RCRefStrand
    output_dict['target_start'] = cmph5_alignment.tStart
    output_dict['target_end'] = cmph5_alignment.tEnd
    output_dict['score'] = -100000
    output_dict['pct_similarity'] = round(cmph5_alignment.similarity, 4)
    output_dict['map_qv'] = cmph5_alignment.MapQV

    return output_dict
