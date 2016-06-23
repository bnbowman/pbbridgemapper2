"""This is main modules for pbbridgemapper. The entry point is the
main function.

Pbbridgemapper finds split alignments for PacBio subreads. It does this
in three steps:

    1. Find unaligned subread affixes of sufficient and create a spoofed
    region table so that only the affixes are marked as inserts.
    2. Call pbalign using the new region table.
    3. Parse the cmp.h5 file produced by pbalign to create the file needed
    by SMRTview.

Pbbridgemapper takes three files as input:
    1. Input fofn. This is a fofn containing bax.h5 files. These are produced
    by SMRTpipe but easy enough independently.
    2. Alignment cmp.h5 file. This is the file produced by the first alignment
    of the data in the input fofn. It can have other alignments as well, and it
    doesn't need to be sorted. We're only looking at the /AlnInfo/AlnIndex
    dataset.
    3. Reference repository. Path to the reference repository directory. This
    would be created by referenceUploader, for example.

There is one general annoyance: subreads are identified by movie name, hole
number, and start and end with respect to the polymerase read. The start and
end values are not available in a cmp.h5; it only has the start and end of the
aligned portion of the subread. So we have to look back at the region table of
the input bax.h5 files.
"""


import errno
import logging
import os
import re
import subprocess

from pbcore.io import BasH5IO, CmpH5IO, FofnIO

import pbbridgemapper.affixes
import pbbridgemapper.smrtview_output
import pbbridgemapper.options


def create_affix_region_tables(input_fofn_filename, cmph5_filename,
                               output_path, min_affix_size):
    """Create the pbbridgemapper rgn.h5 and fofn files.

    Args:
        input_fofn_filename: fofn of bax.h5 filenames
        cmph5_filename: aligned_reads.cmp.h5 for the bax.h5 files
        output_path: path where the fofn and rgn.h5 files will be written
        min_affix_size: smallest affix that will be included in the
                        region table

    Returns:
        output_fofn_filename: file name of the FOFN of pbbridgemapper rgn.h5
                              files
    """

    bash5_filenames = list(FofnIO.readFofn(input_fofn_filename))
    logging.info("Read filenames from input fofn file: %s", bash5_filenames)

    cmph5_reader = CmpH5IO.CmpH5Reader(cmph5_filename)
    logging.info("Opened %s", cmph5_filename)
    output_rgn_filenames = []

    try:
        os.makedirs(os.path.join(output_path, 'pbbridgemapper_regions'))
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise OSError(
                "Could not create regions directory {d}."
                .format(d=os.path.join(output_path, 'pbbridgemapper_regions')))

    for bash5_filename in bash5_filenames:
        logging.debug("Getting affix boundaries from %s", bash5_filename)

        bash5_reader = BasH5IO.BasH5Reader(bash5_filename)
        subread_dict = pbbridgemapper.affixes.subread_dictionary(bash5_reader)
        logging.debug("Created subread dictionary from %d ZMWs",
                      len(subread_dict))

        affix_bounds = pbbridgemapper.affixes.affix_boundaries(
            subread_dict, cmph5_reader, min_affix_size)
        logging.debug("Found %d unmapped affixes",
                      sum([len(k) for k in affix_bounds.itervalues()]))

        original_region_table = (bash5_reader.file
                                 .get('/PulseData/Regions').value)
        affix_region_table = pbbridgemapper.affixes.affix_region_table(
            original_region_table, bash5_reader.movieName, affix_bounds)

        output_rgn_filename = os.path.join(
            output_path, 'pbbridgemapper_regions',
            re.sub(r"ba[sx]\.h5$", "rgn.h5", os.path.basename(bash5_filename)))
        pbbridgemapper.affixes.write_region_table(affix_region_table,
                                                  output_rgn_filename,
                                                  bash5_reader)
        bash5_reader.close()
        logging.info("Wrote pbbridgemapper region table to %s",
                     output_rgn_filename)
        output_rgn_filenames.append(output_rgn_filename)

    # Now the rgn.h5 files have been created, we just need to make the fofn
    output_fofn_filename = os.path.join(
        output_path, re.sub("fofn$", "pbbridgemapper_regions.fofn",
                            os.path.basename(input_fofn_filename)))

    with open(output_fofn_filename, 'w') as output_fofn_file:
        for i in xrange(len(output_rgn_filenames)):
            filename = output_rgn_filenames[i]
            output_fofn_file.write(filename)
            if i < len(output_rgn_filenames) - 1:
                output_fofn_file.write('\n')
    logging.info("Wrote rgn file names to %s", output_fofn_filename)

    return output_fofn_filename


def call_pbalign(input_fofn_filename, pbbridgemapper_fofn_filename,
                 reference_path, output_path, nproc):
    """Call pbalign using the orignal bax.h5 files and the modified region
    file.

    Args:
        input_fofn_filename: FOFN of the input bax.h5 files
        pbbridgemapper_rgn_fofn_filename: FOFN created by
            create_affix_region_tables
        reference_path: path to the reference repository
        output_path: directory where the cmp.h5 will be written
        nproc: number of procs to use for pbalign

    Returns:
        output_cmph5_filename: file name of the cmp.h5 containing the
            alignments of the affixes

    Raises:
        RuntimeError if pbalign does not return with exitcode 0
    """

    output_cmph5_filename = os.path.join(
        output_path, re.sub("fofn$", "pbbridgemapper.cmp.h5",
                            os.path.basename(input_fofn_filename)))

    pbalign_cmd = ['pbalign', input_fofn_filename, reference_path]
    pbalign_cmd.append(output_cmph5_filename)
    pbalign_cmd.extend(['--regionTable', pbbridgemapper_fofn_filename])
    pbalign_cmd.extend(['--nproc', str(nproc)])
    pbalign_cmd.append('--concordant')  # I guess? I don't know

    logging.info("Calling pbalign with command line %s",
                 ' '.join(pbalign_cmd))
    proc = subprocess.Popen(pbalign_cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    stdout, stderr = proc.communicate()
    logging.info("Finished running pbalign.")

    if proc.returncode != 0:
        logging.error("pbalign failed. Stderr was %s", stderr)
        raise RuntimeError("pbalign exited with returncode {e}"
                           .format(e=proc.returncode))

    return output_cmph5_filename


def create_pbbridgemapper_output(input_fofn_filename, affix_cmph5_filename,
                                 primary_cmph5_filename, split_reads_filename,
                                 unique_only):
    """Create the split_reads file the SMRTview wants.
    """

    # First build a dictionary of all the subreads in the input fofn
    bash5_filenames = list(FofnIO.readFofn(input_fofn_filename))
    all_subread_dict = {}
    for bash5_filename in bash5_filenames:
        bash5_reader = BasH5IO.BasH5Reader(bash5_filename)
        subread_dict = pbbridgemapper.affixes.subread_dictionary(bash5_reader)
        all_subread_dict.update(subread_dict)

    # Now iterate through the primary alignments, recording best scoring
    # primary alignments for each subread
    primary_cmph5_reader = CmpH5IO.CmpH5Reader(primary_cmph5_filename)
    for alignment in primary_cmph5_reader:
        key, overlapping_bounds = pbbridgemapper.affixes.find_subread_entry(
            alignment, all_subread_dict)
        if key is None:
            continue

        if all_subread_dict[key][overlapping_bounds] is not None:
            existing_alignment = (all_subread_dict[key][overlapping_bounds]
                                  ['primary'])
        else:
            existing_aligment = None

        if (existing_aligment is None or
                alignment.mapQV > existing_alignment['map_qv']):
            all_subread_dict[key][overlapping_bounds] = {}
            all_subread_dict[key][overlapping_bounds]['primary'] = (
                pbbridgemapper.smrtview_output.alignment_to_output_dict(
                    alignment, overlapping_bounds))

    primary_cmph5_reader.close()

    # Now iterate through the affix alignments
    try:
        affix_cmph5_reader = CmpH5IO.CmpH5Reader(affix_cmph5_filename)
    except CmpH5IO.EmptyCmpH5Error:
        affix_cmph5_reader = []

    for alignment in affix_cmph5_reader:
        key, overlapping_bounds = pbbridgemapper.affixes.find_subread_entry(
            alignment, all_subread_dict)
        if key is None:
            continue

        # Figure out if this is a prefix or suffix alignment
        if all_subread_dict[key][overlapping_bounds] is None:
            continue
        alignment_dict = (pbbridgemapper.smrtview_output
                          .alignment_to_output_dict(alignment,
                                                    overlapping_bounds))
        primary_dict = all_subread_dict[key][overlapping_bounds]['primary']
        primary_subread_start = primary_dict['subread_start']
        primary_subread_end = primary_dict['subread_end']

        if alignment_dict['subread_end'] <= primary_subread_start:
            affix_type = 'prefix'
        elif alignment_dict['subread_start'] >= primary_subread_end:
            affix_type = 'suffix'
        else:
            continue

        if affix_type in all_subread_dict[key][overlapping_bounds]:
            existing_alignment = (all_subread_dict[key][overlapping_bounds]
                                                  [affix_type])
        else:
            existing_alignment = None

        if (existing_alignment is None or
                alignment.mapQV > existing_alignment['map_qv']):
            all_subread_dict[key][overlapping_bounds][affix_type] = \
                alignment_dict

    if unique_only:
        pbbridgemapper.smrtview_output.remove_nonunique_alignments(
            all_subread_dict)
    pbbridgemapper.smrtview_output.write_split_reads_file(all_subread_dict,
                                                          split_reads_filename)


def setup_log(verbose):
    """Set up the logging settings."""

    if verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    log_format = ('[%(levelname)s] %(asctime)-15s [%(funcName)s '
                  '%(lineno)d] %(message)s')

    logging.basicConfig(level=log_level, format=log_format)


def main():
    """The entry point for pbbridgemapper."""

    parser = pbbridgemapper.options.parse_options()
    args = parser.parse_args()
    setup_log(args.verbose)

    rgn_fofn_filename = create_affix_region_tables(
        args.input_fofn, args.aligned_reads_cmp_h5,
        args.output_path, args.min_affix_size)

    affix_cmph5_filename = call_pbalign(
        args.input_fofn, rgn_fofn_filename,
        args.reference_path, args.output_path, args.nproc)

    if args.split_reads_file is None:
        split_reads_basename = re.sub("fofn$", "",
                                      os.path.basename(args.input_fofn))
        split_reads_basename += 'split_reads.bridgemapper.gz'
        split_reads_filename = os.path.join(args.output_path,
                                            split_reads_basename)
    else:
        split_reads_filename = args.split_reads_file

    create_pbbridgemapper_output(
        args.input_fofn, affix_cmph5_filename,
        args.aligned_reads_cmp_h5,
        split_reads_filename,
        args.unique_only)
