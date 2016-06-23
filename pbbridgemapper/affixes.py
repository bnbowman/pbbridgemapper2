"""affixes contains methods for identifying portions of subreads
that are unmapped after the first resequencing pass.
"""
import h5py
import numpy

from pbcore.io import BasH5IO


def find_subread_entry(cmph5_alignment, subread_dict):
    """Find the entry in subread_dict that corresponds to the
    cmph5_alignment. Match the movie name and hole number, then find
    which of the bax.h5 subread bounds overlap with the bounds of the
    alignment.

    Args:
        cmph5_alignment: a CmpH5IO.CmpH5Alignment
        subread_dict: dictionary from affixes.subread_dictionary that will be
                      searched for a match with the cmph5_alignment

    Returns:
        key, overlapping_bounds: the (movie_name, hole_number) key and the
                                 subread bounds that contain the alignment

    Returns None if the alignment is not found in the subread_dict
    """

    key = (cmph5_alignment.movieInfo.Name, cmph5_alignment.HoleNumber)

    if key not in subread_dict:
        # If this alignment isn't in the input fofn, just skip it
        return None, None

    # Figure out which subread this alignment is from by checking for
    # overlap with the bounds from the region table
    overlapping_bounds = None
    for subread_bounds in subread_dict[key].iterkeys():
        if BasH5IO.intersectRanges((cmph5_alignment.rStart,
                                    cmph5_alignment.rEnd),
                                   subread_bounds) is not None:
            overlapping_bounds = subread_bounds
            break

    return key, overlapping_bounds


def affix_boundaries(subread_dict, cmph5_reader, min_affix_size):
    """Identifies the boundaries of unmapped subread affixes.

    Args:
        subread_dict: a dictionary keyed by (movie_name, hole_number)
                        and values of dicts keyed subread bounds from the
                        bax.h5
        cmph5_reader: a CmpH5IO.CmpH5Reader to the aligned_reads.cmph5
        min_affix_size: the minimum size for an affix to be reported

    Returns:
        affix_bounds: a dictionary keyed by (movie_name, hole_number)
                      and values of lists of bounds of affixes
    """

    alignment_extents = {}
    for key in subread_dict:
        alignment_extents[key] = {}
        for subread_bounds in subread_dict[key].iterkeys():
            alignment_extents[key][subread_bounds] = None

    for alignment in cmph5_reader:
        key, overlapping_bounds = find_subread_entry(alignment, subread_dict)
        if key is None:
            continue

        current_alignment_extent = alignment_extents[key][overlapping_bounds]
        if current_alignment_extent is None:
            alignment_extents[key][overlapping_bounds] = (alignment.rStart,
                                                          alignment.rEnd)
        else:
            alignment_extents[key][overlapping_bounds] = (
                min(alignment.rStart, current_alignment_extent[0]),
                max(alignment.rEnd, current_alignment_extent[1]))

    # Alright, we now know the extent of the primary alignments for each
    # subread from the input fofn. Now, record the boundaries of the
    # sufficiently large affixes
    affix_bounds = {}
    for key in subread_dict:
        for subread_bounds in subread_dict[key]:
            alignment_extent = alignment_extents[key][subread_bounds]
            if alignment_extent is None:
                # Skip reads that were entirely unmapped.
                continue

            prefix_bounds = (subread_bounds[0], alignment_extent[0])
            suffix_bounds = (alignment_extent[1], subread_bounds[1])

            if prefix_bounds[1] - prefix_bounds[0] >= min_affix_size:
                affix_bounds.setdefault(key, []).append(prefix_bounds)
            if suffix_bounds[1] - suffix_bounds[0] >= min_affix_size:
                affix_bounds.setdefault(key, []).append(suffix_bounds)

    return affix_bounds


def affix_region_table(original_region_table, movie_name, affix_bounds):
    """Creates a new region table so only the affixes will be mapped.
    Replaces the mapped portions of the original region table with adapter
    entries.

    Args:
        original_region_table: the /PulseData/Regions table from the original
                               bax.h5 file
        movie_name: the movie_name of the original_region_table. Used for
                    lookups in the affix_bounds
        affix_bounds: dictionary with keys (movie_name, hole_number) and values
                      of lists of (affix_start, affix_end)

    Returns:
        region_table: a numpy.array of where subread region that aren't the
                      affixes are marked as adapters
    """

    hole_numbers = numpy.unique(original_region_table[..., 0])
    output_region_table = []

    for hole_number in hole_numbers:
        hole_rows = original_region_table[
            original_region_table[..., 0] == hole_number]
        zmw_length = numpy.max(hole_rows[..., 3])
        hq_row = hole_rows[hole_rows[..., 1] == 2][0]

        # See if this hole has affix to be mapped. If not write an empty hq
        # region and a big adapter region and continue
        try:
            affixes = affix_bounds[(movie_name, hole_number)]
        except KeyError:
            output_region_table.append([hole_number, BasH5IO.HQ_REGION,
                                        0, 0, 0])
            output_region_table.append([hole_number, BasH5IO.ADAPTER_REGION,
                                        0, zmw_length, 0])
            continue

        affixes.sort()
        previous_bound = 0
        for i in xrange(len(affixes)):
            output_region_table.append([hole_number, BasH5IO.ADAPTER_REGION,
                                        previous_bound, affixes[i][0], 900])
            output_region_table.append([hole_number, BasH5IO.INSERT_REGION,
                                        affixes[i][0], affixes[i][1], -1])
            previous_bound = affixes[i][1]
        output_region_table.append([hole_number, BasH5IO.ADAPTER_REGION,
                                    previous_bound, zmw_length, 900])
        output_region_table.append(hq_row.tolist())

    output_region_array = numpy.array(output_region_table,
                                      dtype=original_region_table.dtype)

    return output_region_array


def subread_dictionary(bash5_reader):
    """Return a dictionary of dictionaries. First level keys are
    (movie_name, hole_number). They point to a second dictionary with
    keys (rStart, rEnd) where rStart and rEnd are from the region table
    in the bax.h5

    Args:
        bash5_reader: a BasH5IO.BasH5Reader object from which we will read
                      the region table

    Returns:
        subread_dict: a dictionary, as described above.
    """
    subread_dict = {}
    for zmw in bash5_reader:
        first_key = (bash5_reader.movieName, zmw.holeNumber)
        second_keys = [tuple(k) for k in zmw.insertRegions]

        if first_key not in subread_dict:
            subread_dict[first_key] = {}

        for second_key in second_keys:
            subread_dict[first_key][second_key] = None

    return subread_dict


def write_region_table(region_table, region_filename, bash5_reader):
    """Writes the region_table to region_filename.

    Args:
        region_table: a numpy.array of a region table
        region_filename: name of the file that will be created
        bash5_reader: BasH5IO.BasH5Reader that contains some fields that
                      blasr needs that we'll copy over
    """

    region_file = h5py.File(region_filename, 'w')
    pulsedata_group = region_file.create_group("PulseData")
    pulsedata_group.create_dataset("Regions", data=region_table)

    region_file.copy(bash5_reader.file['ScanData'], region_file['/'])

    for attr_item in bash5_reader.file['PulseData/Regions'].attrs.items():
        attr_key, attr_value = attr_item
        if isinstance(attr_value, basestring):
            new_dtype = h5py.special_dtype(vlen=str)
        elif attr_value.dtype == 'object':
            new_dtype = h5py.special_dtype(vlen=str)
        else:
            new_dtype = attr_value.dtype

        region_file["PulseData/Regions"].attrs.create(
            attr_key, attr_value, dtype=new_dtype)

    region_file.close()
