"""options parses the available options for pbbridgemapper."""
import argparse
import os

import pbbridgemapper
import pbbridgemapper.smrtview_output

def parse_options():
    """Return an argparse.ArgumentParser object with the args for
    pbbridgemapper.
    """

    def canonicalized_file_path(path):
        """Fix non-conforming path names. Copied from Quiver."""
        return os.path.abspath(os.path.expanduser(path))

    def check_gt_zero(value):
        int_value = int(value)
        if int_value <= 0:
            raise argparse.ArgumentTypeError(
                "%s must be an integer greater than zero." % value)
        return int_value

    desc = ("pbbridgemapper finds split alignments of PacBio reads and "
            "reports them in a way that allows for visualization in "
            "SMRTView.")
    
    class FormatterClass(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=FormatterClass)
    
    parser.add_argument(
        '--version',
        action='version',
        version=str(pbbridgemapper.get_full_version()))
    
    parser.add_argument(
        "--output_columns",
        action="version",
        version="\t".join(pbbridgemapper.smrtview_output.OUTPUT_COLUMNS),
        help=argparse.SUPPRESS)

    parser.add_argument(
        "input_fofn",
        action="store",
        type=canonicalized_file_path,
        help="FOFN of input bax.h5 files.")

    parser.add_argument(
        "aligned_reads_cmp_h5",
        action="store",
        type=canonicalized_file_path,
        help="aligned_reads.cmp.h5 file for the files in input_fofn.")

    parser.add_argument(
        "--reference_path", "-r",
        action="store",
        type=canonicalized_file_path,
        required=True,
        help="Path to the reference repository.")

    parser.add_argument(
        "--output_path", "-o",
        action="store",
        type=canonicalized_file_path,
        required=True,
        help="Directory where output files will be written.")

    parser.add_argument(
        "--split_reads_file",
        action="store",
        type=canonicalized_file_path,
        help=("Name of the output split_reads file. If unspecified, "
              "bridgemapper will replace the 'fofn' in the input_fofn with "
              "'split_reads.bridgemapper.gz"))
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        default=False,
        help="Report additional information in the log.")

    parser.add_argument(
        "--nproc",
        action="store",
        type=int,
        default=1,
        help="nproce to pass to pbalign for the affix alignment.")

    parser.add_argument(
        "--min_affix_size",
        action="store",
        type=check_gt_zero,
        default=50,
        help="Minimum length for a prefix or suffix to be reported.")

    parser.add_argument(
        "--unique_only",
        action="store_true",
        default=False,
        help=("Only report prefix and suffix alignments that have"
              "only one valid alignment to the reference."))

    return parser
