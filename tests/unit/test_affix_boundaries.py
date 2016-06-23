import logging
import os

from pbcore.io import BasH5IO, CmpH5IO

import base_test_case
from pbbridgemapper import affixes

log = logging.getLogger(__name__)

class TestAffixBoundaries(base_test_case.BaseTestCase):

    @classmethod
    def setUpClass(cls):
        base_test_case.BaseTestCase.setUpClass()

        bash5_filename = os.path.join(base_test_case.ROOT_DATA_DIR,
                                      "m130522_092457_42208_cTEST1_s1_p0.1.bax.h5")
        cmph5_filename = os.path.join(base_test_case.ROOT_DATA_DIR,
                                      "test_alignment.cmp.h5")
        cls.bash5_reader = BasH5IO.BasH5Reader(bash5_filename)
        cls.cmph5_reader = CmpH5IO.CmpH5Reader(cmph5_filename)
        cls.subread_dict = affixes.subread_dictionary(cls.bash5_reader)
        cls.affix_bounds = affixes.affix_boundaries(cls.subread_dict,
                                                    cls.cmph5_reader, 1)

    def test_no_overlap_alignments(self):
        """Affixes cannot overlap with any aligned part of any subread."""
        
        alignment_dict = {}
        for alignment in self.cmph5_reader:
            key = (alignment.movieInfo.Name, alignment.HoleNumber)
            if key not in alignment_dict:
                alignment_dict[key] = []
            alignment_dict[key].append((alignment.rStart, alignment.rEnd))

        for key in self.affix_bounds:
            try:
                alignments = alignment_dict[key]
            except KeyError:
                pass

            affixes = self.affix_bounds[key]

            for affix in affixes:
                for alignment in alignments:
                    self.assertIsNone(BasH5IO.intersectRanges(affix, alignment))

    def test_contained_subreads(self):
        """Every affix must be contained by subread bounds from the bax.h5"""

        for key in self.affix_bounds:
            subread_bounds_list = self.subread_dict[key].keys()
            affixes = self.affix_bounds[key]

            for affix in affixes:
                found_containing_subread = False
                for subread_bounds in subread_bounds_list:
                    if subread_bounds[0] <= affix[0] < affix[1] <= subread_bounds[1]:
                        found_containing_subread = True
                        break
                self.assertTrue(found_containing_subread)

    def test_min_affix_size(self):
        """Affixes must be longer than min_affix_size."""
        affix_bounds = affixes.affix_boundaries(self.subread_dict,
                                                self.cmph5_reader, 100)

        for key in affix_bounds:
            for affix in affix_bounds[key]:
                self.assertGreaterEqual(affix[1]-affix[0], 100)

