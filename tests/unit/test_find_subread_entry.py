import logging
import os

from pbcore.io import BasH5IO, CmpH5IO

import base_test_case
from pbbridgemapper import affixes

log = logging.getLogger(__name__)

class TestFindSubreadEntry(base_test_case.BaseTestCase):

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
    
    def test_returns_full_subread_bounds(self):
        """Return rStart and rEnd from the baxh5 not the cmph5."""
        
        hole_status_dict = dict(k for k in
            zip(self.bash5_reader.file["PulseData/BaseCalls/ZMW/HoleNumber"],
                self.bash5_reader.file["PulseData/BaseCalls/ZMW/HoleStatus"]))
        for alignment in self.cmph5_reader:
            key, bounds = affixes.find_subread_entry(alignment,
                                                     self.subread_dict)
            
            # Blasr ignores holestatus but bash5_reader does not, so if we
            # don't find a hole, it must have a non-Sequencing status
            if key is None:
                self.assertNotEquals(hole_status_dict[alignment.HoleNumber], 0)
                continue

            self.assertEqual(key[0], self.bash5_reader.movieName)
            self.assertLessEqual(bounds[0], alignment.rStart)
            self.assertGreaterEqual(bounds[1], alignment.rEnd)

