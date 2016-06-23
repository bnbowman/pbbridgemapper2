import logging
import os

from pbcore.io import BasH5IO

import base_test_case
from pbbridgemapper import affixes

log = logging.getLogger(__name__)

class TestSubreadDictionary(base_test_case.BaseTestCase):

    @classmethod
    def setUpClass(cls):
        base_test_case.BaseTestCase.setUpClass()

        bash5_filename = os.path.join(base_test_case.ROOT_DATA_DIR,
                                      "m130522_092457_42208_cTEST1_s1_p0.1.bax.h5")
        cls.bash5_reader = BasH5IO.BasH5Reader(bash5_filename)
        cls.subread_dict = affixes.subread_dictionary(cls.bash5_reader)

    def test_read_all_zmws(self):
        """All zmws in the region table must be read into the dict"""

        keys = self.subread_dict.keys()
        sequencing_zmws = self.bash5_reader.sequencingZmws
        log.info("Number of zmw read: %d", len(keys))
        self.assertEqual(len(keys), len(sequencing_zmws))
        zmws = [k[1] for k in keys]
        for zmw in sequencing_zmws:
            self.assertIn(zmw, zmws)
        self.assertEqual(min(zmws), min(sequencing_zmws))
        self.assertEqual(max(zmws), max(sequencing_zmws))

    def test_hq_region_overlap(self):
        """No subread region should fall outside of an HQ region."""

        region_table = self.bash5_reader.file["PulseData/Regions"].value

        zmw_18_table = region_table[region_table[..., 0] == 18]
        hq_row = zmw_18_table[zmw_18_table[..., 1] == 2][0]
        hq_start = hq_row[2]
        hq_end = hq_row[3]
        
        movie_name = self.bash5_reader.movieName

        subread_bounds_list = self.subread_dict[(movie_name, 18)].keys()

        for subread_bounds in subread_bounds_list:
            self.assertLessEqual(hq_start, subread_bounds[0])
            self.assertGreaterEqual(hq_end, subread_bounds[1])
