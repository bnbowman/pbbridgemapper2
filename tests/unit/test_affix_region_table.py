import logging
import numpy
import os

from pbcore.io import BasH5IO, CmpH5IO

import base_test_case
from pbbridgemapper import affixes

log = logging.getLogger(__name__)

class TestAffixRegionTable(base_test_case.BaseTestCase):

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
        cls.original_region_table = cls.bash5_reader.file['PulseData/Regions'].value
        cls.region_table = affixes.affix_region_table(
            cls.original_region_table, cls.bash5_reader.movieName,
            cls.affix_bounds)
    
    def test_entire_read_covered(self):
        """There should be entries in the region table for each entire
        read.
        """

        zmw_lengths = {}
        for hole_number in numpy.unique(self.original_region_table[..., 0]):
            hole_indices = self.original_region_table[..., 0] == hole_number
            hole_rows = self.original_region_table[hole_indices]
            zmw_length = numpy.max(hole_rows[..., 3])
            zmw_lengths[hole_number] = zmw_length
        
        for hole_number in numpy.unique(self.region_table[..., 0]):
            hole_indices = self.original_region_table[..., 0] == hole_number
            hole_rows = self.original_region_table[hole_indices]
            non_hq_rows = hole_rows[hole_rows[..., 1] != 2]
            
            bounds = non_hq_rows[..., (2,3)].tolist()
            bounds.sort()
            previous_end = 0
            for bound in bounds:
                self.assertEqual(bound[0], previous_end)
                previous_end = bound[1]
            self.assertEqual(previous_end, zmw_lengths[hole_number])
            zmw_lengths.pop(hole_number)

        self.assertEqual(len(zmw_lengths), 0)

    def test_hq_region_identical(self):
        """The HQ region shoud not change between the original and the produced
        region tables unless there are no affixes to map.
        """
        
        for hole_number in numpy.unique(self.original_region_table[..., 0]):
            original_hole_indices = self.original_region_table[..., 0] == hole_number
            hole_indices = self.region_table[..., 0] == hole_number

            original_hole_rows = self.original_region_table[original_hole_indices]
            hole_rows = self.region_table[hole_indices]
            
            original_hq_row = original_hole_rows[original_hole_rows[..., 1] == 2][0]
            hq_row = hole_rows[hole_rows[..., 1] == 2][0]
            
            non_hq_rows = hole_rows[hole_rows[..., 1] != 2]
            
            if len(non_hq_rows) == 1 and non_hq_rows[0][1] == BasH5IO.ADAPTER_REGION:
                continue
            
            for i in range(len(hq_row)):
                self.assertEqual(original_hq_row[i], hq_row[i])

    def test_alternate_insert_adapter(self):
        """Insert and adapter regions should alternate."""
        
        for hole_number in numpy.unique(self.region_table[..., 0]):
            hole_indices = self.original_region_table[..., 0] == hole_number
            hole_rows = self.original_region_table[hole_indices]
            non_hq_rows = hole_rows[hole_rows[..., 1] != 2]
            
            type_bounds = non_hq_rows[..., (1,2,3)].tolist()
            type_bounds.sort(key=lambda x: x[1:])
            
            last_type = None
            for bound in type_bounds:
                bound_type = bound[0]
                if last_type is not None:
                    self.assertNotEqual(bound_type, last_type)
                last_type = bound_type
