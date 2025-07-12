import os
import unittest
from collections import namedtuple

import pandas as pd

from pmutt import cantera


class TestCantera(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)

    def test_get_omkm_range(self):
        # Test that CTI range handles strings that can be directly converted to
        # int
        objs = ['r_0001', 'r_0002', 'r_0003', 'r_0005', 'r_0008',
                'r_0009', 'r_0010']
        CTI_range_out = ('["r_0001", "r_0002", "r_0003", "r_0005", '
                         '"r_0008", "r_0009", "r_0010"]')
        self.assertEqual(cantera._get_omkm_range(objs=objs), CTI_range_out)

        # Test that CTI range handles objects with ID as an attribute
        # that can be directly converted to int
        ObjWithId = namedtuple('ObjWithId', 'id')
        objs_id = [ObjWithId(obj_id) for obj_id in objs]
        self.assertEqual(cantera._get_omkm_range(objs=objs_id), CTI_range_out)

        # Test that CTI range handles objects with name as an attribute
        # that can be directly converted to int
        ObjWithName = namedtuple('ObjWithName', 'name')
        objs_name = [ObjWithName(obj_name) for obj_name in objs]
        self.assertEqual(cantera._get_omkm_range(objs=objs_name),
                         CTI_range_out)

        # Test that CTI range handles strings that have a header and a footer
        objs = [
            'test_0001', 'test_0002', 'test_0003', 'test_0005', 'test_0008',
            'test_0009', 'test_0010'
        ]
        CTI_range_out = ('["test_0001", "test_0002", "test_0003", '
                         '"test_0005", "test_0008", '
                         '"test_0009", "test_0010"]')
        self.assertEqual(cantera._get_omkm_range(objs=objs), CTI_range_out)

        # Test that CTI range handles multiple groups with headers and footers
        objs = [
            'test1_0001', 'test1_0002', 'test1_0003', 'test1_0005',
            'test1_0008', 'test1_0009', 'test1_0010', 'test2_0005',
            'test2_0006', 'test2_0100', '0001', '0002', '0003', '0005'
        ]
        CTI_range_out = ('["test1_0001", "test1_0002", "test1_0003", '
                         '"test1_0005", "test1_0008", "test1_0009", '
                         '"test1_0010", "test2_0005", "test2_0006", '
                         '"test2_0100", "0001", "0002", "0003", "0005"]')
        self.assertEqual(cantera._get_omkm_range(objs=objs), CTI_range_out)


if __name__ == '__main__':
    unittest.main()
