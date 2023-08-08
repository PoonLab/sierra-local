import unittest
import os
from pathlib import Path
from unittest.mock import MagicMock, patch
from requests.exceptions import Timeout

from sierralocal.updater import *


class TestUpdater(unittest.TestCase):

    def testUpdateApobecMutation(self):
        exp_filepath = r'sierralocal/data/apobec_drms.json'

        with patch('sierralocal.updater.update_apobec'):
            res_filepath = update_apobec_mutation()
        res_filepath = '/'.join(res_filepath.split('/')[-3:])

        self.assertIn(exp_filepath, res_filepath)

        # mock_response = MagicMock()
        # mock_requests.get.return_value = mock_response

    def testUpdateSDRMs(self):
        exp_filepath = r'sierralocal/data/sdrms_hiv1.csv'
        with patch('sierralocal.updater.update_sdrms'):
            res_filepath = update_sdrms()
        res_filepath = '/'.join(res_filepath.split('/')[-3:])

        self.assertIn(exp_filepath, res_filepath)

    def testUpdateApobec(self):
        exp_filepath = r'sierralocal/data/apobecs.csv'
        with patch('sierralocal.updater.update_apobec'):
            res_filepath = update_apobec()
        res_filepath = '/'.join(res_filepath.split('/')[-3:])

        self.assertIn(exp_filepath, res_filepath)

    def testUpdateMutationType(self):
        exp_filepath = r'sierralocal/data/mutation-type-pairs_hiv1.csv'
        with patch('sierralocal.updater.update_mutation_type'):
            res_filepath = update_mutation_type()
        res_filepath = '/'.join(res_filepath.split('/')[-3:])

        self.assertIn(exp_filepath, res_filepath)

    def testUpdateIsUnusual(self):
        exp_filepath = r'sierralocal/data/rx-all_subtype-all.csv'
        with patch('sierralocal.updater.update_mutation_type'):
            res_filepath = update_is_unusual()
        res_filepath = '/'.join(res_filepath.split('/')[-3:])

        self.assertIn(exp_filepath, res_filepath)

    def testUpdateHivdb(self):
        with patch('sierralocal.updater.update_hivdb', return_value = r'sierra-local/sierralocal/data/HIVDB_9.4.xml'):
            res_filepath = update_hivdb()
        hivdb_file = res_filepath.split('/')[-1]
        exp_filepath = r'sierralocal/data/' + hivdb_file
        res_filepath = '/'.join(res_filepath.split('/')[-3:])

        self.assertIn(exp_filepath, res_filepath)
