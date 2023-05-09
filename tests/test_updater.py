import unittest
from unittest.mock import MagicMock, patch
from requests.exceptions import Timeout

from sierralocal.updater import update_apobec, update_hivdb


class TestUpdater(unittest.TestCase):

    def testUpdateApobec(self):
        exp_filepath = r'sierralocal\data\apobec_drms.json'

        with patch('sierralocal.updater.update_apobec'):
            res_filepath = update_apobec()

        self.assertIn(exp_filepath, res_filepath)

        # mock_response = MagicMock()
        # mock_requests.get.return_value = mock_response

    def testUpdateHivdb(self):
        exp_filepath = r'sierralocal\data\HIVDB'

        with patch('sierralocal.updater.update_hivdb', return_value = r'sierra-local\sierralocal\data\HIVDB_9.4.xml'):
            res_filepath = update_hivdb()

        self.assertIn(exp_filepath, res_filepath)
