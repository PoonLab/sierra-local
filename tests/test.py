import unittest
# import stuff from project here

class DummyTest(unittest.TestCase):
    def setUp(self):
        # add unit test fixtures
        datum = 123.4
    def testSimple(self):
        expected = 'foobar'
        # do something that should yield expected output
        result = 'foobar'
        self.assertEqual(expected, result)
