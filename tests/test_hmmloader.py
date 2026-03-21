import unittest
import tempfile
from pathlib import Path

# on older python versions, we require importlib.resrouces now
try:
    from importlib.resources import files as resource_files
except ImportError:
    from importlib_resources import files as resource_files  # type: ignore

from cayman.annotate.crazy_annotator import HMMLoader, ThresholdTable
from . import test_data

class Test_CLI(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tempfile = tempfile.NamedTemporaryFile()
        cls.tempfile2 = tempfile.NamedTemporaryFile()
        cls.hmms_path = Path(resource_files(test_data).joinpath("hmms_bin", "cayman.v3.h3m"))  # ty:ignore[invalid-argument-type]
        cls.cutoff_file = Path(resource_files(test_data).joinpath("cutoffs.csv"))  # ty:ignore[invalid-argument-type]

    @classmethod
    def tearDownClass(cls):
        cls.tempfile.close()

    def test_hmm_selection(self):

        blacklist = ThresholdTable.load(self.cutoff_file).hmms_which_will_be_skipped

        hmms = HMMLoader.read_hmms(
            hmmdb_path=self.hmms_path,
            blacklist=blacklist,
            seed=42,
        )
        hmms.write_to_h3m_file(Path(self.tempfile.name))

        hash = HMMLoader.md5_of_file(Path(self.tempfile.name))

        self.assertEqual(hash, "d41d8cd98f00b204e9800998ecf8427e")