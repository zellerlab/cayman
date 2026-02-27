import unittest
import unittest.mock
import errno
import io
import os
import sys
import tempfile
import shutil
import logging
import pickle
from pathlib import Path

try:
    from importlib.resources import files as resource_files
except ImportError:
    from importlib_resources import files as resource_files  # type: ignore

import pyhmmer.plan7

import cayman.annotate.crazy_annotator
from cayman._cli import main, logger
from . import test_data


class Test_CLI(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # logger.handlers.clear()
        # logger.addHandler(logging.StreamHandler(stream=io.StringIO()))
        logging.disable()

        cls._orig_get_hits = cayman.annotate.crazy_annotator.CazyAnnotator.get_hits

        cls.tempfile = tempfile.NamedTemporaryFile(suffix=".csv")
        hmms_dir = str(
            Path(resource_files(test_data).joinpath("hmms")).resolve()
        )
        proteins_fasta = str(
            Path(resource_files(test_data).joinpath("protein.faa")).resolve()
        )
        cutoff_file = str(
            Path(resource_files(test_data).joinpath("cutoffs.csv")).resolve()
        )

        cls.arguments_normal = [
            "annotate_proteome",
            hmms_dir,
            proteins_fasta,
            "--cutoffs",
            cutoff_file,
            "-o",
            cls.tempfile.name,
            "-t",
            str(max(len(os.sched_getaffinity(0))-2, 1))
            if sys.platform == "linux"
            else "1"
        ]

        # cls.arguments_normal_with_pdb = [
        #     "-i",
        #     molecule_path,
        #     "-t",
        #     selected_template_dir,
        #     "-o",
        #     cls.tempfile.name,
        #     "--pdbs",
        #     str(Path(tempfile.gettempdir(), "pdbs").resolve()),
        #     "--include-template",
        # ]

        # cls.bad_argument_6 = [
        #     "-l",
        #     list_path4,
        #     "-o",
        #     cls.tempfile.name,
        # ]  # passing something which is not a PDB in the list

        cls.maxDiff = None

    @classmethod
    def tearDownClass(cls):
        cls.tempfile.close()
        # shutil.rmtree(Path(tempfile.gettempdir(), "pdbs"))

    def _search_hmm_mock(self, hmm, sequences, background):
        testdir = Path("tests/test_data/test_hits")
        testdir.mkdir(parents=True, exist_ok=True)
        fn = Path(testdir, f"{hmm.name}.pkl")
        if not fn.exists():
            hits = self._orig_get_hits(hmm, sequences, background)
            with open(fn, "wb") as f:
                pickle.dump(hits, f)
        with open(fn, "rb") as f:
            return pickle.load(f)

        #if we have the data somewhere, load it and return it
        # print("GEtting hits for", hmm)

        # pipeline = pyhmmer.plan7.Pipeline(background.alphabet, background)
        # return pipeline.search_hmm(hmm, sequences)

    def test_writing_pdbs_query_ref(self):
        with unittest.mock.patch(target="cayman.annotate.crazy_annotator.CazyAnnotator.get_hits", new=self._search_hmm_mock):
            self.assertEqual(main(self.arguments_normal), 0)

        with open(
            self.tempfile.name, "r"
        ) as f:
            actual = f.read()

        with resource_files(test_data).joinpath(
            "result.csv"
        ).open() as f:
            expected = f.read()

        self.assertMultiLineEqual(actual.strip(), expected.strip())