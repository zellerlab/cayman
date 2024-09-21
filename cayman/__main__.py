# pylint: disable=C0103,C0301,C0116

""" module docstring """

import logging
import os
import pathlib
import sys

# pylint: disable=W0611
from gqlib.db.db_import import SmallDatabaseImporter
from gqlib.profilers import RegionQuantifier
from gqlib.runners.alignment_runner import BwaMemRunner
from gqlib.ui.validation import check_bwa_index, check_input_reads

from gqlib import __version__ as gqlib_version
from .handle_args import handle_args
from . import __version__


logger = logging.getLogger(__name__)


def main():

    args = handle_args(sys.argv[1:])
    args.aligner = "bwa"

    logger.info("Version: %s gqlib: %s", __version__, gqlib_version)
    logger.info(
        "Command: %s %s",
        os.path.basename(sys.argv[0]), " ".join(sys.argv[1:])
    )

    print(args)
    args.func(args)

    return None

    


if __name__ == "__main__":
    main()
