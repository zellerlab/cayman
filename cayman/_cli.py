import logging
import os
import sys
from typing import Optional, List

# pylint: disable=W0611

from gqlib import __version__ as gqlib_version
from .handle_args import build_parser, set_log_lvl_from_args
from . import __version__


logger = logging.getLogger(__name__)

def main(argv: Optional[List[str]] = None):

    set_log_lvl_from_args(argv)

    parser = build_parser()
    args = parser.parse_args(args=argv)
    args.aligner = "bwa"

    logger.info("Version: %s gqlib: %s", __version__, gqlib_version)
    logger.info(
        "Command: %s %s",
        os.path.basename(sys.argv[0]), " ".join(sys.argv[1:])
    )

    return args.func(args)