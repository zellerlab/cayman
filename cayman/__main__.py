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

    logger.info("Version: %s gqlib: %s", __version__, gqlib_version)
    logger.info(
        "Command: %s %s",
        os.path.basename(sys.argv[0]), " ".join(sys.argv[1:])
    )

    print(args)

    input_data = check_input_reads(
        args.reads1, args.reads2,
        args.singles, args.orphans,
    )

    if not os.path.exists(args.annotation_db):
        raise ValueError(
            f"{args.annotation_db} is not a valid annotation database"
        )

    if not check_bwa_index(args.bwa_index):
        raise ValueError(f"{args.bwa_index} is not a valid bwa index.")

    if os.path.dirname(args.out_prefix):
        pathlib.Path(os.path.dirname(args.out_prefix)).mkdir(
            exist_ok=True, parents=True
        )

    db_importer = SmallDatabaseImporter(
        logger, args.annotation_db, single_category="cazy", sep=args.db_separator, coords=args.db_coordinates,
    )
    logger.info("Finished loading database.")

    profiler = RegionQuantifier(
        db=db_importer,
        out_prefix=args.out_prefix,
        ambig_mode="1overN",
        reference_type="domain",
    )

    aln_runner = BwaMemRunner(
        args.cpus_for_alignment,
        args.bwa_index,
        sample_id=os.path.basename(args.out_prefix),
    )

    for input_type, *reads in input_data:

        logger.info("Running %s alignment: %s", input_type, ",".join(reads))
        proc, call = aln_runner.run(
            reads,
            single_end_reads=input_type == "single",
        )

        try:
            profiler.count_alignments(
                proc.stdout,
                aln_format="sam",
                min_identity=args.min_identity,
                min_seqlen=args.min_seqlen,
            )

        except Exception as err:
            if isinstance(err, ValueError) and str(err).strip() == "file does not contain alignment data":
                # pylint: disable=W1203
                logger.error(f"Failed to align. Is `{args.aligner}` installed and on the path?")
                logger.error("Aligner call was:")
                logger.error("%s", call)
                sys.exit(1)

            logger.error("Encountered problems digesting the alignment stream:")
            logger.error("%s", err)
            logger.error("Aligner call was:")
            logger.error("%s", call)
            logger.error("Shutting down.")
            sys.exit(1)

    profiler.finalise(restrict_reports=("raw", "rpkm",))


if __name__ == "__main__":
    main()
