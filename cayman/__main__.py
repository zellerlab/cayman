# pylint: disable=C0103,C0301,C0116

""" module docstring """

import logging
import os
import pathlib
import sys

# pylint: disable=W0611
from gqlib.profilers import RegionQuantifier
from gqlib.db.db_import import SmallDatabaseImporter
from gqlib.runners.alignment_runner import BwaMemRunner
from gqlib.ui.validation import check_bwa_index, check_input_reads

from . import __version__
from .handle_args import handle_args


logger = logging.getLogger(__name__)


# # pylint: disable=R0913
# def run_alignment(
#     profiler,
#     input_files,
#     bwa_index,
#     cpus_for_alignment=1,
#     min_identity=None,
#     min_seqlen=None,
#     unmarked_orphans=False,
# ):
#     """ docstring """
#     commands = [
#         f"bwa mem -v 1 -a -t {cpus_for_alignment} "
#         f"-K 10000000 {bwa_index} {' '.join(input_files)}",
#     ]

#     commands = " | ".join(commands)

#     logger.info("Used command: %s", commands)

#     try:
#         with subprocess.Popen(
#             commands, shell=True, stdout=subprocess.PIPE
#         ) as read_processing_proc:
#             profiler.count_alignments(
#                 read_processing_proc.stdout,
#                 aln_format="sam",
#                 min_identity=min_identity,
#                 min_seqlen=min_seqlen,
#                 unmarked_orphans=unmarked_orphans,
#             )
#     except Exception as err:
#         if isinstance(err, ValueError) and str(err).strip() == "file does not contain alignment data":
#             logger.error("Failed to align. Is `bwa mem` installed?")
#             sys.exit(1)
#         logger.error("Caught some exception:")
#         logger.error("%s", err)
#         raise Exception from err


def main():

    args = handle_args(sys.argv[1:])

    logger.info("Version: %s", __version__)
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
        stream = aln_runner.run(
            reads,
            single_end_reads=input_type == "single",
        )

        profiler.count_alignments(
            stream, aln_format="sam", min_identity=args.min_identity, min_seqlen=args.min_seqlen,
        )

        # run_alignment(
        #     profiler,
        #     reads,
        #     args.bwa_index,
        #     cpus_for_alignment=args.cpus_for_alignment,
        #     min_identity=args.min_identity,
        #     min_seqlen=args.min_seqlen,
        #     unmarked_orphans=input_type == "orphan",
        # )

    profiler.finalise(restrict_reports=("raw", "rpkm",))


if __name__ == "__main__":
    main()
