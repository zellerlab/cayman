# pylint: disable=C0103,C0301,C0116

""" module docstring """

import logging
import os
import pathlib
import subprocess
import sys

# pylint: disable=W0611
from gqlib.profilers import RegionQuantifier
from gqlib.db.db_import import DomainBedDatabaseImporter

from . import __version__
from .handle_args import handle_args


logger = logging.getLogger(__name__)


def check_bwa_index(prefix):
    """ docstring """
    suffixes = (".amb", ".ann", ".bwt", ".pac", ".sa")
    return all(os.path.isfile(prefix + suffix) for suffix in suffixes)


# pylint: disable=R0913
def run_alignment(
    profiler,
    input_files,
    bwa_index,
    cpus_for_alignment=1,
    min_identity=None,
    min_seqlen=None,
    unmarked_orphans=False,
):
    """ docstring """
    commands = [
        f"bwa mem -v 1 -a -t {cpus_for_alignment} "
        f"-K 10000000 {bwa_index} {' '.join(input_files)}",
    ]

    commands = " | ".join(commands)

    logger.info("Used command: %s", commands)

    try:
        with subprocess.Popen(
            commands, shell=True, stdout=subprocess.PIPE
        ) as read_processing_proc:
            profiler.count_alignments(
                read_processing_proc.stdout,
                aln_format="sam",
                min_identity=min_identity,
                min_seqlen=min_seqlen,
                unmarked_orphans=unmarked_orphans,
            )
    except Exception as err:
        if isinstance(err, ValueError) and str(err).strip() == "file does not contain alignment data":
            logger.error("Failed to align. Is `bwa mem` installed?")
            sys.exit(1)
        logger.error("Caught some exception:")
        logger.error("%s", err)
        raise Exception from err


def check_input_reads(fwd=None, rev=None, singles=None, orphans=None):
    """ docstring """
    fwd_reads = fwd.split(",") if fwd else None
    rev_reads = rev.split(",") if rev else None
    single_reads = singles.split(",") if singles else None
    orphan_reads = orphans.split(",") if orphans else None

    all_readsets = []

    if fwd_reads and rev_reads:
        if len(fwd_reads) == len(rev_reads):
            all_readsets += zip(
                (["paired"] * len(fwd_reads)),
                fwd_reads, rev_reads
            )
        else:
            raise ValueError(
                f"Found different numbers of forward/R1 {len(fwd_reads)} "
                f"and reverse/R2 {len(rev_reads)} reads."
            )
    elif fwd_reads:
        logger.warning(
            "Found -1 forward/R1 reads but no -2 reverse/R2 reads. "
            "Treating these as single-end reads."
        )
        all_readsets += zip((["single"] * len(fwd_reads)), fwd_reads)
    elif rev_reads:
        raise ValueError(
            "Found -2 reverse/R2 reads but no -1 forward/R1 reads."
        )

    if single_reads:
        all_readsets += zip((["single"] * len(single_reads)), single_reads)
    if orphan_reads:
        all_readsets += zip((["orphan"] * len(orphan_reads)), orphan_reads)

    if not all_readsets:
        raise ValueError("No input reads specified.")

    for _, *reads in all_readsets:
        for r in reads:
            if not os.path.isfile(r):
                raise ValueError(f"{r} does not seem to be a valid read file.")

    return all_readsets


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

    db_importer = DomainBedDatabaseImporter(
        logger, args.annotation_db, single_category="cazy"
    )
    logger.info("Finished loading database.")

    profiler = RegionQuantifier(
        db=db_importer,
        out_prefix=args.out_prefix,
        ambig_mode="1overN",
        reference_type="domain",
    )

    for input_type, *reads in input_data:

        logger.info("Running %s alignment: %s", input_type, ",".join(reads))

        run_alignment(
            profiler,
            reads,
            args.bwa_index,
            cpus_for_alignment=args.cpus_for_alignment,
            min_identity=args.min_identity,
            min_seqlen=args.min_seqlen,
            unmarked_orphans=input_type == "orphan",
        )

    profiler.finalise(restrict_reports=("raw", "rpkm",))


if __name__ == "__main__":
    main()
