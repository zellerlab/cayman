# pylint: disable=C0103,C0301

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
    suffixes = (".amb", ".ann", ".bwt", ".pac", ".sa")
    return all(os.path.isfile(prefix + suffix) for suffix in suffixes)


def main():

    args = handle_args(sys.argv[1:])

    logger.info("Version: %s", __version__)
    logger.info("Command: %s %s", os.path.basename(sys.argv[0]), " ".join(sys.argv[1:]))

    print(args)
    if args.input_files != "-" and not all(os.path.exists(f) for f in args.input_files):
        input_files_str = "\n".join(args.input_files)
        raise ValueError(f"There is an issue with your input files. Please check.\n{input_files_str}")
    if not os.path.exists(args.annotation_db):
        raise ValueError(f"{args.annotation_db} is not a valid annotation database")
    if not check_bwa_index(args.bwa_index):
        raise ValueError(f"{args.bwa_index} is not a valid bwa index.")

    if os.path.dirname(args.out_prefix):
        pathlib.Path(os.path.dirname(args.out_prefix)).mkdir(
            exist_ok=True, parents=True
        )

    db_importer = DomainBedDatabaseImporter(logger, args.annotation_db)
    logger.info("Finished loading database.")

    fq = RegionQuantifier(
        db=db_importer,
        out_prefix=args.out_prefix,
        ambig_mode="1overN",
        unmarked_orphans=args.unmarked_orphans,
        reference_type="domain",
    )

    samtools_io_flags = "-buSh" if args.no_prefilter else "-Sh"

    commands = [
        f"bwa mem -a -t {args.cpus_for_alignment} -K 10000000 {args.bwa_index} {' '.join(args.input_files)}",
        f"samtools view -F 4 {samtools_io_flags} -",
    ]

    if not args.no_prefilter:
        logging.info("Prefiltering activated.")
        commands += [
            f"read_count {args.out_prefix}",
            "samtools view -buSh -",
            f"bedtools intersect -u -ubam -a stdin -b {args.annotation_db}",
        ]

    logger.info("Used command: %s", " | ".join(commands))

    try:
        with subprocess.Popen(" | ".join(commands), shell=True, stdout=subprocess.PIPE) as read_processing_proc:
            fq.process_bamfile(
                read_processing_proc.stdout,
                aln_format="bam",
                min_identity=args.min_identity, min_seqlen=args.min_seqlen,
                external_readcounts=None if args.no_prefilter else (args.out_prefix + ".readcount.json"),
            )

    except Exception as err:
        logger.error("Caught some exception:")
        logger.error("%s", err)
        raise Exception from err


if __name__ == "__main__":
    main()
