# pylint: disable=C0301,C0103
""" module docstring """

import argparse
import logging
import textwrap

from . import __version__


def handle_args(args):

    log_ap = argparse.ArgumentParser(prog="cayman", add_help=False)
    log_ap.add_argument("-l", "--log_level", type=int, choices=range(1, 5), default=logging.INFO)
    log_args, _ = log_ap.parse_known_args(args)

    try:
        logging.basicConfig(
            level=log_args.log_level,
            format='[%(asctime)s] %(message)s'
        )
    except ValueError as invalid_loglevel_err:
        raise ValueError(f"Invalid log level: {log_args.log_level}") from invalid_loglevel_err

    ap = argparse.ArgumentParser(
        prog="cayman",
        formatter_class=argparse.RawTextHelpFormatter,
        parents=(log_ap,),
    )
    ap.add_argument(
        "annotation_db",
        type=str,
        help=textwrap.dedent(
            """\
            Path to an sqlite3 database containing the reference annotation.
			"""
        ),
    )
    ap.add_argument(
        "bwa_index",
        type=str, help="",
    )
    ap.add_argument(
        "input_files",
        type=str,
        nargs="*",
        help=textwrap.dedent(
            """\
            Path to metagenomic reads in fastq format.
            Fastq files can be supplied as a single unpaired file or two paired-end files.
            Input from STDIN can be used with '-'."""
        ),
    )
    ap.add_argument(
        "--out_prefix",
        "-o",
        type=str,
        default="cayman",
        help="Prefix for output files.",
    )
    
    ap.add_argument(
        "--min_identity",
        type=float,
        default=0.97,
        help="Minimum sequence identity [n_match/length] for an alignment to be considered.",
    )

    ap.add_argument(
        "--min_seqlen",
        type=int,
        default=45,
        help="Minimum read length [bp] for an alignment to be considered.",
    )

    # orphan reads will not have flag 0x1 set
    ap.add_argument(
        "--unmarked_orphans",
        action="store_true",
        help="Ensure that alignments from unmarked orphan reads (from preprocessing) are properly accounted for.",
    )

    ap.add_argument(
        "--no_prefilter",
        action="store_true",
        help="",
    )

    ap.add_argument(
        "--cpus_for_alignment", "-t",
        type=int, default=1,
        help="",
    )

    ap.add_argument(
        "--version", "-v", action="version", version="%(prog)s " + __version__
    )
    ap.add_argument("--debug", action="store_true")

    return ap.parse_args(args)
