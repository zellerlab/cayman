# pylint: disable=C0301,C0103
""" module docstring """

import argparse
import logging
import textwrap

from . import __version__
from . import __toolname__
from .functions import run_profile, run_proteome_annotation


def set_log_lvl_from_args(argv):
    log_ap = argparse.ArgumentParser(prog=__toolname__, add_help=False)
    log_ap.add_argument(
        "-l", "--log_level",
        type=int, choices=range(1, 5), default=logging.INFO
    )
    log_args, _ = log_ap.parse_known_args(argv)

    try:
        logging.basicConfig(
            level=log_args.log_level,
            format='[%(asctime)s] %(message)s'
        )
    except ValueError as invalid_loglevel_err:
        raise ValueError(
            f"Invalid log level: {log_args.log_level}"
        ) from invalid_loglevel_err


def build_parser() -> argparse.ArgumentParser:
    """ docstring """

    ap = argparse.ArgumentParser(
        prog=__toolname__,
        formatter_class=argparse.RawTextHelpFormatter,
        # parents=(log_ap,),
    )

    ap.add_argument(
        "-l", "--log_level",
        type=int, choices=range(1, 5), default=logging.INFO
    )

    ap.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )
	
    ap.add_argument("--debug", action="store_true")

    subparsers = ap.add_subparsers(dest="command", required=True)

    annotate_proteome_ap = subparsers.add_parser("annotate_proteome", help="Annotate proteome with CAZy domains.")
    annotate_proteome_ap.add_argument("hmmdb", type=str, help = "path to folder containing HMMs")
    annotate_proteome_ap.add_argument("proteins", type=str, help = "path to protein sequences in fasta format")    
    annotate_proteome_ap.add_argument("--cutoffs", type=str, help = "path to file containing HMM-specific p-value cutoffs")    

    annotate_proteome_ap.add_argument("--output_file", "-o", type=str, default="cazy_annotations.csv")
    annotate_proteome_ap.add_argument("--threads", "-t", type=int, default=1)
    

    annotate_proteome_ap.set_defaults(func=run_proteome_annotation)

    profile_ap = subparsers.add_parser("profile", help="Profile WGS reads.")
    
    profile_ap.add_argument(
        "annotation_db",
        type=str,
        help=textwrap.dedent(
            """\
            Path to a text file containing the domain annotation. This needs to be a 4-column file such as bed4.
            """
        ),
    )

    profile_ap.add_argument(
        "bwa_index",
        type=str,
        help=textwrap.dedent(
            """\
            Path to the bwa reference index.
            """
        ),
    )

    profile_ap.add_argument(
        "--db_format",
        type=str,
        # default="hmmer",
        choices=("bed", "hmmer"),
        help="Format of the annotation database.",
    )

    # ap.add_argument(
    #     "--db_coordinates",
    #     type=str,
    #     default="bed",
    #     choices=("bed", "hmmer"),
    #     help="Coordinate format for text-based annotation databases. bed=[start, end), hmmer=[start, end]"
    # )

    # ap.add_argument(
    #     "--db_separator",
    #     type=str,
    #     default="\t",
    #     help="Separator-character for the annotation database file Default: '\\t'."
    # )


    profile_ap.add_argument(
        "-1",
        dest="reads1",
        nargs="*",
        type=str,
        help="A forward/R1 read fastq file. Multiple files can be separated by spaces."
    )

    profile_ap.add_argument(
        "-2",
        dest="reads2",
        nargs="*",
        type=str,
        help="A reverse/R2 read fastq file. Multiple files can be separated by spaces."
    )

    profile_ap.add_argument(
        "--singles", "-s",
        nargs="*",
        type=str,
        help="A single-end library read fastq file. Multiple files can be separated by spaces." 
    )

    profile_ap.add_argument(
        "--orphans",
        nargs="*",
        type=str,
        help="An orphan read fastq file. Multiple files can be separated by spaces."
    )

    profile_ap.add_argument(
        "--out_prefix",
        "-o",
        type=str,
        default=__toolname__,
        help="Prefix for output files.",
    )

    profile_ap.add_argument(
        "--min_identity",
        type=float,
        default=0.97,
        help="Minimum sequence identity [n_match/length] for an alignment to be considered.",
    )

    profile_ap.add_argument(
        "--min_seqlen",
        type=int,
        default=45,
        help="Minimum read length [bp] for an alignment to be considered.",
    )

    profile_ap.add_argument(
        "--cpus_for_alignment", "-t",
        type=int, default=1,
        help="",
    )

    profile_ap.set_defaults(func=run_profile)

    # ap.add_argument(
    #     "--version", "-v", action="version", version="%(prog)s " + __version__
    # )
    # ap.add_argument("--debug", action="store_true")

    return ap
