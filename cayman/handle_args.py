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


def validate_args(args: argparse.Namespace, logger: logging.Logger) -> argparse.Namespace:

    has_fastq = any(
        map(
            lambda x: x is not None,
            (
                getattr(args, arg) for arg in ("reads1", "reads2", "singles", "orphans") if hasattr(args, arg)
            )
        )
    )

    if tuple(map(bool, (has_fastq, args.bam, args.sam))).count(True) != 1:
        raise ValueError(f"Need exactly one type of input: bam={bool(args.bam)} sam={bool(args.sam)} fastq={bool(has_fastq)}.")

    args.input_type = "fastq" if has_fastq else ("bam" if args.bam else "sam")

    if args.input_type != "fastq" and args.keep_alignment_file:
        raise ValueError("Use the '--keep_alignment_file' flag only with fastQ input if you wish to store the alignment made to a SAM file.")

    if args.db_format is not None:
        logger.warning("Argument --db_format is deprecated and will be removed in a future version. Database format is now automatically detected.")

    return args


def build_parser() -> argparse.ArgumentParser:
    """ docstring """

    ap = argparse.ArgumentParser(
        prog=__toolname__,
        formatter_class=argparse.RawTextHelpFormatter,
        # parents=(log_ap,),
    )

    # ----------------------------- general arguments ----------------------------------
    ap.add_argument(
        "-l", "--log_level",
        type=int, choices=range(1, 5), default=logging.INFO
    )
    ap.add_argument(
        "-v", "--version", action="version", version="%(prog)s " + __version__
    )
    ap.add_argument(
        "--debug", action="store_true"
    )

    subparsers = ap.add_subparsers(dest="command", required=True)

    # ---------------------------- annotate proteome -----------------------------------
    annotate_proteome_ap = subparsers.add_parser(
        name = "annotate_proteome",
        help = "Annotate proteome with CAZy domains."
    )
    annotate_proteome_ap.add_argument(
        "hmmdb", type=str, help = "path to folder containing HMMs"
    )
    annotate_proteome_ap.add_argument(
        "proteins", type=str, help = "path to protein sequences in fasta format"
    )    
    annotate_proteome_ap.add_argument(
        "--cutoffs", type=str, help = "path to file containing HMM-specific p-value cutoffs"
    )    
    annotate_proteome_ap.add_argument(
        "--output_file", "-o", type=str, default="cazy_annotations.csv"
    )
    annotate_proteome_ap.add_argument(
        "--seed", "-s", type=int, default=42,
    )
    annotate_proteome_ap.add_argument(
        "--threads", "-t", type=int, default=1, help="Number of CPU threads for HMM annotation"
    )
    # Define entry point
    annotate_proteome_ap.set_defaults(func=run_proteome_annotation)

    # ------------------------------ profile WGS reads ---------------------------------
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

    # this is called --reference in gqlib and gffquant
    profile_ap.add_argument(
        "--bwa_index",
        type=str,
        required=False,
        help="Path to the bwa reference index.",
    )
    
    # TODO remove this flag and the validate warning
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
            "--bam",
            type=str,
            # nargs="+",
            help=textwrap.dedent(
                """\
                Path to a name-sorted BAM file. Ambiguous alignments need to be flagged as secondary
                alignments with the same read id as their primary alignment.
                (e.g. output from BWA mem -a). All alignments of an ambiguous group need to have MAPQ=0.
                Input from STDIN can be specified with '-'."""
            ),
        )

    profile_ap.add_argument(
        "--sam",
        type=str,
        # nargs="+",
        help=textwrap.dedent(
            """\
            Path to a name-sorted SAM file. Ambiguous alignments need to be flagged as secondary
            alignments with the same read id as their primary alignment.
            (e.g. output from BWA mem -a). All alignments of an ambiguous group need to have MAPQ=0.
            Input from STDIN can be specified with '-'."""
        ),
    )

    # only relevant with alignment input
    ap.add_argument(
        "--import_readcounts",
        type=int,
        help="Import externally derived readcounts to allow readcount-based normalisation for prefiltered bam files.",
    )

    # only relevant with alignment input; orphan reads will not have flag 0x1 set
    ap.add_argument(
        "--unmarked_orphans",
        action="store_true",
        help="Ensure that alignments from unmarked orphan reads (from preprocessing) are properly accounted for.",
    )

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
        "--keep_alignment_file",
        type=str,
        help="Save alignments in sam format to file."
    )
    profile_ap.add_argument(
        "--threads", "-t", type=int, default=1, help="Number of CPU threads for BWA MEM alignment"
    )

    profile_ap.set_defaults(func=run_profile)

    return ap
