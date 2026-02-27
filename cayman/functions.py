# pylint: disable=C0103,C0301,C0116

""" module docstring """

import logging
import os
import pathlib
import errno

# pylint: disable=W0611
from gqlib.db.db_import import SmallDatabaseImporter
from gqlib.profilers import RegionQuantifier
from gqlib.runners.alignment_runner import BwaMemRunner
from gqlib.ui.validation import check_bwa_index, check_input_reads

from .annotate.crazy_annotator import CazyAnnotator


logger = logging.getLogger(__name__)

def run_profile(args):

    if args.db_format is not None:
        logger.warning("Argument --db_format is deprecated and will be removed in a future version. Database format is now automatically detected.")

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

    db_format = None
    with open(args.annotation_db) as db_in:
        for line in db_in:
            line = line.strip()
            if line and line[0] != "#":
                if line.find(",") != -1:
                    db_format = "hmmer"
                    break
                if line.find("\t") != -1:
                    db_format = "bed"
                    break
    if db_format is None:
        logger.error("Cannot determine database format in %s.", args.annotation_db)
        raise ValueError(f"Cannot determine database format in {args.annotation_db}.")
    
    logger.info("Identified database format as `%s`.", db_format)

    db_importer = SmallDatabaseImporter(
        logger, args.annotation_db, single_category="cazy", db_format=db_format,
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
                logger.error("Failed to align. This could have different reasons:")
                logger.error(f"* Is `{args.aligner}` installed and on the path? Type `bwa mem` and see what happens.")
                logger.error("* Syntax errors or missing files. Please try running the aligner call below manually to troubleshoot the problem.")
                logger.error("* Alignment stream was interrupted, perhaps due to a memory issue.")
                
                logger.error("Aligner call was:")
                logger.error("%s", call)
                
                return 1

            logger.error("Encountered problems digesting the alignment stream:")
            logger.error("%s", err)
            logger.error("Aligner call was:")
            logger.error("%s", call)
            logger.error("Shutting down.")
            
            return 1

    profiler.finalise(restrict_reports=("raw", "rpkm",))

    return 0


def run_proteome_annotation(args):

    if args.cutoffs is None:
        args.cutoffs = os.path.join(args.hmmdb, "cutoffs.csv")

    annotator = CazyAnnotator()
    logger.info("Reading HMMs")
    if pathlib.Path(args.hmmdb).is_dir():
        annotator.read_hmms(os.path.join(args.hmmdb, "hmms"))
    elif pathlib.Path(args.hmmdb).is_file():
        annotator.load_hmm(args.hmmdb)
    else:
        logger.error(f"File {str(args.hmmdb)} does not exist")
        return  errno.ENOENT
    logger.info("Reading sequences")
    annotator.read_sequences(path_to_sequences=args.proteins)
    logger.info("Annotating sequences (can take a few minutes; be patient)")
    annotator.annotate_sequences_with_all_hmms(threads=args.threads)
    logger.info("Filtering and merging annotations over folds")
    annotator.curate_annotations(precomputed_hmm_cutoffs=args.cutoffs)
    annotator.annotations_filtered.to_csv(args.output_file, index=False)

    return 0
